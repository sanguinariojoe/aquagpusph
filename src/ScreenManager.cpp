/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Terminal output, with Log automatic copying.
 * (See Aqua::InputOutput::ScreenManager for details)
 */

#include <sys/time.h>
#include <unistd.h>

#include <ScreenManager.h>
#include <FileManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>

#ifdef HAVE_NCURSES
    WINDOW *wnd, *log_wnd;
#endif

namespace Aqua{ namespace InputOutput{

ScreenManager::ScreenManager()
    : _last_row(0)
{
    int i;

    #ifdef HAVE_NCURSES
        wnd = initscr();
        if(!wnd){
            addMessageF(L_INFO, "Failure initializating the screen manager\n");
            return;
        }
        if(has_colors())
            start_color();
        log_wnd = NULL;
    #endif

    gettimeofday(&_start_time, NULL);
}

ScreenManager::~ScreenManager()
{
    #ifdef HAVE_NCURSES
        // Stop ncurses
        endwin();
    #endif
}

void ScreenManager::initFrame()
{
    #ifdef HAVE_NCURSES
        // Clear the entire frame
        if(!wnd)
            return;
        clear();

        // Init the colour pairs
        init_pair(1, COLOR_WHITE, COLOR_BLACK);
        init_pair(2, COLOR_GREEN, COLOR_BLACK);
        init_pair(3, COLOR_BLUE, COLOR_BLACK);
        init_pair(4, COLOR_YELLOW, COLOR_BLACK);
        init_pair(5, COLOR_RED, COLOR_BLACK);
        init_pair(6, COLOR_MAGENTA, COLOR_BLACK);
        init_pair(7, COLOR_CYAN, COLOR_BLACK);

        // Setup the default font
        attrset(A_NORMAL);
        attron(COLOR_PAIR(1));

        // Set the cursor at the start of the window
        _last_row = 0;
    #endif
}

void ScreenManager::endFrame()
{
    #ifdef HAVE_NCURSES
        printLog();
        refreshAll();
    #else
        std::cout << std::flush;
    #endif
}

void ScreenManager::writeReport(std::string input,
                                std::string color,
                                bool bold)
{
    if(!input.size()){
        return;
    }
    #ifndef HAVE_NCURSES
        std::cout << input;
        if(!hasSuffix(input, "\n")){
            std::cout << endl;
        }
    #else
        std::string msg = rtrim(input);
        size_t start = 0, end = 0;
        while((end = str.find("\n", start)) != std::string::npos) {
            writeReport(str.substr(start, end), color, bold);
            start = end;
        }
        if(str.find("\n", start)) != std::string::npos) {
            // Alreadey processed by pieces
            return;
        }

        // Replace the tabulators by spaces
        replaceAll(msg, "\t", " ");

        // Check if the message is larger than the terminal output
        int rows, cols;
        getmaxyx(wnd, rows, cols);
        if(msg.size() > cols){
            // We can try to split it by a blank space, and if it fails just let
            // ncurses select how to divide the string
            size_t last = 0;
            end = 0;
            while((end = str.find(" ", end)) != std::string::npos) {
                if (end > cols)
                    break;
                last = end;
            }
            if (last) {
                writeReport(str.substr(0, last), color, bold);
                writeReport(str.substr(last), color, bold);
                return;
            }
        }

        // Select the font
        attrset(A_NORMAL);
        unsigned int pair_id = 1;
        if(!strcmp(color, "white")){
            pair_id = 1;
        }
        else if(!strcmp(color, "green")){
            pair_id = 2;
        }
        else if(!strcmp(color, "blue")){
            pair_id = 3;
        }
        else if(!strcmp(color, "yellow")){
            pair_id = 4;
        }
        else if(!strcmp(color, "red")){
            pair_id = 5;
        }
        else if(!strcmp(color, "magenta")){
            pair_id = 6;
        }
        else if(!strcmp(color, "cyan")){
            pair_id = 7;
        }
        else{
            std::ostringstream err_msg;
            err_msg << "Invalid message color \"" << color << "\"" << std::endl;
            addMessageF(L_ERROR, err_msg.str());
        }
        attron(COLOR_PAIR(pair_id));
        if(bold){
            attron(A_BOLD);
        }

        // Print the message
        int row = _last_row;
        move(row, 0);
        printw(msg.c_str());
        // The message may require several lines
        _last_row += msg.size() / cols + 1;

        // Refresh
        // printLog();
        // refresh();
    #endif
}

void ScreenManager::addMessage(TLogLevel level, std::string log, const char *func)
{
    std::ostringstream fname;
    fname << "(" << func << "): ";

    // Send the info to the log file (if possible)
    if(FileManager::singleton()){
        FILE *LogFileID = FileManager::singleton()->logFile();
        if(LogFileID){
            if(level == L_INFO)
                fprintf(LogFileID,
                        "<b><font color=\"#000000\">[INFO] %s%s</font></b><br>",
                        fname.str().c_str(),
                        log.c_str());
            else if(level == L_WARNING)
                fprintf(LogFileID,
                        "<b><font color=\"#ff9900\">[WARNING] %s%s</font></b><br>",
                        fname.str().c_str(),
                        log.c_str());
            else if(level == L_ERROR)
                fprintf(LogFileID,
                        "<b><font color=\"#dd0000\">[ERROR] %s%s</font></b><br>",
                        fname.str().c_str(),
                        log.c_str());
            else{
                fprintf(LogFileID,
                        "<font color=\"#000000\">%s%s</font>",
                        fname.str().c_str(),
                        log.c_str());
                if(hasSuffix(log, "\n"))
                    fprintf(LogFileID, "<br>");
            }
            fflush(LogFileID);
        }
    }
    // Compatibility mode for destroyed ScreenManager situations
    if(!ScreenManager::singleton()){
        if(level == L_INFO)
            cout << "INFO ";
        else if(level == L_WARNING)
            cout << "WARNING ";
        else if(level == L_ERROR)
            cout << "ERROR ";
        cout << fname.str() << log << std::flush;
        return;
    }

    // Add the new message to the log register
    std::ostringstream msg;
    msg << fname.str() << log;
    _log_level.insert(_log_level.begin(), level);
    _log.insert(_log.begin(), msg.str());

    // Filter out the messages that never would be printed because the window
    // has not space enough (1 if ncurses is not activated)
    int rows=1, cols;
    #ifdef HAVE_NCURSES
        getmaxyx(wnd, rows, cols);
    #endif
    while(_log_level.size() > (unsigned int)rows){
        _log_level.pop_back();
    }
    while(_log.size() > (unsigned int)rows){
        _log.pop_back();
    }

    printLog();
    refreshAll();
}

void ScreenManager::printDate(TLogLevel level)
{
    char msg[512];
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    sprintf(msg, "%s\n", ctime(&seconds));
    addMessage(level, msg);
}

void ScreenManager::printOpenCLError(cl_int error, TLogLevel level)
{
    switch(error) {
        case CL_DEVICE_NOT_FOUND:
            addMessage(level, "\tCL_DEVICE_NOT_FOUND\n"); break;
        case CL_DEVICE_NOT_AVAILABLE:
            addMessage(level, "\tCL_DEVICE_NOT_AVAILABLE\n"); break;
        case CL_COMPILER_NOT_AVAILABLE:
            addMessage(level, "\tCL_COMPILER_NOT_AVAILABLE\n"); break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            addMessage(level, "\tCL_MEM_OBJECT_ALLOCATION_FAILURE\n"); break;
        case CL_OUT_OF_RESOURCES:
            addMessage(level, "\tCL_OUT_OF_RESOURCES\n"); break;
        case CL_OUT_OF_HOST_MEMORY:
            addMessage(level, "\tCL_OUT_OF_HOST_MEMORY\n"); break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            addMessage(level, "\tCL_PROFILING_INFO_NOT_AVAILABLE\n"); break;
        case CL_MEM_COPY_OVERLAP:
            addMessage(level, "\tCL_MEM_COPY_OVERLAP\n"); break;
        case CL_IMAGE_FORMAT_MISMATCH:
            addMessage(level, "\tCL_IMAGE_FORMAT_MISMATCH\n"); break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            addMessage(level, "\tCL_IMAGE_FORMAT_NOT_SUPPORTED\n"); break;
        case CL_BUILD_PROGRAM_FAILURE:
            addMessage(level, "\tCL_BUILD_PROGRAM_FAILURE\n"); break;
        case CL_MAP_FAILURE:
            addMessage(level, "\tCL_MAP_FAILURE\n"); break;
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            addMessage(level, "\tCL_MISALIGNED_SUB_BUFFER_OFFSET\n"); break;
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            addMessage(level, "\tCL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST\n"); break;
        case CL_INVALID_VALUE:
            addMessage(level, "\tCL_INVALID_VALUE\n"); break;
        case CL_INVALID_DEVICE_TYPE:
            addMessage(level, "\tCL_INVALID_DEVICE_TYPE\n"); break;
        case CL_INVALID_PLATFORM:
            addMessage(level, "\tCL_INVALID_PLATFORM\n"); break;
        case CL_INVALID_DEVICE:
            addMessage(level, "\tCL_INVALID_DEVICE\n"); break;
        case CL_INVALID_CONTEXT:
            addMessage(level, "\tCL_INVALID_CONTEXT\n"); break;
        case CL_INVALID_QUEUE_PROPERTIES:
            addMessage(level, "\tCL_INVALID_QUEUE_PROPERTIES\n"); break;
        case CL_INVALID_COMMAND_QUEUE:
            addMessage(level, "\tCL_INVALID_COMMAND_QUEUE\n"); break;
        case CL_INVALID_HOST_PTR:
            addMessage(level, "\tCL_INVALID_HOST_PTR\n"); break;
        case CL_INVALID_MEM_OBJECT:
            addMessage(level, "\tCL_INVALID_MEM_OBJECT\n"); break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            addMessage(level, "\tCL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n"); break;
        case CL_INVALID_IMAGE_SIZE:
            addMessage(level, "\tCL_INVALID_IMAGE_SIZE\n"); break;
        case CL_INVALID_SAMPLER:
            addMessage(level, "\tCL_INVALID_SAMPLER\n"); break;
        case CL_INVALID_BINARY:
            addMessage(level, "\tCL_INVALID_BINARY\n"); break;
        case CL_INVALID_BUILD_OPTIONS:
            addMessage(level, "\tCL_INVALID_BUILD_OPTIONS\n"); break;
        case CL_INVALID_PROGRAM:
            addMessage(level, "\tCL_INVALID_PROGRAM\n"); break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            addMessage(level, "\tCL_INVALID_PROGRAM_EXECUTABLE\n"); break;
        case CL_INVALID_KERNEL_NAME:
            addMessage(level, "\tCL_INVALID_KERNEL_NAME\n"); break;
        case CL_INVALID_KERNEL_DEFINITION:
            addMessage(level, "\tCL_INVALID_KERNEL_DEFINITION\n"); break;
        case CL_INVALID_KERNEL:
            addMessage(level, "\tCL_INVALID_KERNEL\n"); break;
        case CL_INVALID_ARG_INDEX:
            addMessage(level, "\tCL_INVALID_ARG_INDEX\n"); break;
        case CL_INVALID_ARG_VALUE:
            addMessage(level, "\tCL_INVALID_ARG_VALUE\n"); break;
        case CL_INVALID_ARG_SIZE:
            addMessage(level, "\tCL_INVALID_ARG_SIZE\n"); break;
        case CL_INVALID_KERNEL_ARGS:
            addMessage(level, "\tCL_INVALID_KERNEL_ARGS\n"); break;
        case CL_INVALID_WORK_DIMENSION:
            addMessage(level, "\tCL_INVALID_WORK_DIMENSION\n"); break;
        case CL_INVALID_WORK_GROUP_SIZE:
            addMessage(level, "\tCL_INVALID_WORK_GROUP_SIZE\n"); break;
        case CL_INVALID_WORK_ITEM_SIZE:
            addMessage(level, "\tCL_INVALID_WORK_ITEM_SIZE\n"); break;
        case CL_INVALID_GLOBAL_OFFSET:
            addMessage(level, "\tCL_INVALID_GLOBAL_OFFSET\n"); break;
        case CL_INVALID_EVENT_WAIT_LIST:
            addMessage(level, "\tCL_INVALID_EVENT_WAIT_LIST\n"); break;
        case CL_INVALID_EVENT:
            addMessage(level, "\tCL_INVALID_EVENT\n"); break;
        case CL_INVALID_OPERATION:
            addMessage(level, "\tCL_INVALID_OPERATION\n"); break;
        case CL_INVALID_GL_OBJECT:
            addMessage(level, "\tCL_INVALID_GL_OBJECT\n"); break;
        case CL_INVALID_BUFFER_SIZE:
            addMessage(level, "\tCL_INVALID_BUFFER_SIZE\n"); break;
        case CL_INVALID_MIP_LEVEL:
            addMessage(level, "\tCL_INVALID_MIP_LEVEL\n"); break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
            addMessage(level, "\tCL_INVALID_GLOBAL_WORK_SIZE\n"); break;
        case CL_INVALID_PROPERTY:
            addMessage(level, "\tCL_INVALID_PROPERTY\n"); break;
        default:
            addMessage(level, "\tUnhandled exception\n"); break;
    }
}

void ScreenManager::printLog()
{
    unsigned int i;
    #ifndef HAVE_NCURSES
    for(i = 0; i < _log_level.size(); i++){
        if(_log_level.at(i) == L_INFO)
            cout << "INFO ";
        else if(_log_level.at(i) == L_WARNING)
            cout << "WARNING ";
        else if(_log_level.at(i) == L_ERROR)
            cout << "ERROR ";
        cout << _log.at(i).c_str() << std::flush;
    }
    #else
        // Get a good position candidate
        int row = _last_row + 2, rows, cols;
        getmaxyx(wnd, rows, cols);
        if(row > rows - 2){
            row = rows - 2;
        }

        // Setup the subwindow
        if(log_wnd){
            wclear(log_wnd);
            delwin(log_wnd);
            log_wnd = NULL;
        }
        log_wnd = subwin(wnd, rows - row, cols, row, 0);
        wclear(log_wnd);

        // Print the info
        wattrset(log_wnd, A_NORMAL);
        wattron(log_wnd, A_BOLD);
        wattron(log_wnd, COLOR_PAIR(1));
        wmove(log_wnd, 0, 0);
        wprintw(log_wnd, "--- Log registry ------------------------------");
        unsigned int lines = 1;
        for(i = 0; i < _log_level.size(); i++){
            wattrset(log_wnd, A_NORMAL);
            wattron(log_wnd, COLOR_PAIR(1));
            if(_log_level.at(i) == L_INFO){
                wattron(log_wnd, A_NORMAL);
                wattron(log_wnd, COLOR_PAIR(1));
            }
            else if(_log_level.at(i) == L_WARNING){
                wattron(log_wnd, A_BOLD);
                wattron(log_wnd, COLOR_PAIR(4));
            }
            else if(_log_level.at(i) == L_ERROR){
                wattron(log_wnd, A_BOLD);
                wattron(log_wnd, COLOR_PAIR(5));
            }
            else{
            }
            wmove(log_wnd, lines, 0);
            wprintw(log_wnd, _log.at(i).c_str());
            lines += _log.at(i).size() / cols + 1;
        }
        // wrefresh(log_wnd);
    #endif
}

void ScreenManager::refreshAll()
{
    #ifdef HAVE_NCURSES
        refresh();
        wrefresh(log_wnd);
    #endif
}

}}  // namespace
