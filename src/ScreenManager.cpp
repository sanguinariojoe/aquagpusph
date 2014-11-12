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

#include <CL/cl.h>

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
            addMessageF(1, "Failure initializating the screen manager\n");
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
    unsigned int i;
    #ifdef HAVE_NCURSES
        // Stop ncurses
        endwin();
    #endif
    // Free allocated memory
    for(i = 0; i < _log.size(); i++){
        delete[] _log.at(i);
    }
    _log_level.clear();
    _log.clear();
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

        // Setup the default font
        attron(A_NORMAL);
        attron(COLOR_PAIR(1));
    #endif
}

void ScreenManager::endFrame()
{
    #ifdef HAVE_NCURSES
        printLog();
        refresh();
    #else
        fflush(stdout);
    #endif
}

void ScreenManager::addMessage(int level, const char *log, const char *func)
{
    char fname[256]; strcpy(fname, "");
    if(func){
        sprintf(fname, "(%s): ", func);
    }
    // Send the info to the log file (if possible)
    if(FileManager::singleton()){
        FILE *LogFileID = FileManager::singleton()->logFile();
        if(LogFileID){
            if(level == 1)        // Log
                fprintf(LogFileID, "<b><font color=\"#000000\">[INFO] %s%s</font></b><br>", fname, log);
            else if(level == 2)        // Warning
                fprintf(LogFileID, "<b><font color=\"#ff9900\">[WARNING] %s%s</font></b><br>", fname, log);
            else if(level == 3)        // Error
                fprintf(LogFileID, "<b><font color=\"#dd0000\">[ERROR] %s%s</font></b><br>", fname, log);
            else{
                fprintf(LogFileID, "<font color=\"#000000\">%s%s</font>", fname, log);
                if(log[strlen(log)-1] == '\n')
                fprintf(LogFileID, "<br>");
            }
            fflush(LogFileID);
        }
    }
    // Compatibility mode for destroyed ScreenManager situations
    if(!ScreenManager::singleton()){
        if(level == 1)
            printf("INFO ");
        else if(level == 2)
            printf("WARNING ");
        else if(level == 3)
            printf("ERROR ");
        printf("%s%s", fname, log);
        fflush(stdout);
        return;
    }

    // Add the new message to the log register
    char *msg = new char[strlen(fname) + strlen(log) + 1];
    strcpy(msg, fname);
    strcat(msg, log);
    _log_level.insert(_log_level.begin(), level);
    _log.insert(_log.begin(), msg);

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
}

void ScreenManager::printDate(int level)
{
    char msg[512];
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    sprintf(msg, "%s\n", ctime(&seconds));
    addMessage(level, msg);
}

void ScreenManager::printOpenCLError(int error, int level)
{
    char msg[128];
    strcpy(msg, "\tUnhandled exception\n");
    if(error == CL_DEVICE_NOT_FOUND)
        strcpy(msg, "\tCL_DEVICE_NOT_FOUND\n");
    else if(error == CL_DEVICE_NOT_AVAILABLE)
        strcpy(msg, "\tCL_DEVICE_NOT_AVAILABLE\n");
    else if(error == CL_COMPILER_NOT_AVAILABLE)
        strcpy(msg, "\tCL_COMPILER_NOT_AVAILABLE\n");
    else if(error == CL_MEM_OBJECT_ALLOCATION_FAILURE)
        strcpy(msg, "\tCL_MEM_OBJECT_ALLOCATION_FAILURE\n");
    else if(error == CL_OUT_OF_RESOURCES)
        strcpy(msg, "\tCL_OUT_OF_RESOURCES\n");
    else if(error == CL_OUT_OF_HOST_MEMORY)
        strcpy(msg, "\tCL_OUT_OF_HOST_MEMORY\n");
    else if(error == CL_PROFILING_INFO_NOT_AVAILABLE)
        strcpy(msg, "\tCL_PROFILING_INFO_NOT_AVAILABLE\n");
    else if(error == CL_MEM_COPY_OVERLAP)
        strcpy(msg, "\tCL_MEM_COPY_OVERLAP\n");
    else if(error == CL_IMAGE_FORMAT_MISMATCH)
        strcpy(msg, "\tCL_IMAGE_FORMAT_MISMATCH\n");
    else if(error == CL_IMAGE_FORMAT_NOT_SUPPORTED)
        strcpy(msg, "\tCL_IMAGE_FORMAT_NOT_SUPPORTED\n");
    else if(error == CL_BUILD_PROGRAM_FAILURE)
        strcpy(msg, "\tCL_BUILD_PROGRAM_FAILURE\n");
    else if(error == CL_MAP_FAILURE)
        strcpy(msg, "\tCL_MAP_FAILURE\n");
    else if(error == CL_MISALIGNED_SUB_BUFFER_OFFSET)
        strcpy(msg, "\tCL_MISALIGNED_SUB_BUFFER_OFFSET\n");
    else if(error == CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
        strcpy(msg, "\tCL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST\n");
    else if(error == CL_COMPILE_PROGRAM_FAILURE)
        strcpy(msg, "\tCL_COMPILE_PROGRAM_FAILURE\n");
    else if(error == CL_LINKER_NOT_AVAILABLE)
        strcpy(msg, "\tCL_LINKER_NOT_AVAILABLE\n");
    else if(error == CL_LINK_PROGRAM_FAILURE)
        strcpy(msg, "\tCL_LINK_PROGRAM_FAILURE\n");
    else if(error == CL_DEVICE_PARTITION_FAILED)
        strcpy(msg, "\tCL_DEVICE_PARTITION_FAILED\n");
    else if(error == CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
        strcpy(msg, "\tCL_KERNEL_ARG_INFO_NOT_AVAILABLE\n");
    else if(error == CL_INVALID_VALUE)
        strcpy(msg, "\tCL_INVALID_VALUE\n");
    else if(error == CL_INVALID_DEVICE_TYPE)
        strcpy(msg, "\tCL_INVALID_DEVICE_TYPE\n");
    else if(error == CL_INVALID_PLATFORM)
        strcpy(msg, "\tCL_INVALID_PLATFORM\n");
    else if(error == CL_INVALID_DEVICE)
        strcpy(msg, "\tCL_INVALID_DEVICE\n");
    else if(error == CL_INVALID_CONTEXT)
        strcpy(msg, "\tCL_INVALID_CONTEXT\n");
    else if(error == CL_INVALID_QUEUE_PROPERTIES)
        strcpy(msg, "\tCL_INVALID_QUEUE_PROPERTIES\n");
    else if(error == CL_INVALID_COMMAND_QUEUE)
        strcpy(msg, "\tCL_INVALID_COMMAND_QUEUE\n");
    else if(error == CL_INVALID_HOST_PTR)
        strcpy(msg, "\tCL_INVALID_HOST_PTR\n");
    else if(error == CL_INVALID_MEM_OBJECT)
        strcpy(msg, "\tCL_INVALID_MEM_OBJECT\n");
    else if(error == CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
        strcpy(msg, "\tCL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
    else if(error == CL_INVALID_IMAGE_SIZE)
        strcpy(msg, "\tCL_INVALID_IMAGE_SIZE\n");
    else if(error == CL_INVALID_SAMPLER)
        strcpy(msg, "\tCL_INVALID_SAMPLER\n");
    else if(error == CL_INVALID_BINARY)
        strcpy(msg, "\tCL_INVALID_BINARY\n");
    else if(error == CL_INVALID_BUILD_OPTIONS)
        strcpy(msg, "\tCL_INVALID_BUILD_OPTIONS\n");
    else if(error == CL_INVALID_PROGRAM)
        strcpy(msg, "\tCL_INVALID_PROGRAM\n");
    else if(error == CL_INVALID_PROGRAM_EXECUTABLE)
        strcpy(msg, "\tCL_INVALID_PROGRAM_EXECUTABLE\n");
    else if(error == CL_INVALID_KERNEL_NAME)
        strcpy(msg, "\tCL_INVALID_KERNEL_NAME\n");
    else if(error == CL_INVALID_KERNEL_DEFINITION)
        strcpy(msg, "\tCL_INVALID_KERNEL_DEFINITION\n");
    else if(error == CL_INVALID_KERNEL)
        strcpy(msg, "\tCL_INVALID_KERNEL\n");
    else if(error == CL_INVALID_ARG_INDEX)
        strcpy(msg, "\tCL_INVALID_ARG_INDEX\n");
    else if(error == CL_INVALID_ARG_VALUE)
        strcpy(msg, "\tCL_INVALID_ARG_VALUE\n");
    else if(error == CL_INVALID_ARG_SIZE)
        strcpy(msg, "\tCL_INVALID_ARG_SIZE\n");
    else if(error == CL_INVALID_KERNEL_ARGS)
        strcpy(msg, "\tCL_INVALID_KERNEL_ARGS\n");
    else if(error == CL_INVALID_WORK_DIMENSION)
        strcpy(msg, "\tCL_INVALID_WORK_DIMENSION\n");
    else if(error == CL_INVALID_WORK_GROUP_SIZE)
        strcpy(msg, "\tCL_INVALID_WORK_GROUP_SIZE\n");
    else if(error == CL_INVALID_WORK_ITEM_SIZE)
        strcpy(msg, "\tCL_INVALID_WORK_ITEM_SIZE\n");
    else if(error == CL_INVALID_GLOBAL_OFFSET)
        strcpy(msg, "\tCL_INVALID_GLOBAL_OFFSET\n");
    else if(error == CL_INVALID_EVENT_WAIT_LIST)
        strcpy(msg, "\tCL_INVALID_EVENT_WAIT_LIST\n");
    else if(error == CL_INVALID_EVENT)
        strcpy(msg, "\tCL_INVALID_EVENT\n");
    else if(error == CL_INVALID_OPERATION)
        strcpy(msg, "\tCL_INVALID_OPERATION\n");
    else if(error == CL_INVALID_GL_OBJECT)
        strcpy(msg, "\tCL_INVALID_GL_OBJECT\n");
    else if(error == CL_INVALID_BUFFER_SIZE)
        strcpy(msg, "\tCL_INVALID_BUFFER_SIZE\n");
    else if(error == CL_INVALID_MIP_LEVEL)
        strcpy(msg, "\tCL_INVALID_MIP_LEVEL\n");
    else if(error == CL_INVALID_GLOBAL_WORK_SIZE)
        strcpy(msg, "\tCL_INVALID_GLOBAL_WORK_SIZE\n");
    else if(error == CL_INVALID_PROPERTY)
        strcpy(msg, "\tCL_INVALID_PROPERTY\n");
    else if(error == CL_INVALID_IMAGE_DESCRIPTOR)
        strcpy(msg, "\tCL_INVALID_IMAGE_DESCRIPTOR\n");
    else if(error == CL_INVALID_COMPILER_OPTIONS)
        strcpy(msg, "\tCL_INVALID_COMPILER_OPTIONS\n");
    else if(error == CL_INVALID_LINKER_OPTIONS)
        strcpy(msg, "\tCL_INVALID_LINKER_OPTIONS\n");
    else if(error == CL_INVALID_DEVICE_PARTITION_COUNT)
        strcpy(msg, "\tCL_INVALID_DEVICE_PARTITION_COUNT\n");
    addMessage(level, msg);
}

void ScreenManager::printLog()
{
    unsigned int i;
    #ifndef HAVE_NCURSES
    for(i = 0; i < _log_level.size(); i++){
        if(_log_level.at(i) == 1)
            printf("INFO %s", _log.at(i));
        else if(_log_level.at(i) == 2)
            printf("WARNING %s", _log.at(i));
        else if(_log_level.at(i) == 3)
            printf("ERROR %s", _log.at(i));
        else
            printf("%s", _log.at(i));
        fflush(stdout);
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
            delwin(log_wnd);
            log_wnd = NULL;
        }
        log_wnd = subwin(wnd, rows - row, cols, row, 0);
        wclear(log_wnd);

        // Print the info
        wattron(log_wnd, A_NORMAL);
        wattron(log_wnd, COLOR_PAIR(1));
        wmove(log_wnd, 0, 2);
        wprintw(log_wnd, "--- Log registry ------------------------------");
        for(i = 0; i < _log_level.size(); i++){
            if(_log_level.at(i) == 1){
                wattron(log_wnd, A_BOLD);
                wattron(log_wnd, COLOR_PAIR(1));
            }
            else if(_log_level.at(i) == 2){
                wattron(log_wnd, A_BOLD);
                wattron(log_wnd, COLOR_PAIR(4));
            }
            else if(_log_level.at(i) == 3){
                wattron(log_wnd, A_BOLD);
                wattron(log_wnd, COLOR_PAIR(5));
            }
            else{
            }
            wmove(log_wnd, i + 1, 2);
            wprintw(log_wnd, _log.at(i));
        }
        wrefresh(log_wnd);
    #endif
}

}}  // namespace
