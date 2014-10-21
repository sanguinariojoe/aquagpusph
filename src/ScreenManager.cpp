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
#include <CalcServer.h>
#include <TimeManager.h>

#ifdef HAVE_NCURSES
    WINDOW *wnd;
#endif

namespace Aqua{ namespace InputOutput{

ScreenManager::ScreenManager()
{
    int i;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    TimeManager *T = TimeManager::singleton();

    _n_log = 20;
    #ifdef HAVE_NCURSES
        wnd = initscr();
        if(!wnd){
            addMessageF(1, "Failure initializating the screen manager\n");
            return;
        }
        if(has_colors())
            start_color();
        int rows, cols;
        getmaxyx(wnd,rows,cols);
        // 17 lines for the fix text (including the log registry title)
        _n_log = rows - 17;
    #endif

    gettimeofday(&_start_time, NULL);

    _c_log = new int[_n_log];
    _m_log = new char*[_n_log];
    for(i=0;i<_n_log;i++)
    {
        _c_log[i] = 0;
        _m_log[i] = new char[128];
    }

    _old_frame = -1;
}

ScreenManager::~ScreenManager()
{
    int i;
    // Stop ncurses
    #ifdef HAVE_NCURSES
        endwin();
    #endif
    // Free allocated memory
    delete[] _c_log;
    for(i=0;i<_n_log;i++)
        delete[] _m_log[i];
    delete[] _m_log;
}

void ScreenManager::update()
{
    CalcServer::CalcServer *c = CalcServer::CalcServer::singleton();
    TimeManager *T = TimeManager::singleton();
    if(!c->verbose_level)
    {
        #ifdef HAVE_NCURSES
            move(10, 2);
            printw("Verbose level = 0, no data will be printed.");
            refresh();
        #endif
        return;
    }

    if((c->verbose_level == 1) && (_old_frame == T->frame()))
    {
        return;
    }
    _old_frame = T->frame();

    float Percentage = 0;
    int iPercentage = 0;
    long seconds;
    unsigned int i,line;

    #ifdef HAVE_NCURSES
        if(!wnd)
            return;
        clear();
        init_pair(1, COLOR_WHITE, COLOR_BLACK);
        init_pair(2, COLOR_GREEN, COLOR_BLACK);
        init_pair(3, COLOR_BLUE, COLOR_BLACK);
        init_pair(4, COLOR_YELLOW, COLOR_BLACK);
        init_pair(5, COLOR_RED, COLOR_BLACK);
    #endif

    #ifdef HAVE_NCURSES
        attron(A_NORMAL | COLOR_PAIR(2));
        move(1, 2);
        printw("Time: %g s (dt = %g s).", T->time() - T->dt(), T->dt());
    #else
        printf("Time: %g s (dt = %g s).\n", T->time() - T->dt(), T->dt());
    #endif

    float _time = T->time() - T->startTime();
    float mMaxTime = T->maxTime() - T->startTime();
    int mFrame = T->frame() - T->startFrame();
    int mMaxFrame = T->maxFrame() - T->startFrame();
    float mDivisor = 1.f;
    Percentage = 0.f;
    int TimeRules = 0;
    if(T->maxTime() >= 0.f) {
        Percentage = _time / mMaxTime;
        mDivisor = mMaxTime;
        TimeRules = 1;
    }
    if(T->maxFrame() >= 0) {
        if(!TimeRules)
            mDivisor = mMaxFrame-1;
        float fPercentage = (float)(mFrame-1) / (mMaxFrame-1);
        if(fPercentage > Percentage) {
            Percentage = fPercentage;
            mDivisor = mMaxFrame-1;
            TimeRules = 0;
        }
    }
    float subPercentage = 0.f;
    if(!TimeRules) {
        if((T->outputFPS() >= 0) && (T->time() > 0.f)){
            subPercentage = (T->time() - T->outputTime()) / (1.f/T->outputFPS());
        }
        if(T->outputIPF() >= 0) {
            float subfPercentage = (float)(T->step() - T->outputStep()) / T->outputIPF();
            if(subfPercentage > subPercentage)
                subPercentage = subfPercentage;
        }
    }
    Percentage += subPercentage / mDivisor;
    iPercentage = (int)(Percentage * 100.0);
    if(iPercentage < 0)
        iPercentage = 0;
    #ifdef HAVE_NCURSES
        move(1, 45);
        printw("[%d", iPercentage);
        move(1, 50);
        printw("]");
    #else
        printf("\t[%d]", iPercentage);
    #endif

    #ifdef HAVE_NCURSES
        move(2, 2);
        printw("Frame: %d.", T->frame());
        move(2, 18);
        printw("Step: %u.", T->step());
    #else
        printf("\tFrame: %d.", T->frame());
        printf("\tStep: %u.", T->step());
    #endif

    gettimeofday(&_actual_time, NULL);
    seconds = _actual_time.tv_sec  - _start_time.tv_sec;
    seconds = (long unsigned int)(seconds / Percentage) - seconds;
    #ifdef HAVE_NCURSES
        move(2, 35);
        printw("[ETA: %lu s]", seconds);
    #else
        printf("\t[ETA: %lu s]\n", seconds);
    #endif
    // Print profiling info if needed
    #ifdef HAVE_GPUPROFILE
        float Sum = c->predictor->profileTime() +
                    c->grid->profileTime() +
                    c->link_list->profileTime() +
                    c->rates->profileTime() +
                    c->elastic_bounce->profileTime() +
                    c->de_Leffe->profileTime() +
                    c->corrector->profileTime() +
                    c->domain->profileTime() +
                    c->time_step->profileTime() +
                    c->dens_int->profileTime() +
                    c->sensors->profileTime();
        for(i=0;i<c->motions.size();i++){
            Sum += c->motions.at(i)->profileTime();
        }
        char *str = new char[512]; strcpy(str, "");
        char *aux = new char[512]; strcpy(str, "");
        #ifndef HAVE_NCURSES
            strcat(str, "\t");
        #endif
        unsigned int relTime;
        relTime = 100*c->predictor->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->predictor->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->grid->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->grid->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->link_list->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->link_list->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->rates->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u (%g) ", c->rates->name(), relTime, c->rates->profileTime());
            strcat(str, aux);
        }
        relTime = 100*c->elastic_bounce->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->elastic_bounce->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->de_Leffe->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->de_Leffe->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->dens_int->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->dens_int->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->corrector->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->corrector->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->domain->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->domain->name(), relTime);
            strcat(str, aux);
        }
        relTime = 100*c->time_step->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->time_step->name(), relTime);
            strcat(str, aux);
        }
        relTime = 0;
        for(i=0;i<c->motions.size();i++){
            relTime += c->motions.at(i)->profileTime();
        }
        relTime = 100*relTime/Sum;
        if(relTime){
            sprintf(aux, "Motions=%u ", relTime);
            strcat(str, aux);
        }
        relTime = 100.f*c->sensors->profileTime()/Sum;
        if(relTime){
            sprintf(aux, "%s=%u ", c->sensors->name(), relTime);
            strcat(str, aux);
        }
        #ifdef HAVE_NCURSES
            move(3, 2);
            printw(str);
        #else
            strcat(str, "\n");
            printf(str);
        #endif
        delete[] str; str=0;
        delete[] aux; aux=0;
    #endif

    #ifdef HAVE_NCURSES
        attron(A_NORMAL | COLOR_PAIR(5));
        move(5, 2);
        printw("Number of particles: %u", c->n);
        move(6, 2);
        printw("Number of cells (allocated in memory): %u", c->num_cells_allocated);
        move(7, 2);
        printw("Allocated memory: %lu bytes", (int)c->allocated_mem);
    #else
        printf("\tNumber of particles: %u,", c->n);
        printf("\tAllocated memory: %lu bytes\n", (int)c->allocated_mem);
    #endif

    #ifdef HAVE_NCURSES
        attron(A_NORMAL | COLOR_PAIR(3));
    #endif

    c->energy();
    c->bounds();
    #ifdef HAVE_NCURSES
        move(10, 2);
        printw("U = %g (J)",c->eint);
        move(10, 32);
        printw("Ekin = %g (J)",c->ekin);
        move(10, 62);
        printw("E = %g (J)",c->etot);
        move(11, 2);
        printw("xmin = %g (m)",c->min_fluid_bound.x);
        move(11, 32);
        printw("ymin = %g (m)",c->min_fluid_bound.y);
        #ifdef HAVE_3D
            move(11, 62);
            printw("zmin = %g (m)",c->min_fluid_bound.z);
        #endif
        move(12, 2);
        printw("xmax = %g (m)",c->max_fluid_bound.x);
        move(12, 32);
        printw("ymax = %g (m)",c->max_fluid_bound.y);
        #ifdef HAVE_3D
            move(12, 62);
            printw("zmax = %g (m)",c->max_fluid_bound.z);
        #endif
        move(13, 2);
        printw("Vmin = %g (m/s)",c->min_v);
        move(13, 32);
        printw("Vmax = %g (m/s)",c->max_v);
    #else
        printf("\tU = %g (J), ", c->eint);
        printf("Ekin = %g (J), ", c->ekin);
        printf("E = %g (J)\n", c->etot);
        printf("\txmin = %g (m), ",c->min_fluid_bound.x);
        printf("ymin = %g (m)",c->min_fluid_bound.y);
        #ifdef HAVE_3D
            printf(", zmin = %g (m)\n",c->min_fluid_bound.z);
        #else
            printf("\n");
        #endif
        printf("\txmax = %g (m), ",c->max_fluid_bound.x);
        printf("ymax = %g (m)",c->max_fluid_bound.y);
        #ifdef HAVE_3D
            printf(", zmax = %g (m)\n",c->max_fluid_bound.z);
        #else
            printf("\n");
        #endif
        printf("\tvmin = %g (m/s), ",c->min_v);
        printf("vmax = %g (m/s)\n",c->max_v);
    #endif

    line=16;
    i=0;
    #ifdef HAVE_NCURSES
        move(line, 2);
        attron(A_NORMAL | COLOR_PAIR(1));
        printw("--- Log registry ------------------------------");
        line++;
        while((i<_n_log) && (_c_log[i]>0))
        {
            move(line, 2);
            if(_c_log[i] == 1)        // Info
                attron(A_BOLD   | COLOR_PAIR(1));
            else if(_c_log[i] == 2)   // Warning
                attron(A_BOLD   | COLOR_PAIR(4));
            else if(_c_log[i] == 3)   // Error
                attron(A_BOLD   | COLOR_PAIR(5));
            else                    // Auxiliar messages
                attron(A_NORMAL | COLOR_PAIR(1));
            printw(_m_log[i]);
            line++;
            i++;
        }
    #endif

    #ifdef HAVE_NCURSES
        refresh();
    #else
        fflush(stdout);
    #endif
}

void ScreenManager::addMessage(int Level, const char *log, const char *func)
{
    char fname[256]; strcpy(fname, "");
    if(func){
        sprintf(fname, "(%s): ", func);
    }
    // Compatibility mode when screen manager have been closed
    if(!ScreenManager::singleton()){
        if(Level == 1)
            printf("INFO ");
        else if(Level == 2)
            printf("WARNING ");
        else if(Level == 3)
            printf("ERROR ");
        printf("%s%s", fname, log);
        fflush(stdout);
        if(FileManager::singleton()){
            FILE *LogFileID = FileManager::singleton()->logFile();
            if(LogFileID){
                if(Level == 1)        // Log
                    fprintf(LogFileID, "<b><font color=\"#000000\">[INFO] %s%s</font></b><br>", fname, log);
                else if(Level == 2)        // Warning
                    fprintf(LogFileID, "<b><font color=\"#ff9900\">[WARNING] %s%s</font></b><br>", fname, log);
                else if(Level == 3)        // Error
                    fprintf(LogFileID, "<b><font color=\"#dd0000\">[ERROR] %s%s</font></b><br>", fname, log);
                else{
                    fprintf(LogFileID, "<font color=\"#000000\">%s%s</font>", fname, log);
                    if(log[strlen(log)-1] == '\n')
                    fprintf(LogFileID, "<br>");
                }
                fflush(LogFileID);
            }
        }
        return;
    }

    int i=0;
    TimeManager *T = TimeManager::singleton();
    while((i<_n_log) && (_c_log[i]>0))
    {
        i++;
    }
    if(i >= _n_log)
    {
        for(i=0;i<_n_log-1;i++)
        {
            _c_log[i] = _c_log[i+1];
            strcpy(_m_log[i], _m_log[i+1]);
        }
        _c_log[_n_log-1] = Level;
        strcpy(_m_log[_n_log-1], log);
    }
    else
    {
        _c_log[i] = Level;
        strcpy(_m_log[i], log);
    }
    // Write the message into a file
    FILE *LogFileID = FileManager::singleton()->logFile();
    if(LogFileID)
    {
        if(_c_log[i] == 1)    // Log
            fprintf(LogFileID, "<b><font color=\"#000000\">[INFO (T=%gs, Step=%d)] %s%s</font></b><br>",
                    T->time() - T->dt(), T->step(), fname, log);
        else if(_c_log[i] == 2)    // Warning
            fprintf(LogFileID, "<b><font color=\"#ff9900\">[WARNING (T=%gs, Step=%d)] %s%s</font></b><br>",
                    T->time() - T->dt(), T->step(), fname, log);
        else if(_c_log[i] == 3)    // Error
            fprintf(LogFileID, "<b><font color=\"#dd0000\">[ERROR (T=%gs, Step=%d)] %s%s</font></b><br>",
                    T->time() - T->dt(), T->step(), fname, log);
        else{
            fprintf(LogFileID, "<font color=\"#000000\">%s</font>", fname, log);
            if(log[strlen(log)-1] == '\n')
            fprintf(LogFileID, "<br>");
        }
        fflush(LogFileID);
    }
    #ifndef HAVE_NCURSES
        if(Level == 1)
            printf("\t\tINFO %s%s", fname, log);
        else if(Level == 2)
            printf("\t\tWARNING %s%s", fname, log);
        else if(Level == 3)
            printf("\t\tERROR %s%s", fname, log);
        else
            printf("\t\t%s%s", fname, log);
        fflush(stdout);
    #endif
}

void ScreenManager::printDate()
{
    char msg[512];
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    sprintf(msg, "%s\n", ctime(&seconds));
    addMessage(0, msg);
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

}}  // namespace
