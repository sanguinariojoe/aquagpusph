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
 * @brief Log file manager.
 * (See Aqua::InputOutput::Log for details)
 */

#include <stdlib.h>
#include <string.h>

#include <InputOutput/Log.h>
#include <ScreenManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Log::Log()
    : _file(NULL)
{
    if(create())
        exit(EXIT_FAILURE);
}

Log::~Log()
{
    close();
}

bool Log::save()
{
    if(!_file){
        return true;
    }

    return false;
}

bool Log::create()
{
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();

    if(file("log.%d.html", 0)){
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessageF(0, "\tHow do you received this message?.\n");
        return true;
    }
    _file = fopen(file(), "w");
    if(!_file){
        sprintf(msg,
                "Failure creating the log file \"%s\"\n",
                file());
        return true;
    }

    fprintf(_file, "<html>\n");
    fprintf(_file, "<head><title>AQUAgpusph log file.</title></head>\n");
    fprintf(_file, "<body bgcolor=\"#f0ffff\">\n");
    fprintf(_file, "<h1 align=\"center\">AQUAgpusph log file.</h1>\n");
    // Starting data
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    fprintf(_file, "<p align=\"left\">%s</p>\n", ctime(&seconds));
    fprintf(_file, "<hr><br>\n");
    fflush(_file);

    return false;
}

bool Log::close()
{
    if(!_file)
        return true;

    unsigned int i;
    char msg[512];
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    float fluid_mass = 0.f;
    int err_code = CL_SUCCESS;
    strcpy(msg, "");

    fprintf(_file, "<br><hr>\n");
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    fprintf(_file,
            "<b><font color=\"#000000\">End of simulation</font></b><br>\n");
    fprintf(_file, "<p align=\"left\">%s</p>\n", ctime(&seconds));
    fprintf(_file, "</body>\n");
    fprintf(_file, "</html>\n");

    fflush(_file);
    fclose(_file);
    _file = NULL;

    return false;
}

}}  // namespace
