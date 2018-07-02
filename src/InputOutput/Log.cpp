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

#include <InputOutput/Log.h>
#include <ScreenManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Log::Log()
{
    create();
}

Log::~Log()
{
    close();
}

void Log::save()
{
    if(!_file){
        LOG(L_ERROR, "The Log file was not correctly created");
        throw std::runtime_error("Invalid Log file");
    }
}

void Log::create()
{
    file("log.%d.html", 0);
    _file.open(file().c_str());
    if(!_file.is_open()){
        std::ostringstream msg;
        msg << "Failure creating the log file \"" << file()
            << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Failure creating the Log file");
    }

    _file << "<html>" << std::endl;
    _file << "<head><title>AQUAgpusph log file.</title></head>" << std::endl;
    _file << "<body bgcolor=\"#f0ffff\">" << std::endl;
    _file << "<h1 align=\"center\">AQUAgpusph log file.</h1>" << std::endl;
    // Starting data
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    _file << "<p align=\"left\">" << ctime(&seconds) << "</p>" << std::endl;
    _file << "<hr><br>" << std::endl;
    _file.flush();
}

void Log::close()
{
    if(!_file.is_open()) {
        LOG(L_ERROR, "The Log file was not correctly created");
        throw std::runtime_error("Invalid Log file");
    }

    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    float fluid_mass = 0.f;
    int err_code = CL_SUCCESS;

    _file << "<br><hr>" << std::endl;
    struct timeval now_time;
    gettimeofday(&now_time, NULL);
    const time_t seconds = now_time.tv_sec;
    _file << "<b><font color=\"#000000\">End of simulation</font></b><br>"
          << std::endl;
    _file << "<p align=\"left\">" << ctime(&seconds) << "</p>" << std::endl;
    _file << "</body>" << std::endl;
    _file << "</html>" << std::endl;
    _file.close();
}

}}  // namespace
