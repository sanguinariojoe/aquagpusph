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
 * @brief Performance report.
 * (See Aqua::CalcServer::Reports::Performance for details)
 */

#include <algorithm>
#include <iomanip>
#include <CalcServer/Reports/Performance.h>
#include <CalcServer.h>
#include <InputOutput/Logger.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Performance::Performance(const std::string tool_name,
                         const std::string color,
                         bool bold,
                         const std::string output_file)
    : Report(tool_name, "dummy_fields_string")
    , _color(color)
    , _bold(bold)
    , _output_file(output_file)
    , _first_execution(true)
{
    gettimeofday(&_tic, NULL);
}

Performance::~Performance()
{
    if(_f.is_open()) _f.close();
}

void Performance::setup()
{
    unsigned int i;

    std::ostringstream msg;
    msg << "Loading the report \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    // Set the color in lowercase
    std::transform(_color.begin(), _color.end(), _color.begin(), ::tolower);

    // Open the output file
    if(_output_file.compare("")) {
        _f.open(_output_file.c_str(), std::ios::out);
        // Write the header
        _f << "# t elapsed average(elapsed) variance(elapsed) "
           << "overhead average(overhead) variance(overhead) progress ETA"
           << std::endl;
    }
}

size_t Performance::computeAllocatedMemory(){
    unsigned int i;
    size_t allocated_mem = 0;
    CalcServer *C = CalcServer::singleton();

    // Get the allocated memory in the variables
    InputOutput::Variables *vars = C->variables();
    allocated_mem += vars->allocatedMemory();

    // Gwet the additionally allocated memory in the tools
    std::vector<Tool*> tools = C->tools();
    for(i = 0; i < tools.size(); i++){
        allocated_mem += tools.at(i)->allocatedMemory();
    }

    return allocated_mem;
}

void Performance::_execute()
{
    unsigned int i;
    CalcServer *C = CalcServer::singleton();
    std::stringstream data;

    size_t allocated_MB = computeAllocatedMemory() / (1024 * 1024);
    data << "Performance:" << std::endl
         << "Memory=" << std::setw(2) << allocated_MB << "MB" << std::endl;

    // Add the tools time elapsed
    std::vector<Tool*> tools = C->tools();
    float elapsed = 0.f;
    float elapsed_ave = 0.f;
    for(i = 0; i < tools.size(); i++){
        // Exclude the tool itself
        if(this == tools.at(i)){
            continue;
        }
        elapsed += tools.at(i)->elapsedTime(false);
        elapsed_ave += tools.at(i)->elapsedTime();
    }

    timeval tac;
    gettimeofday(&tac, NULL);
    float elapsed_seconds;
    elapsed_seconds = (float)(tac.tv_sec - _tic.tv_sec);
    elapsed_seconds += (float)(tac.tv_usec - _tic.tv_usec) * 1E-6f;
    gettimeofday(&_tic, NULL);
    if(_first_execution){
        _first_execution = false;
        // Purge out all the tool building required time
        elapsed_seconds = elapsed;
    }
    addElapsedTime(elapsed_seconds);

    data << "Elapsed=" << std::setw(16) << elapsedTime()
         << "s (+-" << std::setw(16) << elapsedTimeVariance()
         << "s)" << std::endl;
    data << "Overhead=" << std::setw(16) << elapsedTime() - elapsed_ave
         << "s" << std::endl;

    // Compute the progress
    InputOutput::Variables *vars = C->variables();
    float progress = 0.f;
    float t = *(float *)vars->get("t")->get();
    float end_t = *(float *)vars->get("end_t")->get();
    progress = max(progress, t / end_t);
    unsigned int iter = *(unsigned int *)vars->get("iter")->get();
    unsigned int end_iter = *(unsigned int *)vars->get("end_iter")->get();
    progress = max(progress, (float)iter / end_iter);
    unsigned int frame = *(unsigned int *)vars->get("frame")->get();
    unsigned int end_frame = *(unsigned int *)vars->get("end_frame")->get();
    progress = max(progress, (float)frame / end_frame);

    // And the estimated time to arrive
    float total_elapsed = elapsedTime() * used_times();
    float ETA = total_elapsed * (1.f / progress - 1.f);

    data << "Percentage=" << std::setw(16) << std::setprecision(2)
         << progress * 100.f << "\tETA=" << std::setw(16) << ETA << std::endl;

    // Replace the trailing space by a line break
    if(data.str().back() == ' ') {
        data.seekp(-1, data.cur);
        data << std::endl;
    }

    InputOutput::Logger::singleton()->writeReport(data.str(),
                                                         _color,
                                                         _bold);

    // Write the output file
    if(_f.is_open()){
        _f << t << " "
           << elapsedTime(false) << " " << elapsedTime() << " "
           << elapsedTimeVariance() << " "
           << elapsedTime(false) - elapsed << " "
           << elapsedTime() - elapsed_ave << " "
           << progress * 100.f << " " << ETA << std::endl;
    }
}

}}} // namespace
