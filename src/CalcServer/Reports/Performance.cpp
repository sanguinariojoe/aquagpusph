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

#include <CalcServer/Reports/Performance.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Performance::Performance(const char* tool_name,
                         const char* color,
                         bool bold,
                         const char* output_file)
    : Report(tool_name, "dummy_fields_string")
    , _color(NULL)
    , _bold(bold)
    , _output_file(NULL)
    , _f(NULL)
    , _first_execution(true)
{
    _color = new char[strlen(color) + 1];
    strcpy(_color, color);
    _output_file = new char[strlen(output_file) + 1];
    strcpy(_output_file, output_file);

    gettimeofday(&_tic, NULL);
}

Performance::~Performance()
{
    if(_color) delete[] _color; _color=NULL;
    if(_output_file) delete[] _output_file; _output_file=NULL;
    if(_f) fclose(_f); _f = NULL;
}

bool Performance::setup()
{
    unsigned int i;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the report \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    // Set the color in lowercase
    for(i = 0; i < strlen(_color); i++){
        _color[i] = tolower(_color[i]);
    }

    // Open the output file
    if(strcmp(_output_file, "")){
        _f = fopen(_output_file, "w");
        if(!_f){
            sprintf(msg,
                    "The file \"%s\" cannot be written\n",
                    _output_file);
            S->addMessageF(L_ERROR, msg);
            return true;
        }
        // Write the header
        fprintf(_f,
                "# t elapsed average(elapsed) variance(elapsed) overhead average(overhead) variance(overhead) progress ETA\n");
    }

    return false;
}

size_t Performance::computeAllocatedMemory(){
    unsigned int i;
    size_t allocated_mem = 0;
    CalcServer *C = CalcServer::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    // Get the allocated memory in the variables
    InputOutput::Variables vars = C->variables();
    allocated_mem += vars.allocatedMemory();

    // Gwet the additionally allocated memory in the tools
    std::vector<Tool*> tools = C->tools();
    for(i = 0; i < tools.size(); i++){
        allocated_mem += tools.at(i)->allocatedMemory();
    }

    return allocated_mem;
}

bool Performance::_execute()
{
    unsigned int i;
    CalcServer *C = CalcServer::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char data[4096];

    size_t allocated_MB = computeAllocatedMemory() / (1024 * 1024);
    sprintf(data, "Performance:\nMemory=%16luMB\n", allocated_MB);

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

    sprintf(data + strlen(data), "Elapsed=%16gs (+-%16gs)\n",
            elapsedTime(),
            elapsedTimeVariance());

    sprintf(data + strlen(data), "Overhead=%16gs\n",
            elapsedTime() - elapsed_ave);

    // Compute the progress
    InputOutput::Variables vars = C->variables();
    float progress = 0.f;
    float t = *(float *)vars.get("t")->get();
    float end_t = *(float *)vars.get("end_t")->get();
    progress = max(progress, t / end_t);
    unsigned int iter = *(unsigned int *)vars.get("iter")->get();
    unsigned int end_iter = *(unsigned int *)vars.get("end_iter")->get();
    progress = max(progress, (float)iter / end_iter);
    unsigned int frame = *(unsigned int *)vars.get("frame")->get();
    unsigned int end_frame = *(unsigned int *)vars.get("end_frame")->get();
    progress = max(progress, (float)frame / end_frame);

    // And the estimated time to arrive
    float total_elapsed = elapsedTime() * used_times();
    float ETA = total_elapsed * (1.f / progress - 1.f);

    sprintf(data + strlen(data), "Percentage=%16.2f\tETA=%16gs\n",
            progress * 100.f,
            ETA);

    // Replace the trailing space by a line break
    if(data[strlen(data) - 1] == ' ')
        data[strlen(data) - 1] = '\n';

    S->writeReport(data, _color, _bold);

    // Write the output file
    if(_f){
        fprintf(_f,
                "%16g %16g %16g %16g %16g %16g %16.2f %16g\n",
                t,
                elapsedTime(false),
                elapsedTime(),
                elapsedTimeVariance(),
                elapsedTime(false) - elapsed,
                elapsedTime() - elapsed_ave,
                progress * 100.f,
                ETA);
    }

    return false;
}

}}} // namespace
