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
 * @brief Simulation time flow events manager.
 * (See Aqua::InputOutput::TimeManager for details)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <TimeManager.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

TimeManager::TimeManager(ProblemSetup sim_data)
    : _step(NULL)
    , _time(NULL)
    , _dt(NULL)
    , _frame(NULL)
    , _start_time(0.f)
    , _start_frame(0)
    , _time_max(NULL)
    , _steps_max(NULL)
    , _frames_max(NULL)
    , _log_time(0.f)
    , _log_fps(-1.f)
    , _log_step(0)
    , _log_ipf(-1)
    , _output_time(0.f)
    , _output_fps(-1.f)
    , _output_step(0)
    , _output_ipf(-1)
{
    unsigned int i;
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    Variables* vars = C->variables();
    // Check the variables validity
    const unsigned int var_num = 7;
    const char* var_names[var_num] = {"t",
                                      "dt",
                                      "iter",
                                      "frame",
                                      "end_t",
                                      "end_iter",
                                      "end_frame"};
    const char* var_types[var_num] = {"float",
                                      "float",
                                      "unsigned int",
                                      "unsigned int",
                                      "float",
                                      "unsigned int",
                                      "unsigned int"};
    for(i = 0; i < var_num; i++){
        if(strcmp(vars->get(var_names[i])->type(), var_types[i])){
            sprintf(msg,
                    "Expected a variable \"%s\" of type \"%s\", but \"%s\" one was found\n",
                    var_names[i],
                    var_types[i],
                    vars->get(var_names[i])->type());
            S->addMessageF(L_ERROR, msg);
            exit(EXIT_FAILURE);
        }
    }

    _time = (float *)vars->get("t")->get();
    _dt = (float *)vars->get("dt")->get();
    _step = (unsigned int *)vars->get("iter")->get();
    _frame = (unsigned int *)vars->get("frame")->get();
    _time_max = (float *)vars->get("end_t")->get();
    _steps_max = (unsigned int *)vars->get("end_iter")->get();
    _frames_max = (unsigned int *)vars->get("end_frame")->get();

    unsigned int mode = sim_data.time_opts.sim_end_mode;
    if(mode & __FRAME_MODE__) {
        *_frames_max = sim_data.time_opts.sim_end_frame;
    }
    if(mode & __ITER_MODE__) {
        *_steps_max = sim_data.time_opts.sim_end_step;
    }
    if(mode & __TIME_MODE__) {
        *_time_max = sim_data.time_opts.sim_end_time;
    }

    mode = sim_data.time_opts.output_mode;
    if(mode >= __IPF_MODE__)
    {
        mode -= __IPF_MODE__;
        _output_ipf = sim_data.time_opts.output_ipf;
    }
    if(mode >= __FPS_MODE__)
    {
        mode -= __FPS_MODE__;
        _output_fps = sim_data.time_opts.output_fps;
    }

    *_step = 0;
    *_dt = 0.f;
    *_time = 0.f;
    _start_time = 0.f;
    *_frame = 0;

    if(*_time > 0.f){
        _log_time = *_time;
        _log_step = *_step;
        _output_time = *_time;
        _output_step = *_step;
    }

    S->addMessageF(L_INFO, "Time manager built OK.\n");
}

TimeManager::~TimeManager()
{
}

void TimeManager::update(float sim_dt)
{
    dt(sim_dt);
    step(step() + 1);
    time(time() + dt());
}

bool TimeManager::mustStop()
{
    if(time() >= maxTime())
        return true;
    if(step() >= maxStep())
        return true;
    if(frame() >= maxFrame())
        return true;
    return false;
}

bool TimeManager::mustPrintOutput()
{
    if(time() < 0.f){
        step(0);
        return false;
    }
    if( ( (_output_fps >= 0.f) || (_output_ipf >= 0.f) ) && (frame()==0) && (step()==1) ) {
        _output_time = time();
        _output_step = step();
        frame(frame() + 1);
        return true;
    }
    if( (_output_fps > 0.f) && (time() - _output_time >= 1.f/_output_fps) ) {
        _output_time += 1.f/_output_fps;
        _output_step = step();
        frame(frame() + 1);
        return true;
    }
    if( (_output_ipf > 0) && (step() - _output_step >= _output_ipf) ) {
        _output_time = time();
        _output_step = step();
        frame(frame() + 1);
        return true;
    }
    // We are interested into print an output in the simulation end state
    if(mustStop()){
        return true;
    }

    return false;
}

}}  // namespace
