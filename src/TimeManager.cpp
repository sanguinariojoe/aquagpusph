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
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

TimeManager::TimeManager()
    : _step(NULL)
    , _time(NULL)
    , _dt(NULL)
    , _frame(0)
    , _start_time(0.f)
    , _start_frame(0)
    , _time_max(-1.f)
    , _steps_max(-1)
    , _frames_max(-1)
    , _log_time(0.f)
    , _log_fps(-1.f)
    , _log_step(0)
    , _log_ipf(-1)
    , _output_time(0.f)
    , _output_fps(-1.f)
    , _output_step(0)
    , _output_ipf(-1)
{
    char msg[1024];
    ProblemSetup *P = ProblemSetup::singleton();
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    Variables* vars = C->variables();
    if(strcmp(vars->get("iter")->type(), "unsigned int")){
        sprintf(msg,
                "Expected a variable \"%s\" of type \"%s\", but \"%s\" one was found\n",
                "iter",
                "unsigned int",
                vars->get("iter")->type());
        S->addMessageF(3, msg);
        exit(EXIT_FAILURE);
    }
    _step = (unsigned int *)vars->get("iter")->get();

    if(strcmp(vars->get("t")->type(), "float")){
        sprintf(msg,
                "Expected a variable \"%s\" of type \"%s\", but \"%s\" one was found\n",
                "t",
                "float",
                vars->get("t")->type());
        S->addMessageF(3, msg);
        exit(EXIT_FAILURE);
    }
    _time = (float *)vars->get("t")->get();

    if(strcmp(vars->get("dt")->type(), "float")){
        sprintf(msg,
                "Expected a variable \"%s\" of type \"%s\", but \"%s\" one was found\n",
                "dt",
                "float",
                vars->get("dt")->type());
        S->addMessageF(3, msg);
        exit(EXIT_FAILURE);
    }
    _dt = (float *)vars->get("dt")->get();

    unsigned int mode = P->time_opts.sim_end_mode;
    if(mode & __FRAME_MODE__) {
        _frames_max = P->time_opts.sim_end_frame;
    }
    if(mode & __ITER_MODE__) {
        _steps_max = P->time_opts.sim_end_step;
    }
    if(mode & __TIME_MODE__) {
        _time_max = P->time_opts.sim_end_time;
    }

    mode = P->time_opts.log_mode;
    if(mode >= __IPF_MODE__)
    {
        mode -= __IPF_MODE__;
        _log_ipf = P->time_opts.log_ipf;
    }
    if(mode >= __FPS_MODE__)
    {
        mode -= __FPS_MODE__;
        _log_fps = P->time_opts.log_fps;
    }

    mode = P->time_opts.output_mode;
    if(mode >= __IPF_MODE__)
    {
        mode -= __IPF_MODE__;
        _output_ipf = P->time_opts.output_ipf;
    }
    if(mode >= __FPS_MODE__)
    {
        mode -= __FPS_MODE__;
        _output_fps = P->time_opts.output_fps;
    }

    *_step = P->time_opts.step0;
    *_dt = P->time_opts.dt0;
    *_time = P->time_opts.t0;
    _start_time = P->time_opts.t0;
    _frame = P->time_opts.frame0;

    if(*_time > 0.f){
        _log_time = *_time;
        _log_step = *_step;
        _output_time = *_time;
        _output_step = *_step;
    }

    S->addMessageF(1, "Time manager built OK.\n");
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
    if( (_time_max >= 0.f) && (time() >= _time_max) )
        return true;
    if( (_steps_max >= 0) && (step() >= _steps_max) )
        return true;
    if( (_frames_max >= 0) && (_frame >= _frames_max) )
        return true;
    return false;
}

bool TimeManager::mustPrintLog()
{
    if( ( (_log_fps >= 0.f) || (_log_ipf >= 0.f) ) && (_frame==1) && (step()==1) ) {
        _log_time = time();
        _log_step = step();
        return true;
    }
    if( (_log_fps >= 0.f) && (time() - _log_time >= 1.f/_log_fps) ) {
        _log_time += 1.f/_log_fps;
        _log_step = step();
        return true;
    }
    if( (_log_ipf > 0) && (step() - _log_step >= _log_ipf) ) {
        _log_time = time();
        _log_step = step();
         return true;
    }
    return false;
}

bool TimeManager::mustPrintOutput()
{
    if(time() < 0.f){
        step(0);
        return false;
    }
    if( ( (_output_fps >= 0.f) || (_output_ipf >= 0.f) ) && (_frame==0) && (step()==1) ) {
        _output_time = time();
        _output_step = step();
        _frame++;
        return true;
    }
    if( (_output_fps > 0.f) && (time() - _output_time >= 1.f/_output_fps) ) {
        _output_time += 1.f/_output_fps;
        _output_step = step();
        _frame++;
        return true;
    }
    if( (_output_ipf > 0) && (step() - _output_step >= _output_ipf) ) {
        _output_time = time();
        _output_step = step();
        _frame++;
        return true;
    }
    // We are interested into print an output in the simulation end state
    if(mustStop()){
        return true;
    }

    return false;
}

}}  // namespace
