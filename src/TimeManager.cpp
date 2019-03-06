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

#include <string>
#include <TimeManager.h>
#include <CalcServer.h>
#include <InputOutput/Logger.h>

namespace Aqua{ namespace InputOutput{

TimeManager::TimeManager(const ProblemSetup& sim_data)
    : _step(NULL)
    , _time(NULL)
    , _dt(NULL)
    , _frame_var(NULL)
    , _frame(NULL)
    , _time_max(NULL)
    , _steps_max(NULL)
    , _frames_max(NULL)
    , _output_time(0.f)
    , _output_fps(-1.f)
    , _output_step(0)
    , _output_ipf(-1)
{
    unsigned int i;
    Variables *vars = CalcServer::CalcServer::singleton()->variables();
    // Check the variables validity
    std::map<std::string, std::string> var_types {
        {"t", "float"},
        {"dt", "float"},
        {"iter", "unsigned int"},
        {"frame", "unsigned int"},
        {"end_t", "float"},
        {"end_iter", "unsigned int"},
        {"end_frame", "unsigned int"}
    };
    for (auto& type : var_types) {
        if(vars->get(type.first)->type().compare(type.second)){
            std::ostringstream msg;
            msg << "Variable \"" << type.first
                << "\" has been redeclared as \"" << type.second
                << "\" instead of the expected type, \""
                << vars->get(type.first)->type() << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
    }

    _time = (float *)vars->get("t")->get();
    _dt = (float *)vars->get("dt")->get();
    _step = (unsigned int *)vars->get("iter")->get();
    _frame_var = (UIntVariable*)vars->get("frame");
    _frame = (unsigned int *)(_frame_var->get());
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

    _output_time = *_time;
    _output_step = *_step;
}

TimeManager::~TimeManager()
{
}

const bool TimeManager::mustStop() const
{
    if((time() >= maxTime()) ||
       (step() >= maxStep()) ||
       (frame() >= maxFrame())) {
        return true;
    }

    return false;
}

const bool TimeManager::mustPrintOutput()
{
    if(time() < 0.f){
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
