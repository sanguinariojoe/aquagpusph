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

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the problem setup header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

TimeManager::TimeManager()
	: iStep(0)
	, mTime(0.f)
	, iFrame(0)
	, zeroTime(0.f)
	, zeroFrame(0)
	, mSimMaxTime(-1.f)
	, mSimMaxSteps(-1)
	, mSimMaxFrames(-1)
	, mDt(0.f)
	, mSigma(0.f)
	, mLogTime(0.f)
	, mLogFPS(-1.f)
	, mLogStep(0)
	, mLogIPF(-1)
	, mReTime(0.f)
	, mReFPS(-1.f)
	, mReStep(0)
	, mReIPF(-1)
	, mBoundsTime(0.f)
	, mBoundsFPS(-1.f)
	, mBoundsStep(0)
	, mBoundsIPF(-1)
	, mOutputTime(0.f)
	, mOutputFPS(-1.f)
	, mOutputStep(0)
	, mOutputIPF(-1)
{
	ProblemSetup *P = ProblemSetup::singleton();
	ScreenManager *S = ScreenManager::singleton();
	//! Take simulation end criteria
	unsigned int Mode = P->time_opts.sim_end_mode;
	if(Mode & __FRAME_MODE__) {
		mSimMaxFrames = P->time_opts.sim_end_frame;
	}
	if(Mode & __ITER_MODE__) {
		mSimMaxSteps = P->time_opts.sim_end_step;
	}
	if(Mode & __TIME_MODE__) {
		mSimMaxTime = P->time_opts.sim_end_time;
	}
	//! Take log file print criteria
	Mode = P->time_opts.log_mode;
	if(Mode >= __IPF_MODE__)
	{
		Mode -= __IPF_MODE__;
		mLogIPF = P->time_opts.log_ipf;
	}
	if(Mode >= __FPS_MODE__)
	{
		Mode -= __FPS_MODE__;
		mLogFPS = P->time_opts.log_fps;
	}
	//! Take Energy file print criteria
	Mode = P->time_opts.energy_mode;
	if(Mode >= __IPF_MODE__)
	{
		Mode -= __IPF_MODE__;
		mReIPF = P->time_opts.energy_ipf;
	}
	if(Mode >= __FPS_MODE__)
	{
		Mode -= __FPS_MODE__;
		mReFPS = P->time_opts.energy_fps;
	}
	//! Take Bounds file print criteria
	Mode = P->time_opts.bounds_mode;
	if(Mode >= __IPF_MODE__)
	{
		Mode -= __IPF_MODE__;
		mBoundsIPF = P->time_opts.bounds_ipf;
	}
	if(Mode >= __FPS_MODE__)
	{
		Mode -= __FPS_MODE__;
		mBoundsFPS = P->time_opts.bounds_fps;
	}
	//! Take Output file print criteria
	Mode = P->time_opts.output_mode;
	if(Mode >= __IPF_MODE__)
	{
		Mode -= __IPF_MODE__;
		mOutputIPF = P->time_opts.output_ipf;
	}
	if(Mode >= __FPS_MODE__)
	{
		Mode -= __FPS_MODE__;
		mOutputFPS = P->time_opts.output_fps;
	}
	//! Set stabilization time
	mTime -= P->time_opts.stabilization_time;
	zeroTime -= P->time_opts.stabilization_time;
	S->addMessage(1, "(TimeManager::TimeManager): Time manager built OK.\n");
}

TimeManager::~TimeManager()
{
}

void TimeManager::update(float DT)
{
	mDt = DT;
	iStep++;
	mTime += mDt;
}

bool TimeManager::mustStop()
{
	if( (mSimMaxTime >= 0.f) && (mTime >= mSimMaxTime) )
		return true;
	if( (mSimMaxSteps >= 0) && (iStep >= mSimMaxSteps) )
		return true;
	if( (mSimMaxFrames >= 0) && (iFrame >= mSimMaxFrames) )
		return true;
	return false;
}

bool TimeManager::mustPrintLog()
{
	if( ( (mLogFPS >= 0.f) || (mLogIPF >= 0.f) ) && (iFrame==1) && (iStep==1) )
	{
		mLogTime = mTime;
		mLogStep = iStep;
		return true;
	}
	if( (mLogFPS >= 0.f) && (mTime - mLogTime >= 1.f/mLogFPS) )
	{
		mLogTime += 1.f/mLogFPS;
		mLogStep = iStep;
		return true;
	}
	if( (mLogIPF > 0) && (iStep - mLogStep >= mLogIPF) )
	{
		mLogTime = mTime;
		mLogStep = iStep;
		return true;
	}
	return false;
}

bool TimeManager::mustPrintEnergy()
{
	if( ( (mReFPS >= 0.f) || (mReIPF >= 0.f) ) && (iFrame==1) && (iStep==1) )
	{
		mReTime = mTime;
		mReStep = iStep;
		return true;
	}
	if( (mReFPS >= 0.f) && (mTime - mReTime >= 1.f/mReFPS) )
	{
		mReTime += 1.f/mReFPS;
		mReStep = iStep;
		return true;
	}
	if( (mReIPF > 0) && (iStep - mReStep >= mReIPF) )
	{
		mReTime = mTime;
		mReStep = iStep;
		return true;
	}
	return false;
}

bool TimeManager::mustPrintBounds()
{
	if( ( (mBoundsFPS >= 0.f) || (mBoundsIPF >= 0.f) ) && (iFrame==1) && (iStep==1) )
	{
		mBoundsTime = mTime;
		mBoundsStep = iStep;
		return true;
	}
	if( (mBoundsFPS >= 0.f) && (mTime - mBoundsTime >= 1.f/mBoundsFPS) )
	{
		mBoundsTime += 1.f/mBoundsFPS;
		mBoundsStep = iStep;
		return true;
	}
	if( (mBoundsIPF > 0) && (iStep - mBoundsStep >= mBoundsIPF) )
	{
		mBoundsTime = mTime;
		mBoundsStep = iStep;
		return true;
	}
	return false;
}

bool TimeManager::mustPrintOutput()
{
	if(mTime < 0.f){
	    iStep = 0;
	    return false;
	}
	if( ( (mOutputFPS >= 0.f) || (mOutputIPF >= 0.f) ) && (iFrame==0) && (iStep==1) )
	{
		mOutputTime = mTime;
		mOutputStep = iStep;
		iFrame++;
		return true;
	}
	if( (mOutputFPS > 0.f) && (mTime - mOutputTime >= 1.f/mOutputFPS) )
	{
		mOutputTime += 1.f/mOutputFPS;
		mOutputStep = iStep;
		iFrame++;
		return true;
	}
	if( (mOutputIPF > 0) && (iStep - mOutputStep >= mOutputIPF) )
	{
		mOutputTime = mTime;
		mOutputStep = iStep;
		iFrame++;
		return true;
	}
	// Special case about the simulation was end
	if(mustStop()){
        return true;
	}

	return false;
}

}}  // namespace
