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

#ifndef TIMEMANAGER_H_INCLUDED
#define TIMEMANAGER_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------------
// Include Singleton abstract class
// ----------------------------------------------------------------------------
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class TimeManager TimeManager.h TimeManager.h
 * @brief Simulation timing manager.
 */
struct TimeManager : public Aqua::Singleton<Aqua::InputOutput::TimeManager>
{
public:
	/** Constructor
	 */
	TimeManager();

	/** Destructor
	 */
	~TimeManager();

	/** Update method.
	 * @param DT Simulation time elapsed in this frame.
	 * @remarks This method must be called every step.
	 */
	void update(float DT);

	/** Must stop the simulation.
	 * @return true if simulation end has been reached.
	 */
	bool mustStop();

	/** Must print a log file.
	 * @return true if a log file must be printed.
	 * @remarks This method asumes that a file will be printed if returns 1, so
	 * control flag will be reset, aiming to the next time that need print a file.
	 */
	bool mustPrintLog();

	/** Must print an energy report.
	 * @return true if an energy file must be printed.
	 * @remarks This method asumes that a file will be printed if returns 1, so
	 * control flag will be reset, aiming to the next time that need print a file.
	 */
	bool mustPrintEnergy();

	/** Must print a Bounds file.
	 * @return true if a bounds file must be printed.
	 * @remarks This method asumes that a file will be printed if returns 1, so
	 * control flag will be reset, aiming to the next time that need print a file.
	 */
	bool mustPrintBounds();

	/** Must print a output file.
	 * @return true if a output file must be printed.
	 * @remarks This method asumes that a file will be printed if returns 1, so
	 * control flag will be reset, aiming to the next time that need print a file.
	 */
	bool mustPrintOutput();

	/** Set the simulation step.
	 * @param s Simulation step.
	 */
	void step(int s){iStep = s;}
	/** Get the simulation step.
	 * @return Simulation step.
	 */
	unsigned int step(){return iStep;}
	/** Set the simulation time.
	 * @param t Simulation time.
	 */
	void time(float t){mTime = t;}
	/** Get the simulation time.
	 * @return Simulation time.
	 */
	float time(){return mTime;}
	/** Set the simulation frame.
	 * @param frame Simulation frame.
	 */
	void frame(int frame){iFrame = frame;}
	/** Get the simulation frame.
	 * @return Simulation frame.
	 */
	unsigned int frame(){return iFrame;}
	/** Set the simulation time step.
	 * @param t Simulation time step.
	 */
	void dt(float t){mDt = t;}
	/** Get the simulation time step.
	 * @return Simulation time step.
	 */
	float dt(){return mDt;}
	/** Set the simulation inital frame.
	 * @param start intial frame.
	 */
	void startFrame(int start){zeroFrame=start;}
	/** Get the simulation inital frame.
	 * @return Simulation intial frame.
	 */
	int startFrame(){return zeroFrame;}
	/** Set the simulation inital time.
	 * @param start Simulation intial time.
	 */
	void startTime(float start){zeroTime=start;}
	/** Get the simulation inital time.
	 * @return Simulation intial time.
	 */
	float startTime(){return zeroTime;}

	/** Set the output step.
	 * @param s last output step.
	 */
	void outputStep(int s){mOutputStep = s;}
	/** Get the output step.
	 * @return Last output step.
	 */
	int outputStep(){return mOutputStep;}
	/** Get the output iteration per frame.
	 * @return Iterations per frame.
	 */
	int outputIPF(){return mOutputIPF;}
	/** Set the output time.
	 * @param t last output time.
	 */
	void outputTime(float t){mOutputTime = t;}
	/** Get the output time.
	 * @return Last output time.
	 */
	float outputTime(){return mOutputTime;}
	/** Get the output frames per second.
	 * @return Frames per second.
	 */
	float outputFPS(){return mOutputFPS;}

	/** Get maximum simulation time.
	 * @return maximum simulation time.
	 */
	float maxTime(){return mSimMaxTime;}
	/** Get maximum simulation frame.
	 * @return maximum simulation frame.
	 */
	int maxFrame(){return mSimMaxFrames;}

private:
	/// Actual step
	unsigned int iStep;
	/// Actual time
	float mTime;
	/// Actual frame
	unsigned int iFrame;
	/// Start frame
	float zeroTime;
	/// Start frame
	int zeroFrame;
	/// Maximum time into simulation (-1 if simulation don't stop by time criteria)
	float mSimMaxTime;
	/// Maximum number of steps into simulation (-1 if simulation don't stop by steps criteria)
	int mSimMaxSteps;
	/// Maximum number of frames into simulation (-1 if simulation don't stop by frames criteria)
	int mSimMaxFrames;
	/// Time step
	float mDt;
	/// Convection time step viscosity
	float mSigma;

	/// Time when last log file printed
	float mLogTime;
	/// FPS for Log files (-1 if Log file must not be printed)
	float mLogFPS;
	/// Step when last log file printed
	int mLogStep;
	/// IPF for Log files (-1 if Log file must not be printed)
	int mLogIPF;

	/// Time when last Energy file printed
	float mReTime;
	/// FPS for Energy files (-1 if Energy file must not be printed)
	float mReFPS;
	/// Step when last Energy file printed
	int mReStep;
	/// IPF for Energy files (-1 if Energy file must not be printed)
	int mReIPF;

	/// Time when last Bounds file printed
	float mBoundsTime;
	/// FPS for Bounds files (-1 if Bounds file must not be printed)
	float mBoundsFPS;
	/// Step when last Bounds file printed
	int mBoundsStep;
	/// IPF for Bounds files (-1 if Bounds file must not be printed)
	int mBoundsIPF;

	/// Time when last Output file printed
	float mOutputTime;
	/// FPS for Output files (-1 if Output file must not be printed)
	float mOutputFPS;
	/// Step when last Output file printed
	int mOutputStep;
	/// IPF for Output files (-1 if Output file must not be printed)
	int mOutputIPF;
};

}}  // namespace

#endif // TIMEMANAGER_H_INCLUDED
