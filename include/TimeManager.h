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

#include <sphPrerequisites.h>
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class TimeManager TimeManager.h TimeManager.h
 * @brief Control the simulation time, including the events triggering like the
 * ouput files updating.
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
	 * @param dt Simulation time elapsed in this frame.
	 */
	void update(float dt);

	/** Get if the simulation must be end.
	 * @return true if simulation has finished, false otherwise.
	 */
	bool mustStop();

	/** Get if the log file must be updated.
	 * @return true if the log should be updated, false otherwise.
	 * @remarks This method is assuming that when it is asked, the file will be
	 * updated, so asking it again will return false until the next file
	 * updating event.
	 */
	bool mustPrintLog();

	/** Get if the energy file must be updated.
	 * @return true if the energy should be updated, false otherwise.
	 * @remarks This method is assuming that when it is asked, the file will be
	 * updated, so asking it again will return false until the next file
	 * updating event.
	 */
	bool mustPrintEnergy();

	/** Get if the fluid bounds file must be updated.
	 * @return true if the bounds file should be updated, false otherwise.
	 * @remarks This method is assuming that when it is asked, the file will be
	 * updated, so asking it again will return false until the next file
	 * updating event.
	 */
	bool mustPrintBounds();

	/** Get if a general simulation output must be printed.
	 * @return true if an output should be printed, false otherwise.
	 * @remarks This method is assuming that when it is asked, the file will be
	 * updated, so asking it again will return false until the next file
	 * updating event.
	 */
	bool mustPrintOutput();

	/** Set the simulation step.
	 * @param s Simulation step.
	 */
	void step(unsigned int s){_step = s;}
	/** Get the simulation step.
	 * @return Simulation step.
	 */
	unsigned int step(){return _step;}
	/** Set the simulation time.
	 * @param t Simulation time.
	 */
	void time(float t){_time = t;}
	/** Get the simulation time.
	 * @return Simulation time.
	 */
	float time(){return _time;}
	/** Set the simulation frame.
	 * @param frame Simulation frame.
	 */
	void frame(unsigned int frame){_frame = frame;}
	/** Get the simulation frame.
	 * @return Simulation frame.
	 */
	unsigned int frame(){return _frame;}
	/** Set the simulation time step.
	 * @param dt Simulation time step.
	 */
	void dt(float dt){_dt = dt;}
	/** Get the simulation time step.
	 * @return Simulation time step.
	 */
	float dt(){return _dt;}
	/** Set the simulation starting frame.
	 * @param start starting frame.
	 */
	void startFrame(int start){_start_frame=start;}
	/** Get the simulation starting frame.
	 * @return Simulation starting frame.
	 */
	int startFrame(){return _start_frame;}
	/** Set the simulation starting time.
	 * @param start Simulation starting time.
	 */
	void startTime(float start){_start_time=start;}
	/** Get the simulation inital time.
	 * @return Simulation intial time.
	 */
	float startTime(){return _start_time;}

	/** Set the last output event step.
	 * @param s last output event step.
	 */
	void outputStep(int s){_output_step = s;}
	/** Get the last output event step.
	 * @return Last output event step.
	 */
	int outputStep(){return _output_step;}
	/** Get the iterations per output frame.
	 * @return Iterations per output frame.
	 */
	int outputIPF(){return _output_ipf;}
	/** Set the last output event time instant.
	 * @param t Last output event time instant.
	 */
	void outputTime(float t){_output_time = t;}
	/** Get the output time.
	 * @return Last output time.
	 */
	float outputTime(){return _output_time;}
	/** Get the output frames per second.
	 * @return Frames per second.
	 */
	float outputFPS(){return _output_fps;}

	/** Get maximum simulation time.
	 * @return maximum simulation time.
	 */
	float maxTime(){return _time_max;}
	/** Get maximum simulation frame.
	 * @return maximum simulation frame.
	 */
	int maxFrame(){return _frames_max;}

private:
	/// Actual step
	unsigned int _step;
	/// Actual time
	float _time;
	/// Actual frame
	unsigned int _frame;
	/// Start frame
	float _start_time;
	/// Start frame
	int _start_frame;
	/// Maximum time into simulation (-1 if simulation don't stop by time criteria)
	float _time_max;
	/// Maximum number of steps into simulation (-1 if simulation don't stop by steps criteria)
	int _steps_max;
	/// Maximum number of frames into simulation (-1 if simulation don't stop by frames criteria)
	int _frames_max;
	/// Time step
	float _dt;

	/// Time when last log file printed
	float _log_time;
	/// FPS for Log files (-1 if Log file must not be printed)
	float _log_fps;
	/// Step when last log file printed
	int _log_step;
	/// IPF for Log files (-1 if Log file must not be printed)
	int _log_ipf;

	/// Time when last Energy file printed
	float _en_time;
	/// FPS for Energy files (-1 if Energy file must not be printed)
	float _en_fps;
	/// Step when last Energy file printed
	int _en_step;
	/// IPF for Energy files (-1 if Energy file must not be printed)
	int _en_ipf;

	/// Time when last Bounds file printed
	float _bounds_time;
	/// FPS for Bounds files (-1 if Bounds file must not be printed)
	float _bounds_fps;
	/// Step when last Bounds file printed
	int _bounds_step;
	/// IPF for Bounds files (-1 if Bounds file must not be printed)
	int _bounds_ipf;

	/// Time when last Output file printed
	float _output_time;
	/// FPS for Output files (-1 if Output file must not be printed)
	float _output_fps;
	/// Step when last Output file printed
	int _output_step;
	/// IPF for Output files (-1 if Output file must not be printed)
	int _output_ipf;
};

}}  // namespace

#endif // TIMEMANAGER_H_INCLUDED
