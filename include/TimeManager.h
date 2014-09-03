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

#ifndef TIMEMANAGER_H_INCLUDED
#define TIMEMANAGER_H_INCLUDED

#include <sphPrerequisites.h>
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class TimeManager TimeManager.h TimeManager.h
 * @brief Simulation time flow manager.
 *
 * This class controls the time variable \f$ t \f$, and more specifically the
 * time events triggers, like the output updates, the number of output frames
 * performed, etc...
 *
 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
 */
struct TimeManager : public Aqua::Singleton<Aqua::InputOutput::TimeManager>
{
public:
	/// Constructor
	TimeManager();

	/// Destructor
	~TimeManager();

	/** @brief Pass to the next time step.
	 * @param dt Simulation time elapsed in this frame.
	 */
	void update(float dt);

	/** @brief Check if the simulation must be finished.
	 * @return true if simulation should finish, false otherwise.
	 */
	bool mustStop();

	/** @brief Check if the log file must be updated.
	 * @return true if the log should be updated, false otherwise.
	 * @warning This method is returning true just one time per time step (i.e.
     * until update() is called again).
	 */
	bool mustPrintLog();

	/** @brief Check if the energy file must be updated.
	 * @return true if the energy should be updated, false otherwise.
	 * @warning This method is returning true just one time per time step (i.e.
     * until update() is called again).
	 */
	bool mustPrintEnergy();

	/** @brief Check if the fluid bounds file must be updated.
	 * @return true if the bounds file should be updated, false otherwise.
	 * @warning This method is returning true just one time per time step (i.e.
     * until update() is called again).
	 */
	bool mustPrintBounds();

	/** @brief Check if a general simulation output must be printed.
	 * @return true if an output should be printed, false otherwise.
	 * @warning This method is returning true just one time per time step (i.e.
     * until update() is called again).
	 */
	bool mustPrintOutput();

	/** @brief Set the simulation time step index.
	 * @param s Simulation time step index.
	 */
	void step(unsigned int s){_step = s;}
	/** @brief Get the simulation time step index.
	 * @return Simulation time step index.
	 */
	unsigned int step(){return _step;}
	/** @brief Set the simulation time instant.
	 * @param t Simulation time instant.
	 */
	void time(float t){_time = t;}
	/** @brief Get the simulation time instant.
	 * @return Simulation time instant.
	 */
	float time(){return _time;}
	/** @brief Set the simulation frame.
	 *
	 * The frame is the index of the current particles output.
	 *
	 * @param frame Simulation frame.
	 */
	void frame(unsigned int frame){_frame = frame;}
	/** @brief Get the simulation frame.
	 *
	 * The frame is the index of the current particles output.
	 *
	 * @return Simulation frame.
	 */
	unsigned int frame(){return _frame;}
	/** @brief Set the simulation time step \f$ \Delta t \f$.
	 * @param dt Simulation time step \f$ \Delta t \f$.
	 */
	void dt(float dt){_dt = dt;}
	/** @brief Get the simulation time step \f$ \Delta t \f$.
	 * @return Simulation time step \f$ \Delta t \f$.
	 */
	float dt(){return _dt;}
	/** @brief Set the simulation starting frame.
	 *
	 * The frame is the index of the current particles output.
	 *
	 * @param start starting frame.
	 */
	void startFrame(int start){_start_frame=start;}
	/** @brief Get the simulation starting frame.
	 *
	 * The frame is the index of the current particles output.
	 *
	 * @return Simulation starting frame.
	 */
	int startFrame(){return _start_frame;}
	/** @brief Set the simulation starting time \f$ t_0 \f$.
	 * @param start Simulation starting time \f$ t_0 \f$.
	 */
	void startTime(float start){_start_time=start;}
	/** @brief Get the simulation starting time \f$ t_0 \f$.
	 * @return Simulation starting time \f$ t_0 \f$.
	 */
	float startTime(){return _start_time;}

	/** @brief Set the last output event time step index.
	 * @param s last output event time step index.
	 */
	void outputStep(int s){_output_step = s;}
	/** @brief Get the last output event time step index.
	 * @return Last output event time step index.
	 */
	int outputStep(){return _output_step;}
	/** @brief Get the iterations per output frame.
	 * @return Iterations per frame.
	 */
	int outputIPF(){return _output_ipf;}
	/** @brief Set the last output event time instant.
	 * @param t Last output event time instant.
	 */
	void outputTime(float t){_output_time = t;}
	/** @brief Get the last output event time instant.
	 * @return Last output event time instant.
	 */
	float outputTime(){return _output_time;}
	/** @brief Get the output frames per second.
	 * @return Frames per second.
	 */
	float outputFPS(){return _output_fps;}

	/** @brief Get the total simulation time to compute.
	 * @return Total simulation time to compute.
	 */
	float maxTime(){return _time_max;}
	/** @brief Get the number of frames to compute.
	 * @return Number of frames to compute.
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
