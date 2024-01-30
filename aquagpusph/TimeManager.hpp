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

#include "sphPrerequisites.hpp"
#include "ProblemSetup.hpp"

namespace Aqua {
namespace InputOutput {

/** @class TimeManager TimeManager.h TimeManager.h
 * @brief Simulation time flow manager.
 *
 * This class controls the time variable \f$ t \f$, and more specifically the
 * time events triggers, like the output updates, the number of output frames
 * performed, etc...
 *
 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
 */
struct TimeManager
{
  public:
	/** Constructor
	 *
	 * @param sim_data Simulation data
	 */
	TimeManager(const ProblemSetup& sim_data);

	/// Destructor
	~TimeManager();

	/** @brief Check if the simulation must be finished.
	 * @return true if simulation should finish, false otherwise.
	 */
	bool mustStop();

	/** @brief Check if a general simulation output must be printed.
	 * @return true if an output should be printed, false otherwise.
	 * @warning This method is returning true just one time per time step (i.e.
	 * until update() is called again).
	 */
	bool mustPrintOutput();

	/** @brief Get the simulation time step index.
	 * @return Simulation time step index.
	 */
	inline unsigned int step() const { return *_step; }

	/** @brief Get the simulation time instant.
	 * @return Simulation time instant.
	 */
	inline float time() const { return *_time; }

	/** @brief Get the simulation frame.
	 *
	 * The frame is the index of the current particles output.
	 *
	 * @return Simulation frame.
	 */
	inline unsigned int frame() const { return *_frame; }

	/** @brief Get the simulation time step \f$ \Delta t \f$.
	 * @return Simulation time step \f$ \Delta t \f$.
	 */
	inline float dt() const { return *_dt; }

	/** @brief Get the total simulation time to compute.
	 * @return Total simulation time to compute.
	 */
	inline float maxTime() const { return *_time_max; }

	/** @brief Get the number of frames to compute.
	 * @return Number of frames to compute.
	 */
	inline unsigned int maxStep() const { return *_steps_max; }

	/** @brief Get the number of frames to compute.
	 * @return Number of frames to compute.
	 */
	inline unsigned int maxFrame() const { return *_frames_max; }

	/** @brief Set the last output event time step index.
	 * @param s last output event time step index.
	 */
	inline void outputStep(int s) { _output_step = s; }

	/** @brief Get the last output event time step index.
	 * @return Last output event time step index.
	 */
	inline int outputStep() const { return _output_step; }

	/** @brief Get the iterations per output frame.
	 * @return Iterations per frame.
	 */
	inline int outputIPF() const { return _output_ipf; }

	/** @brief Set the last output event time instant.
	 * @param t Last output event time instant.
	 */
	inline void outputTime(float t) { _output_time = t; }

	/** @brief Get the last output event time instant.
	 * @return Last output event time instant.
	 */
	inline float outputTime() const { return _output_time; }

	/** @brief Get the output frames per second.
	 * @return Frames per second.
	 */
	inline float outputFPS() const { return _output_fps; }

  private:
	/** @brief Ask to synchronize all the dependant variables, so they can be
	 * safety read and wrote
	 */
	void sync();

	/// Actual step
	unsigned int* _step;
	/// Actual time
	float* _time;
	/// Time step
	float* _dt;
	/// Actual frame
	unsigned int* _frame;
	/// Maximum time into simulation (-1 if simulation don't stop by time
	/// criteria)
	float* _time_max;
	/// Maximum number of steps into simulation (-1 if simulation don't stop by
	/// steps criteria)
	unsigned int* _steps_max;
	/// Maximum number of frames into simulation (-1 if simulation don't stop by
	/// frames criteria)
	unsigned int* _frames_max;

	/// Time when last Output file printed
	float _output_time;
	/// FPS for Output files (-1 if Output file must not be printed)
	float _output_fps;
	/// Step when last Output file printed
	int _output_step;
	/// IPF for Output files (-1 if Output file must not be printed)
	int _output_ipf;
};

}
} // namespace

#endif // TIMEMANAGER_H_INCLUDED
