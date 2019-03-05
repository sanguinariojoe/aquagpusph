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
#include <ProblemSetup.h>
#include <Variable.h>

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

    /** @brief Check if the simulation must be finished
     * @return true if simulation is over, false otherwise
     */
    const bool mustStop() const;

    /** @brief Check if simulation output must be printed
     * @return true if an output should be printed, false otherwise
     * @warning This method is returning true just once time per time step (i.e.
     * further calls during the same time step is returning always false)
     */
    const bool mustPrintOutput();

    /** @brief Get the simulation time step index.
     * @return Simulation time step index.
     */
    inline const unsigned int step() const {return *_step;}

    /** @brief Get the simulation time instant
     * @return Simulation time instant
     */
    inline const float time() const {return *_time;}

    /** @brief Get the simulation frame
     * @return Simulation frame.
     * @note Frame stands for output files, while step stands for the number of
     * time iterations
     */
    inline const unsigned int frame() const {return *_frame;}

    /** @brief Get the simulation time step \f$ \Delta t \f$
     * @return Simulation time step \f$ \Delta t \f$
     */
    inline const float dt() const {return *_dt;}

    /** @brief Get the last output event time step index
     * @return Last output event time step index
     */
    inline const int outputStep() const {return _output_step;}

    /** @brief Get the last output event time instant
     * @return Last output event time instant
     */
    inline const float outputTime() const {return _output_time;}

    /** @brief Get the output frames per second
     * @return Frames per second
     */
    inline const float outputFPS() const {return _output_fps;}

    /** @brief Get the total simulation time to compute
     * @return Total simulation time to compute
     */
    inline const float maxTime() const {return *_time_max;}

    /** @brief Get the number of frames to compute
     * @return Number of frames to compute
     */
    inline const unsigned int maxStep() const {return *_steps_max;}

    /** @brief Get the number of frames to compute
     * @return Number of frames to compute
     */
    inline const unsigned int maxFrame() const {return *_frames_max;}

private:
    /** @brief Get the simulation frame
     * @return Simulation frame.
     * @note Frame stands for output files, while step stands for the number of
     * time iterations
     */
    void frame(const unsigned int& frame);

    /// Current time step
    unsigned int *_step;
    /// Current time
    float *_time;
    /// Time step
    float *_dt;
    /// Current frame (frame stands for the number of outputs)
    unsigned int *_frame;
    /// Current frame variable
    FloatVariable *_frame_var;
    /// Maximum simulation time (-1 if simulation don't stop by time criteria)
    float *_time_max;
    /// Maximum number of steps (-1 if simulation don't stop by steps criteria)
    unsigned int *_steps_max;
    /// Maximum number of frames (-1 if simulation don't stop by frames criteria)
    unsigned int *_frames_max;

    /// Time when last Output file was printed
    float _output_time;
    /// FPS for Output files (-1 to neglect this criteria)
    float _output_fps;
    /// Step when last Output file was printed
    int _output_step;
    /// IPF for Output files (-1 if neglect this criteria)
    int _output_ipf;
};

}}  // namespace

#endif // TIMEMANAGER_H_INCLUDED
