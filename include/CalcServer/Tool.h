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

#ifndef TOOL_H_INCLUDED
#define TOOL_H_INCLUDED

#include <sphPrerequisites.h>
#include <math.h>

namespace Aqua{ namespace CalcServer{

/** @class Tool Tool.h CalcServer/Tool.h
 * @brief Tools base class. The way that AQUAgpusph compute each problem is set
 * through a set of tools that are computed sequentially. Several tools can be
 * considered, for instance:
 *   -# A single OpenCL kernel
 *   -# A more complex OpenCL tool like Reductions or LinkList
 *   -# Python scripts
 *   -# Variables set
 */
class Tool
{
public:
    /** Destructor
     */
    virtual ~Tool();

    /** Set the tool name.
     * @param tool_name Tool name.
     */
    void name(const char* tool_name);

    /** Get the tool name.
     * @return Tool name.
     */
    const char* name(){return (const char*)_name;}

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    virtual bool setup(){return false;}

    /** Execute the tool measuring the elapsed time.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

    /** Get the allocated memory for this tool.
     * @return allocated memory by this tool.
     */
    size_t allocatedMemory() const {return _allocated_memory;}

    /** Get the number of times that this tool has been called.
     * @return Number of times this tool has been called.
     */
    unsigned int used_times() const {return _n_iters;}

    /** Get the average time consumed by the tool.
     * @return Average time consumed.
     */
    float elapsedTime() const {return _average_elapsed_time;}

    /** Get the time consumed variance.
     * @return Time consumed variance.
     */
    float elapsedTimeVariance() const {
        return _squared_elapsed_time - pow(_average_elapsed_time, 2);
    }

    /** Get the time consumed standard deviation.
     * @return Time consumed standard deviation.
     */
    float elapsedTimeDeviation() const {return sqrt(elapsedTimeVariance());}

protected:
    /** Constructor.
     * @param tool_name Name of the tool. Useful to identify errors.
     */
    Tool(const char* tool_name);

    /** Set the allocated memory for this tool.
     * @param mem_size allocated memory by this tool.
     */
    void allocatedMemory(size_t mem_size){_allocated_memory = mem_size;}

    /** Execute the tool.
     * @return false if all gone right, true otherwise.
     */
    virtual bool _execute(){return false;}

private:
    /** @brief Add new data to the average and squared elapsed times computation.
     * @param elapsed_time Elapsed time of the last time that this tool has been called.
     */
    void addElapsedTime(float elapsed_time);

    /// Kernel name
    char* _name;

    /// Total auxiliar memory allocated in the device
    size_t _allocated_memory;

    /// Times that this tool has been called
    unsigned int _n_iters;

    /// Average elapsed time
    float _average_elapsed_time;

    /// Average squared elapsed time
    float _squared_elapsed_time;
};

}}  // namespace

#endif // TOOL_H_INCLUDED
