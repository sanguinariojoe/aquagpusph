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
 * @brief Tools virtual environment to allow the user to define/manipulate the
 * tools used to carry out the simulation.
 * (see Aqua::CalcServer::Tool)
 */

#ifndef TOOL_H_INCLUDED
#define TOOL_H_INCLUDED

#include <sphPrerequisites.h>
#include <Variable.h>
#include <math.h>
#include <vector>

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
    /** Constructor.
     * @param tool_name Name of the tool. Useful to identify errors.
     * @param once Run this tool just once. Useful to make initializations.
     */
    Tool(const std::string tool_name, bool once=false);

    /** Destructor
     */
    virtual ~Tool();

    /** Set the tool name.
     * @param tool_name Tool name.
     */
    void name(const std::string tool_name){_name = tool_name;};

    /** Get the tool name.
     * @return Tool name.
     */
    const std::string name(){return _name;}

    /** Initialize the tool.
     */
    virtual void setup(){return;}

    /** @brief Execute the tool measuring the elapsed time.
     *
     * Actually this method is just ensuring that the tool can be executed,
     * e.g. the tool has been already executed, but it is asked to be ran just
     * once.
     * If the tool can be executed, then _execute() method is called, measuring
     * the time required to carry out the task.
     * @return false if all gone right, true otherwise.
     * @note Usually you don't want to overload this method, but the _execute()
     * protected one.
     */
    virtual void execute();

    /** Get the allocated memory for this tool.
     * @return allocated memory by this tool.
     */
    size_t allocatedMemory() const {return _allocated_memory;}

    /** Get the number of times that this tool has been called.
     * @return Number of times this tool has been called.
     */
    unsigned int used_times() const {return _n_iters;}

    /** Get the time consumed by the tool.
     * @param averaged true if the avergaed time step is required, false
     * otherwise.
     * @return time consumed.
     */
    float elapsedTime(bool averaged=true) const {
        if(!averaged)
            return _elapsed_time;
        return _average_elapsed_time;
    }

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
    /** Set the allocated memory for this tool.
     * @param mem_size allocated memory by this tool.
     */
    void allocatedMemory(size_t mem_size){_allocated_memory = mem_size;}

    /** Execute the tool.
     */
    virtual void _execute(){return;}

    /** @brief Add new data to the average and squared elapsed times
     * @param elapsed_time Elapsed time
     */
    void addElapsedTime(float elapsed_time);

    /** @brief Set the depedencies of the tool
     *
     * The dependencies are the variables that this tool is either reading or
     * writing.
     *
     * @param var_names Names of the dependencies
     */
    void setDependencies(std::vector<std::string> var_names);

    /** @brief Set the depedencies of the tool
     *
     * The dependencies are the variables that this tool is either reading or
     * writing.
     *
     * @param vars Dependencies
     */
    void setDependencies(std::vector<InputOutput::Variable*> vars);

    /** @brief Get the depedencies of the tool
     *
     * @return Dependencies
     */
    const std::vector<InputOutput::Variable*> getDependencies();

    /** @brief Get the list of events that this tool shall wait for
     *
     * @return C++ vector of events
     * @warning The events returned have been retained, so call clReleaseEvent()
     * after using them
     */
    const std::vector<cl_event> getEvents();

private:
    /// Kernel name
    std::string _name;

    /// true if the tool shall be run just once, false otherwise
    bool _once;

    /// Total auxiliar memory allocated in the device
    size_t _allocated_memory;

    /// Times that this tool has been called
    unsigned int _n_iters;

    /// Average elapsed time
    float _elapsed_time;

    /// Average elapsed time
    float _average_elapsed_time;

    /// Average squared elapsed time
    float _squared_elapsed_time;

    /// List of dependencies
    std::vector<InputOutput::Variable*> _vars;

    /// Internal storage to can safely return memory with getEvents()
    std::vector<cl_event> _events;
};

}}  // namespace

#endif // TOOL_H_INCLUDED
