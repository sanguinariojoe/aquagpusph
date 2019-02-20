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
 * @brief OpenCL kernel kernel based tool.
 * (see Aqua::CalcServer::Kernel for details)
 */

#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

#include <sphPrerequisites.h>

#include <vector>
#include <CL/cl.h>
#include <CalcServer/Tool.h>

namespace Aqua{ namespace CalcServer{

/** @class Kernel Kernel.h CalcServer/Kernel.h
 * @brief A tool consisting in an OpenCL kernel execution. The variables used
 * in the OpenCL kernel are automatically detected.
 */
class Kernel : public Aqua::CalcServer::Tool
{
private:
    /** @class event_wait_list Kernel.h CalcServer/Kernel.h
     * @brief Helper to store events wait list
     *
     * This is useful to can send the data to the event listener.
     */
    class event_wait_list
    {
    public:
        /** Constructor
         * @param events Events to be waited.
         */
        event_wait_list(std::vector<cl_event> events);

        /** Destructor
         */
        ~event_wait_list();
    
        /// Number of events
        cl_uint num;

        /// List of events
        cl_event *list;
    };


public:
    /** Constructor.
     * @param tool_name Tool name.
     * @param kernel_path Kernel path.
     * @param n Number of threads to launch.
     * @param once Run this tool just once. Useful to make initializations.
     */
    Kernel(const std::string tool_name,
           const std::string kernel_path,
           const std::string entry_point="entry",
           const std::string n="N",
           bool once=false);

    /** Destructor
     */
    virtual ~Kernel();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    void setup();

    /** Get the kernel file path.
     * @return Tool kernel file path.
     */
    const std::string path(){return (const std::string)_path;}

    /** Get the work group size
     * @return Work group size
     */
    size_t workGroupSize() const {return _work_group_size;}

    /** Get the work group size
     * @return Work group size
     */
    size_t globalWorkSize() const {return _global_work_size;}

protected:
    /** Execute the tool.
     * @return false if all gone right, true otherwise.
     */
    void _execute();

protected:
    /** Compile the OpenCL program
     * @param entry_point Program entry point method.
     * @param flags Compiling additional flags.
     * @param header Header to be append at the start of the source code.
     * @return false if all gone right, true otherwise.
     */
    void compile(const std::string entry_point="entry",
                 const std::string flags="",
                 const std::string header="");

    /** Compute the variables required by the program
     * @param entry_point Program entry point method.
     * @return false if all gone right, true otherwise.
     */
    void variables(const std::string entry_point="main");

    /** @brief Set the variables to the OpenCL kernel.
     * 
     * The method detects if a variable should be updated or if it already set either.
     */
    void setVariables();

    /** Compute the global work size
     */
    void computeGlobalWorkSize();

private:
    /// Kernel path
    std::string _path;

    /// Kernel entry point
    std::string _entry_point;

    /// Number of threads expression
    std::string _n;

    /// OpenCL kernel
    cl_kernel _kernel;

    /// work group size
    size_t _work_group_size;

    /// global work size
    size_t _global_work_size;

    /// List of required variables
    std::vector<std::string> _var_names;
    /// List of variable values
    std::vector<void*> _var_values;
};

}}  // namespace

#endif // KERNEL_H_INCLUDED
