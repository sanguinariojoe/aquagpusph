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

#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

#include <sphPrerequisites.h>

#include <deque>
#include <CL/cl.h>
#include <clang-c/Index.h>
#include <clang-c/Platform.h>

#include <CalcServer/Tool.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{

/** @class Kernel Kernel.h CalcServer/Kernel.h
 * @brief A tool consisting in an OpenCL kernel execution. The variables used
 * in the OpenCL kernel are automatically detected.
 */
class Kernel : public Aqua::CalcServer::Tool
{
public:
    /** Constructor.
     * @param tool_name Tool name.
     * @param kernel_path Kernel path.
     * @param Number of threads to launch.
     */
    Kernel(const char* tool_name, const char* kernel_path, const char* n="N");

    /** Destructor
     */
    virtual ~Kernel();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

    /** Set the kernel file path.
     * @param kernel_path kernel file path.
     */
    void path(const char* kernel_path);

    /** Get the kernel file path.
     * @return Tool kernel file path.
     */
    const char* path(){return (const char*)_path;}

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
    bool _execute();

protected:
    /** Compile the OpenCL program
     * @param entry_point Program entry point method.
     * @param flags Compiling additional flags.
     * @param header Header to be append at the start of the source code.
     * @return false if all gone right, true otherwise.
     */
    bool compile(const char* entry_point="main",
                 const char* flags="",
                 const char* header="");

    /** Compute the variables required by the program
     * @param entry_point Program entry point method.
     * @return false if all gone right, true otherwise.
     */
    bool variables(const char* entry_point="main");

    /** Set the variables to the OpenCL kernel. The method detects if a variable
     * should be updated or if it already set either.
     * @return false if all gone right, true otherwise.
     */
    bool setVariables();

    /** Compute the global work size
     * @return false if all gone right, true otherwise.
     */
    bool computeGlobalWorkSize();

private:
    /// Kernel path
    char* _path;

    /// OpenCL kernel
    cl_kernel _kernel;

    /// work group size
    size_t _work_group_size;

    /// global work size
    size_t _global_work_size;

    /// Number of threads expression
    char* _n;

    /// List of required variables
    std::deque<char*> _var_names;
    /// List of variable values
    std::deque<void*> _var_values;
};

}}  // namespace

#endif // KERNEL_H_INCLUDED
