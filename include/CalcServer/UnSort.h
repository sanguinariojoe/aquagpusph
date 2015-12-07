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
 * @brief UnSort Recover the original id of each particle.
 * (See Aqua::CalcServer::UnSort for details)
 * @note Hardcoded versions of the files CalcServer/UnSort.cl.in and
 * CalcServer/UnSort.hcl.in are internally included as a text array.
 */

#ifndef UNSORT_H_INCLUDED
#define UNSORT_H_INCLUDED

#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class UnSort UnSort.h CalcServer/UnSort.h
 * @brief UnSort Recover the original id of each particle. This tool is not
 * designed for the common usage but as an auxiliar tool for the savers,
 * therefore it will not be selectable for the users.
 */
class UnSort : public Aqua::CalcServer::Tool
{
public:
    /** Constructor.
     * @param name Tool name.
     * @param var_name Variable to unsort.
     */
    UnSort(const char *name, const char *var_name);

    /** Destructor.
     */
    ~UnSort();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

    /** Get the memory object where the unsorted data is stored.
     * @return The memory object where the unsorted data is stored.
     */
    cl_mem output(){return _output;}

protected:
    /** Compute the reduction.
     * @return false if all gone right, true otherwise.
     */
    bool _execute();

private:
    /** Get the input variable
     * @return false if all gone right, true otherwise
     */
    bool variables();

    /** Create the output memory object
     * @return false if all gone right, true otherwise.
     */
    bool setupMem();

    /** Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /** Compile the source code and generate the corresponding kernel
     * @param source Source code to be compiled.
     * @return Kernel instance, NULL if error happened.
     */
    cl_kernel compile(const char* source);

    /** Update the input looking for changed value.
     * @return false if all gone right, true otherwise.
     */
    bool setVariables();

    /// Input variable name
    char* _var_name;

    /// ID variable
    InputOutput::ArrayVariable *_id_var;

    /// Input variable
    InputOutput::ArrayVariable *_var;

    /// ID Memory object sent
    cl_mem _id_input;

    /// Memory object sent
    cl_mem _input;

    /// Output memory object
    cl_mem _output;

    /// OpenCL kernel
    cl_kernel _kernel;

    /// Global work sizes in each step
    size_t _global_work_size;
    /// Local work sizes in each step
    size_t _local_work_size;
    /// Number of elements
    unsigned int _n;
};

}}  // namespace

#endif // UNSORT_H_INCLUDED
