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
 * @brief Reductions, like scans, prefix sums, maximum or minimum, etc...
 * (See Aqua::CalcServer::Reduction for details)
 */

#ifndef REDUCTION_H_INCLUDED
#define REDUCTION_H_INCLUDED

#include <deque>

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Reduction Reduction.h CalcServer/Reduction.h
 * @brief Reductions, like scans, prefix sums, maximum or minimum, etc...
 * @see Reduction.cl
 */
class Reduction : public Aqua::CalcServer::Kernel
{
public:
    /** @brief Reduction definition.
     * @param input Input data memory object.
     * @param n Input data elements.
     * @param type The data type array.
     * @param null_val The value considered as the null one, i.e.
     * Infinity for min operation, (float2)(0.f,0.f) for 2D vec prefix sum, etc.
     * @param operation The reduction operation.
     * For instance:
     *   - "a += b;"
     *   - "a.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;"
     */
    Reduction(cl_mem input,
              unsigned int n,
              const char* type,
              const char* null_val,
              const char* operation);

    /// Destructor.
    ~Reduction();

    /** @brief Perform the work.
     * @return Output memory object, NULL if error is detected.
     */
    cl_mem execute();

    /** @brief Number of steps needed.
     *
     * To reduce the array to just one variable several steps may be needed,
     * depending on the number of work groups that should be launched at each
     * pass.
     *
     * @return Number of steps needed.
     */
    unsigned int nSteps(){return _global_work_sizes.size();}

    /** @brief Memory object where the final result is stored.
     * @return Memory object.
     */
    cl_mem resultMem(){return _mems.at(_mems.size() - 1);}

    /** @brief Memory buffer to be reduced.
     * @param input New input data array.
     * @return false if all gone right, true otherwise.
     * @warning The new data array must be of the same size and type of the
     * previously used one in the construction.
     * Otherwise, please destroy this object and call the constructor again.
     */
    bool setInput(cl_mem input);

private:
    /** @brief Setup the OpenCL stuff
     * @param type The data type array.
     * @param null_val The value considered as the null one, i.e.
     * Infinity for min operation, (float2)(0.f,0.f) for 2D vec sum
     * reduction, etc.
     * @param operation The reduction operation. For instance:
     *   - "a += b;"
     *   - "a.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;"
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL(const char* type, const char* null_val, const char* operation);

    /// OpenCL script path
    char* _path;

    /// OpenCL program
    cl_program _program;
    /// OpenCL kernel
    std::deque<cl_kernel> _kernels;

    /// Global work sizes in each step
    std::deque<size_t> _global_work_sizes;
    /// Local work sizes in each step
    std::deque<size_t> _local_work_sizes;
    /// Number of work groups in each step
    std::deque<size_t> _number_groups;
    /// Number of input elements for each step
    std::deque<size_t> _n;

    /// Memory objects
    std::deque<cl_mem> _mems;
    /// Input memory object
    cl_mem _input;
};

}}  // namespace

#endif // REDUCTION_H_INCLUDED
