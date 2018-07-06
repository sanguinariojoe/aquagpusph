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
 * @note Hardcoded versions of the files CalcServer/Reduction.cl.in and
 * CalcServer/Reduction.hcl.in are internally included as a text array.
 */

#ifndef REDUCTION_H_INCLUDED
#define REDUCTION_H_INCLUDED

#include <vector>

#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Reduction Reduction.h CalcServer/Reduction.h
 * @brief Reductions, like scans, prefix sums, maximum or minimum, etc...
 * @see Reduction.cl
 * @note Hardcoded versions of the files CalcServer/Reduction.cl.in and
 * CalcServer/Reduction.hcl.in are internally included as a text array.
 */
class Reduction : public Aqua::CalcServer::Tool
{
public:
    /** @brief Reduction definition.
     * @param name Tool name.
     * @param input_name Variable to be reduced name.
     * @param output_name Variable where the reduced value will be stored.
     * @param operation The reduction operation.
     * For instance:
     *   - "c += b;"
     *   - "c.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;"
     * @param null_val The value considered as the null one, i.e. INFINITY for
     * float min value reduction, or (vec2)(0.f,0.f) for a 2D vec sum reduction.
     * @note Some helpers are available for null_val:
     *   - VEC_ZERO: Zeroes vector.
     *   - VEC_ONE: Ones vector, in 3D cases the last component will be zero.
     *   - VEC_ALL_ONE: Equal to VEC_ONE, but in 3D cases the last component will be one as well.
     *   - VEC_INFINITY: INFINITY components vector, in 3D cases the last component will be zero.
     *   - VEC_ALL_INFINITY: Equal to VEC_INFINITY, but in 3D cases the last component will be INFINITY as well.
     *   - VEC_NEG_INFINITY: -VEC_INFINITY
     *   - VEC_ALL_NEG_INFINITY: -VEC_ALL_INFINITY.
     */
    Reduction(const std::string name,
              const std::string input_name,
              const std::string output_name,
              const std::string operation,
              const std::string null_val);

    /// Destructor.
    ~Reduction();

    /** @brief Initialize the tool.
     *
     * This method should be called after the constructor, such that it could
     * report errors that the application may handle quitting in a safe way.
     */
    void setup();

    /** @brief Number of steps needed.
     *
     * To reduce the array to just one variable several steps may be needed,
     * depending on the number of work groups that should be launched at each
     * pass.
     *
     * @return Number of steps needed.
     */
    unsigned int nSteps(){return _global_work_sizes.size();}

protected:
    /** @brief Perform the work.
     * @return Output memory object, NULL if error is detected.
     */
    void _execute();

private:
    /** @brief Extract the input and output variables from the provided data in
     * Reduction().
     * @see Aqua::InputOutput::Variables
     */
    void variables();

    /** @brief Setup the OpenCL stuff
     */
    void setupOpenCL();

    /** @brief Compile the source code and generate the corresponding kernel.
     * @param source Source code to be compiled.
     * @param local_work_size Desired local work size.
     * @return Kernel instance.
     */
    cl_kernel compile(const std::string source, size_t local_work_size);

    /** Update the input variables.
     *
     * This function is looking for changed value to send them again to the
     * computational device.
     */
    void setVariables();

    /// Input variable name
    std::string _input_name;
    /// Output variable name
    std::string _output_name;
    /// Operation to be computed
    std::string _operation;
    /// Considered null val
    std::string _null_val;

    /// Input variable
    InputOutput::ArrayVariable *_input_var;
    /// Output variable
    InputOutput::Variable *_output_var;

    /// Input array
    cl_mem _input;

    /// OpenCL kernels
    std::vector<cl_kernel> _kernels;

    /// Global work sizes in each step
    std::vector<size_t> _global_work_sizes;
    /// Local work sizes in each step
    std::vector<size_t> _local_work_sizes;
    /// Number of work groups in each step
    std::vector<size_t> _number_groups;
    /// Number of input elements for each step
    std::vector<size_t> _n;

    /// Memory objects
    std::vector<cl_mem> _mems;
};

}}  // namespace

#endif // REDUCTION_H_INCLUDED
