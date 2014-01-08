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

#ifndef REDUCTION_H_INCLUDED
#define REDUCTION_H_INCLUDED

#include <deque>

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Reduction Reduction.h CalcServer/Reduction.h
 * @brief Array reduction tool. It could every prefix scan operation that you
 * want.
 */
class Reduction : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 * @param input Input data memory object.
	 * @param n Input data elements.
	 * @param type The data type array.
	 * @param null_val The value considered as the null one, i.e.
	 * Infinity for min operation, (float2)(0.f,0.f) for 2D vec sum
	 * reduction, etc.
	 * @param operation The reduction operation. For instance:
	 *   - "a += b;"
	 *   - "a.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;"
	 */
	Reduction(cl_mem input,
              unsigned int n,
              const char* type,
              const char* null_val,
              const char* operation);

	/** Destructor.
	 */
	~Reduction();

	/** Compute the prefix scan.
	 * @return Output memory object, NULL if error is detected.
	 */
	cl_mem execute();

    /** Number of steps needed
     * @return Number of steps needed.
     */
    unsigned int nSteps(){return _global_work_sizes.size();}

    /** Return the memory object which stores the final result
     * @return Memory object.
     */
    cl_mem resultMem(){return _mems.at(_mems.size() - 1);}

    /** Change the data array input.
     * @param New input data array.
     * @return false if all gone right, true otherwise.
     * @warning The new data array must be of the same size and type of the
     * previously used one in the construction. Otherwise, please destroy this
     * object and call the constructor again.
     */
    bool setInput(cl_mem input);

private:
	/** Setup the OpenCL stuff
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
