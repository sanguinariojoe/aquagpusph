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

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{

/** @class Reduction Reduction.h CalcServer/Reduction.h
 * @brief Variables array reduction. \n
 * The reduction is produced over other.
 */
class Reduction : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 * @param input Input data memory object.
	 * @param N Input data elements.
	 * @param type The data type array.
	 * @param identity The type considered the identity, i.e.
	 * Infinity for min operation, (float2)(0.f,0.f) for 2D vec
	 * sum reduction, etc.
	 * @param operation The reduction operation. For instance
	 * "a += b;" or
	 * "a.x = (a.x < b.x) ? a.x : b.x; a.y = (a.y < b.y) ? a.y : b.y;".
	 */
	Reduction(cl_mem input, unsigned int N, const char* type, const char* identity, const char* operation);

	/** Destructor.
	 */
	~Reduction();

	/** Reduction boundary computation.
	 * @return Output memory object, \n NULL if error is detected.
	 */
	cl_mem execute();

    /** Number of steps needed
     * @return Numkber of steps needed.
     */
    unsigned int nSteps(){return mGSize.size();}

    /** Return the memory object that stores the result
     * @return Memory object.
     */
    cl_mem resultMem(){return mMems.at(mMems.size() - 1);}

    /** Change the data array input.
     * @param New input data array.
     * @return false if all gone right, true otherwise.
     * @warning The new data array must be of the same size and type of the
     * used one in the construction. Otherwise, please destroy this object
     * and call the constructor again.
     */
    bool setInput(cl_mem input);
private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL(const char* type, const char* identity, const char* operation);

	/// OpenCL script path
	char* _path;

	/// OpenCL program
	cl_program _program;
	/// OpenCL kernel
	std::deque<cl_kernel> kernels;

    /// Global work sizes in each step
    std::deque<size_t> mGSize;
    /// Local work sizes in each step
    std::deque<size_t> mLSize;
    /// Number of work groups in each step
    std::deque<size_t> mNGroups;
    /// Number of input elements for each step
    std::deque<size_t> mN;

    /// Memory objects
    std::deque<cl_mem> mMems;
    /// Input memory object
    cl_mem mInput;
};

}}  // namespace

#endif // REDUCTION_H_INCLUDED
