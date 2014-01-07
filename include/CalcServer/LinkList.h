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

#ifndef LINKLIST_H_INCLUDED
#define LINKLIST_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include Radix sort tool
// ----------------------------------------------------------------------------
#include <CalcServer/RadixSort.h>

namespace Aqua{ namespace CalcServer{

/** @class LinkList LinkList.h CalcServer/LinkList.h
 * @brief In cell particles localization. When particle cells have been determined
 * a chain of them is set. To learn more see Grid.
 */
class LinkList : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	LinkList();

	/** Destructor.
	 */
	~LinkList();

	/** Performs link-list.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();
	/** Auxiliar method to alloc memory for link-list (when needed).
	 * @return false if all gone right. \n true otherwise.
	 */
	bool allocLinkList();

	/// OpenCL script path
	char *_path;

	/// OpenCL program
	cl_program program;
	/// OpenCL icell kernel
	cl_kernel clLcellKernel;
	/// OpenCL ihoc init kernel
	cl_kernel clIhocKernel;
	/// OpenCL link-list kernel
	cl_kernel clLLKernel;
	/// Global work size.
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;

	/// Radix sort module
	RadixSort *mRadixSort;
};

}}  // namespace

#endif // LINKLIST_H_INCLUDED
