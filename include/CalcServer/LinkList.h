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

#include <CalcServer/Kernel.h>
#include <CalcServer/RadixSort.h>

namespace Aqua{ namespace CalcServer{

/** @class LinkList LinkList.h CalcServer/LinkList.h
 * @brief PIC (Particle In Cell) based neighbours localization. This method
 * will sort the particles in function of its cells, shuch that it can be
 * granted that the next particle in each cell will be the next one. This
 * process optimize the interactions stage with high performance increases.
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

	/** Performs the link-list.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	bool setupOpenCL();

	/** Auxiliar method to alloc memory for the link-list (when needed).
	 * @return false if all gone right, true otherwise.
	 */
	bool allocLinkList();

	/// OpenCL script path
	char *_path;

	/// OpenCL program
	cl_program _program;
	/// OpenCL icell kernel
	cl_kernel _icell_kernel;
	/// OpenCL ihoc init kernel
	cl_kernel _ihoc_kernel;
	/// OpenCL link-list kernel
	cl_kernel _ll_kernel;
	/// Global work size.
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;

	/// Radix sort module
	RadixSort *_radix_sort;
};

}}  // namespace

#endif // LINKLIST_H_INCLUDED
