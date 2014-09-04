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
 * @brief Link-list computation tool.
 * (See Aqua::CalcServer::LinkList for details)
 */

#ifndef LINKLIST_H_INCLUDED
#define LINKLIST_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/RadixSort.h>

namespace Aqua{ namespace CalcServer{

/** @class LinkList LinkList.h CalcServer/LinkList.h
 * @brief Link-List computation tool.
 *
 * Link-list is a PIC (Particle In Cell) based neighs localization.
 *
 * In this technique the particles are located in a lattice grid such that it
 * can be asserted that the non null interactions can only happens between
 * particles in neigh cells.
 *
 * After that the particles are sorted by its cell indexes, such that particles
 * in the same cell can be found just increasing/decreasing the particle index.
 *
 * Finally a "Head of Chain" array is generated to identify the first particle
 * in each cell.
 *
 * Such information allows to get all the particles in a cell just taking the
 * "Head of Chain", and increasing the index while a particle out of the cell
 * is not detected.
 *
 * Please visit the user manual, chapter 7.5 to learn more about the link-list.
 *
 * @see LinkList.cl
 * @see Aqua::CalcServer::Grid
 * @see Aqua::CalcServer::RadixSort
 */
class LinkList : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    LinkList();

    /// Destructor.
    ~LinkList();

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** @brief Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /** @brief Helper method to allocate memory for the link-list (when needed).
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
