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
 * @brief Complex tool to perform the link-list based on the "pos" array.
 * (See Aqua::CalcServer::LinkList for details)
 * @note Hardcoded versions of the files CalcServer/LinkList.cl.in and
 * CalcServer/LinkList.hcl.in are internally included as a text array.
 */

#ifndef LINKLIST_H_INCLUDED
#define LINKLIST_H_INCLUDED

#include <sphPrerequisites.h>
#include <vector>
#include <CalcServer/Tool.h>
#include <CalcServer/Reduction.h>
#include <CalcServer/RadixSort.h>

namespace Aqua{ namespace CalcServer{

/** @class LinkList LinkList.h CalcServer/LinkList.h
 * @brief Complex tool to perform the link-list based on the "pos" array. This
 * tool include the following steps:
 *   -# Minimum and maximum positions computations
 *   -# Number of cells calculation
 *   -# "ihoc" array allocation
 *   -# "ihoc" and "icell" calculations
 *   -# Radix sort of "icell", computing permutation array "id_sorted" and "id_unsorted" as well.
 * @note Hardcoded versions of the files CalcServer/LinkList.cl.in and
 * CalcServer/LinkList.hcl.in are internally included as a text array.
 */
class LinkList : public Aqua::CalcServer::Tool
{
public:
    /** Constructor.
     * @param tool_name Tool name.
     * @param input Input array to be used as the particles positions.
     * @param once Run this tool just once. Useful to make initializations.
     */
    LinkList(const std::string tool_name,
             const std::string input="pos",
             bool once=false);

    /** Destructor
     */
    ~LinkList();

    /** Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accesing the dependencies
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /** Setup the OpenCL stuff
     */
    void setupOpenCL();

    /** Compute the number of cells
     */
    void nCells();

    /** Allocate the "ihoc" array
     */
    void allocate();

    /** Update the input and output looking for changed values.
     */
    void setVariables();

    /// Input variable name
    std::string _input_name;

    /// Cells length
    float _cell_length;

    /// Number of cells
    uivec4 _n_cells;

    /// Minimum position computation tool
    Reduction *_min_pos;

    /// Maximum position computation tool
    Reduction *_max_pos;

    /// Sorting by cells computation tool
    RadixSort *_sort;

    /// "ihoc" array initialization
    cl_kernel _ihoc;
    /// "ihoc" array initialization local work size
    size_t _ihoc_lws;
    /// "ihoc" array initialization global work size
    size_t _ihoc_gws;
    /// "ihoc" array initialization sent arguments
    std::vector<void*> _ihoc_args;

    /// "icell" array computation
    cl_kernel _icell;
    /// "icell" array computation local work size
    size_t _icell_lws;
    /// "icell" array computation global work size
    size_t _icell_gws;
    /// "icell" array computation sent arguments
    std::vector<void*> _icell_args;

    /// "ihoc" array computation
    cl_kernel _ll;
    /// "ihoc" array computation local work size
    size_t _ll_lws;
    /// "ihoc" array computation global work size
    size_t _ll_gws;
    /// "ihoc" array computation sent arguments
    std::vector<void*> _ll_args;
};

}}  // namespace

#endif // LINKLIST_H_INCLUDED
