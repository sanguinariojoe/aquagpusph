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
 * @brief Array copy tool.
 * (see Aqua::CalcServer::Copy for details)
 */

#ifndef COPY_H_INCLUDED
#define COPY_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Copy Copy.h CalcServer/Copy.h
 * @brief Copy an array component by component.
 */
class Copy : public Aqua::CalcServer::Tool
{
public:
    /** Constructor.
     * @param name Tool name.
     * @param input_name Variable to copy.
     * @param output_name Variable to set.
     * @param once Run this tool just once. Useful to make initializations.
     */
    Copy(const std::string name,
         const std::string input_name,
         const std::string output_name,
         bool once=false);

    /** Destructor.
     */
    ~Copy();

    /** Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accessing the dependencies
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /** Get the input and output variables
     */
    void variables();

    /// Input variable name
    std::string _input_name;
    /// Output variable name
    std::string _output_name;

    /// Input variable
    InputOutput::ArrayVariable *_input_var;
    /// Output variable
    InputOutput::ArrayVariable *_output_var;
};

}}  // namespace

#endif // COPY_H_INCLUDED
