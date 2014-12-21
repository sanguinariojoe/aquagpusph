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

#ifndef COPY_H_INCLUDED
#define COPY_H_INCLUDED

#include <CalcServer.h>
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
     */
    Copy(const char *name, const char *input_name, const char *output_name);

    /** Destructor.
     */
    ~Copy();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

protected:
    /** Compute the reduction.
     * @return false if all gone right, true otherwise.
     */
    bool _execute();

private:
    /** Get the input and output variables
     * @return false if all gone right, true otherwise
     */
    bool variables();

    /// Input variable name
    char* _input_name;
    /// Output variable name
    char* _output_name;

    /// Input variable
    InputOutput::ArrayVariable *_input_var;
    /// Output variable
    InputOutput::ArrayVariable *_output_var;
};

}}  // namespace

#endif // COPY_H_INCLUDED
