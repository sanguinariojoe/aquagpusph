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
 * @brief Set a scalar variable.
 * (See Aqua::CalcServer::SetScalar for details)
 */

#ifndef SETSCALAR_H_INCLUDED
#define SETSCALAR_H_INCLUDED

#include <CalcServer.h>
#include <CalcServer/Tool.h>

namespace Aqua{ namespace CalcServer{

/** @class SetScalar SetScalar.h CalcServer/SetScalar.h
 * @brief Set a scalar variable.
 */
class SetScalar : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor
     * @param name Tool name.
     * @param var_name Variable to set.
     * @param value Value to set.
     * @param once Run this tool just once. Useful to make initializations.
     */
    SetScalar(const std::string& name,
              const std::string& var_name,
              const std::string& value,
              const bool once=false);

    /// Destructor
    ~SetScalar();

    /** @brief Initialize the tool
     */
    void setup();

protected:
    /** @brief Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accessing the dependencies
     */
    const cl_event _execute(const std::vector<cl_event>& events);

private:
    /** @brief Get the input variable
     */
    void variable();

    /// Input variable name
    std::string _var_name;
    /// Value to set
    std::string _value;

    /// Input variable
    InputOutput::Variable *_var;
};

}}  // namespace

#endif // SETSCALAR_H_INCLUDED
