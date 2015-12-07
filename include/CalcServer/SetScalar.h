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
    /** @brief Constructor.
     * @param name Tool name.
     * @param var_name Variable to set.
     * @param value Value to set.
     */
    SetScalar(const char *name, const char *var_name, const char *value);

    /// Destructor.
    ~SetScalar();

    /** @brief Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

protected:
    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool _execute();

private:
    /** @brief Get the input variable
     * @return false if all gone right, true otherwise
     */
    bool variable();

    /// Input variable name
    char* _var_name;
    /// Value to set
    char* _value;

    /// Input variable
    InputOutput::Variable *_var;
};

}}  // namespace

#endif // SETSCALAR_H_INCLUDED
