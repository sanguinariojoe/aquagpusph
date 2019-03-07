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
 * @brief Check that a condition holds true, or trhow a fatal error otherwise.
 * (See Aqua::CalcServer::Assert for details)
 */

#ifndef ASSERT_H_INCLUDED
#define ASSERT_H_INCLUDED

#include <CalcServer/Tool.h>

namespace Aqua{ namespace CalcServer{

/** @class Assert Assert.h CalcServer/Assert.h
 * @brief Check that a condition holds true, or trhow a fatal error otherwise.
 *
 * If the result of evaluating the condition expression is equal to 0, the
 * result will be considered false, and therefore a fatal error will be raised.
 * Any other value will be considered as true, letting the simulation to
 * normally continue.
 */
class Assert : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor
     * @param name Tool name.
     * @param condition Condition to evaluate. If the result is 0, false will be
     * considered, and fatal error will be raised, otherwise the simulation will
     * continue.
     * @param once Run this tool just once. Useful to make initializations.
     */
    Assert(const std::string& name,
           const std::string& condition,
           const bool once=false);

    /// Destructor
    ~Assert();

    /** @brief Initialize the tool.
     */
    void setup();

protected:
    /** @brief Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accessing the dependencies
     */
    const cl_event _execute(const std::vector<cl_event>& events);

private:
    /// Condition expression to evaluate
    std::string _condition;
};

}}  // namespace

#endif // ASSERT_H_INCLUDED
