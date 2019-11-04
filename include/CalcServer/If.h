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
 * @brief Check a condition to enable/disable all the tools in its scope
 * (See Aqua::CalcServer::If and Aqua::CalcServer::EndIf for details)
 */

#ifndef IF_H_INCLUDED
#define IF_H_INCLUDED

#include <CalcServer/Tool.h>

namespace Aqua{ namespace CalcServer{

/** @class If If.h CalcServer/If.h
 * @brief Check a condition to enable/disable all the tools in its scope
 *
 * If the result of evaluating the condition expression is equal to 0, the
 * result will be considered false, and therefore all the tools until the next
 * EndIf tool will be disabled
 */
class If : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor.
     * @param name Tool name.
     * @param condition Condition to evaluate. If the result is 0, false will be
     * considered and all the subsequent tools will be disabled until an EndIf
     * tool is reached.
     * @param once Run this tool just once. Useful to make initializations.
     */
    If(const std::string name,
       const std::string condition,
       bool once=false);

    /// Destructor.
    ~If();

    /** @brief Initialize the tool.
     */
    void setup();

    /** Open a new scope
     * @return 1
     * @note The scope shall be close at some point by an EndIf tool
     */
    const int scope_modifier(){return 1;}

    /** Return the result of the condition evaluation
     * @return true if the condition evaluation returned any value different
     * than 0, false otherwise
     */
    const bool scope_enabled(){return _result;}
protected:
    /** Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accessing the dependencies
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /// Condition expression to evaluate
    std::string _condition;
    /// Condition result
    bool _result;
};

/** @class If If.h CalcServer/If.h
 * @brief Close the scope open by a previous If call
 */
class EndIf : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor.
     * @param name Tool name.
     * @param once Run this tool just once. Useful to make initializations.
     */
    EndIf(const std::string name, bool once=false);

    /// Destructor.
    ~EndIf();

    /** @brief Initialize the tool.
     */
    void setup();

    /** Close an already open scope
     * @return -1
     */
    const int scope_modifier(){return -1;}
};

}}  // namespace

#endif // IF_H_INCLUDED
