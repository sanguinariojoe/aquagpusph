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
 * @brief Just a small demo/template to create installable tools
 */

#ifndef INSTALLABLEDEMO_H_INCLUDED
#define INSTALLABLEDEMO_H_INCLUDED

#include <aquagpusph/CalcServer/Tool.hpp>

namespace Aqua{ namespace CalcServer{

/** @class InstallableDemo InstallableDemo.h
 * @brief Just a demo, which is printing a message in the screen.
 */
class InstallableDemo : public Aqua::CalcServer::Tool
{
public:
    /** Constructor, with same arguments than Tool.
     * @param name Tool name.
     * @param once Run this tool just once. Useful to make initializations.
     */
    InstallableDemo(const std::string tool_name, bool once);

    /** Destructor.
     */
    ~InstallableDemo();

    /** Initialize the tool.
     */
    void setup(){ Tool::setup(); };

protected:
    /** Execute the tool.
     */
    cl_event _execute(const std::vector<cl_event> events);
};

}}  // namespace

#endif // INSTALLABLEDEMO_H_INCLUDED 
