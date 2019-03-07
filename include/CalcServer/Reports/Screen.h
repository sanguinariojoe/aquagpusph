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
 * @brief On screen runtime output.
 * (See Aqua::CalcServer::Reports::Screen for details)
 */

#ifndef SCREEN_H_INCLUDED
#define SCREEN_H_INCLUDED

#include <CalcServer/Reports/Report.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

/** @class Screen Screen.h CalcServer/Screen.h
 * @brief On screen runtime output.
 *
 * A runtime output is an output value that:
 *    -# Is composed by a relatively low amount of memory
 *    -# Its computation is not taking too much time
 * Therefore it could be computed and printed oftenly.
 *
 * @see Aqua::InputOutput::Logger
 */
class Screen : public Aqua::CalcServer::Reports::Report
{
public:
    /** @brief Constructor
     * @param tool_name Tool name.
     * @param fields Fields to be printed.
     * The fields are separated by commas or semicolons, and the spaces are just
     * ignored.
     * The semicolons will also force a line break in the report.
     * @param color Output color.
     * @param bold true if the text should be highlighted as bold text, false
     * otherwise.
     */
    Screen(const std::string& tool_name,
           const std::string& fields,
           const std::string color="white",
           const bool bold=false);

    /** @brief Destructor
     */
    ~Screen();

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
    /// Output color
    std::string _color;
    /// Output bold or normal flag
    bool _bold;
};

}}} // namespace

#endif // SCREEN_H_INCLUDED
