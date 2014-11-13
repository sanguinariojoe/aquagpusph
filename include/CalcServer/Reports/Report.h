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
 * @brief Runtime output base class.
 * (See Aqua::CalcServer::Reports::Report for details)
 */

#ifndef REPORTS_REPORT_H_INCLUDED
#define REPORTS_REPORT_H_INCLUDED

#include <sphPrerequisites.h>

#include <deque>
#include <CalcServer/Tool.h>
#include <Variable.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{
/// @namespace Aqua::CalcServer::Reports Runtime outputs name space.
namespace Reports{

/** @class Report Report.h CalcServer/Report.h
 * @brief Runtime outputs base class.
 *
 * A runtime output is an output value that:
 *    -# Is composed by a relatively low amount of memory
 *    -# Its computation is not taking too much time
 * Therefore it could be computed and printed oftenly.
 *
 * It is tipically applied to print some relevant screen information or plot
 * friendly tabulated files.
 */
class Report : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor.
     * @param tool_name Tool name.
     * @param fields Fields to be printed.
     * The fields are separated by commas or semicolons, and the spaces are just
     * ignored.
     * The semicolons will also force a line break in the report.
     */
    Report(const char* tool_name, const char* fields);

    /** @brief Destructor
     */
    virtual ~Report();

    /** @brief Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    virtual bool setup();

    /** @brief Execute the tool.
     * @return false if all gone right, true otherwise.
     */
    virtual bool execute(){return false;}

    /** @brief Return the text string of the data to be printed.
     * @param with_title true if the report title should be inserted, false
     * otherwise.
     * @return Text string to be printed either in a file or in the screen.
     */
    const char* data(bool with_title=true);
protected:
    /** Compute the fields by lines
     * @return false if all gone right, true otherwise.
     */
    bool processFields(const char* fields);

    /** Get data string length
     * @param with_title true if the report title should be inserted, false
     * otherwise.
     * @return Length of the output string.
     */
    size_t dataLength(bool with_title=true);
private:
    /// Input fields string
    char *_fields;
    /// Output data string
    char *_data;
    /// Number of variables per line
    std::deque<unsigned int> _vars_per_line;
    /// List of variables to be printed
    std::deque<InputOutput::Variable*> _vars;
};

}}} // namespace

#endif // REPORTS_REPORT_H_INCLUDED
