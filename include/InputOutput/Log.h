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
 * @brief Log file manager.
 * (See Aqua::InputOutput::Log for details)
 */

#ifndef LOG_H_INCLUDED
#define LOG_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Report.h>

namespace Aqua{
namespace InputOutput{

/** @class Log Log.h InputOutput/Log.h
 * @brief Log file saver.
 *
 * The log file is a HTML file where all the relevant events in the simulation
 * are written.
 *
 * The output file will be the first non existent file called `"log.%d.html"`,
 * where `"%d"` is replaced by a unsigned integer.
 */
class Log : public Report
{
public:
    /// Constructor
    Log();

    /// Destructor
    ~Log();

    /** @brief Save the data.
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Get the log file handler
     * @return The log file handler.
     */
    FILE* fileHandler(){return _file;}
private:
    /** @brief Create the log file
     * @return false if all gone right, true otherwise.
     */
    bool create();

    /** @brief Close the log file
     * @return false if all gone right, true otherwise.
     */
    bool close();

    /// Output file
    FILE *_file;

};  // class InputOutput

}}  // namespaces

#endif // LOG_H_INCLUDED
