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

#ifndef LOG_H_INCLUDED
#define LOG_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Report.h>

namespace Aqua{
namespace InputOutput{

/** \class Log Log.h InputOutput/Log.h
 * Log file loader/saver. The log file is saved in an HTML file called
 * log.%d.html where %d is replaced by the first integer such that the file
 * does not exist yet.
 */
class Log : public Report
{
public:
    /** Constructor
     */
    Log();

    /** Destructor
     */
    ~Log();

    /** Save the data. The data
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** Useless method
     * @return false.
     */
    bool load(){return false;}

    /** Get the file handler
     * @return The log file handler.
     */
    FILE* fileHandler(){return _file;}
private:
    /** Create the log file
     * @return false if all gone right, true otherwise.
     */
    bool create();

    /** Close the log file
     * @return false if all gone right, true otherwise.
     */
    bool close();

    /// Output file
    FILE *_file;

};  // class InputOutput

}}  // namespaces

#endif // LOG_H_INCLUDED
