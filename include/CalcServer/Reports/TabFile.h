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
 * @brief Runtime output file.
 * (See Aqua::CalcServer::Reports::TabFile for details)
 */

#ifndef TABFILE_H_INCLUDED
#define TABFILE_H_INCLUDED

#include <fstream>
#include <CalcServer/Reports/Report.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

/** @class TabFile TabFile.h CalcServer/TabFile.h
 * @brief Runtime output file.
 *
 * A runtime output is an output value that:
 *    -# Is composed by a relatively low amount of memory
 *    -# Its computation is not taking too much time
 * Therefore it could be computed and printed oftenly.
 */
class TabFile : public Aqua::CalcServer::Reports::Report
{
public:
    /** @brief Constructor.
     * @param tool_name Tool name.
     * @param fields Fields to be printed.
     * The fields are separated by commas or semicolons, and the spaces are just
     * ignored.
     * The semicolons will also force a line break in the report.
     * @param output_file File to be written.
     * @remarks The output file will be cleared.
     */
    TabFile(const std::string tool_name,
            const std::string fields,
            const std::string output_file);

    /** @brief Destructor
     */
    ~TabFile();

    /** @brief Initialize the tool.
     */
    void setup();

protected:
    /** @brief Execute the tool.
     */
    void _execute();

private:
    /// Output file name
    std::string _output_file;
    /// Output file handler
    std::ofstream _f;
};

}}} // namespace

#endif // TABFILE_H_INCLUDED
