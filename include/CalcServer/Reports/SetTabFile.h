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
 * @brief Runtime particles set output file.
 * (See Aqua::CalcServer::Reports::SetSetTabFile for details)
 */

#ifndef SETSetTabFile_H_INCLUDED
#define SETSetTabFile_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sphPrerequisites.h>
#include <CalcServer/Reports/Report.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

/** @class SetTabFile SetTabFile.h CalcServer/SetTabFile.h
 * @brief Particles set runtime output.
 *
 * A runtime output is an output value that:
 *    -# Is composed by a relatively low amount of memory
 *    -# Its computation is not taking too much time
 * Therefore it could be computed and printed oftenly.
 *
 * This tool is printing the selected properties of a set of particles in a
 * tabulated file with the following columns:
 *    - Time instant
 *    - First particle, first property
 *    - First particle, second property
 *    - ...
 *    - First particle, last property
 *    - Second particle, first property
 *    - ...
 *    - Second particle, last property
 *    - ...
 *    - Last particle, first property
 *    - ...
 *    - Last particle, last property
 *
 * And therefore \f$ n_{prop} \cdot n_{parts} \f$ fields should be doownloaded
 * and printed in plain text, so be careful about what particles sets and fields
 * are requested.
 */
class SetTabFile : public Aqua::CalcServer::Reports::Report
{
public:
    /** @brief Constructor.
     * @param tool_name Tool name.
     * @param fields Fields to be printed.
     * The fields are separated by commas or semicolons, and the spaces are just
     * ignored.
     * The semicolons will also force a line break in the report.
     * @param first First particle managed by this report (unsorted indexes).
     * @param n Number of particles managed by this report (unsorted indexes).
     * @param output_file File to be written.
     * @param ipf Iterations per frame, 0 to just ignore this printing criteria.
     * @param fps Frames per second, 0 to just ignore this printing criteria.
     * @remarks The output file will be cleared.
     */
    SetTabFile(const std::string tool_name,
               const std::string fields,
               unsigned int first,
               unsigned int n,
               const std::string output_file,
               unsigned int ipf=1,
               float fps=0.f);

    /** @brief Destructor
     */
    ~SetTabFile();

    /** @brief Initialize the tool.
     */
    void setup();

protected:
    /** @brief Execute the tool.
     */
    void _execute();

    /** @brief Get the particle index bounds of the "set of particles" managed
     * by this class.
     * @return The index bounds (first and last particle).
     */
    uivec2 bounds(){return _bounds;}

    /** Download the data from the device, and store it.
     * @param vars Fields to download.
     * @return host allocated memory. A clear list if errors happened.
     * @note The returned data must be manually cleared.
     */
    std::vector<void*> download(std::vector<InputOutput::Variable*> vars);

private:
    /** Remove the content of the data list.
     * @param data List of memory allocated arrays to be cleared.
     */
    void clearList(std::vector<void*> *data);

    /// Particles managed bounds
    uivec2 _bounds;

    /// Output file name
    std::string _output_file;
    /// Output file handler
    std::ofstream _f;
};

}}} // namespace

#endif // SETSetTabFile_H_INCLUDED
