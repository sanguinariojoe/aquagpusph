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
 * @brief Fluid bounds report.
 * (See Aqua::InputOutput::Bounds for details)
 */

#ifndef REPORT_BOUNDS_H_INCLUDED
#define REPORT_BOUNDS_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Report.h>

namespace Aqua{
namespace InputOutput{

/** @class Bounds Bounds.h InputOutput/Bounds.h
 * @brief Fluid bounds report saver.
 *
 * In this report the fluid bounds box, as well as the minimum and maximum
 * velocities is saved.
 *
 * The bounds box is defined as the smallest box where all the fluid particles
 * are included inside.
 *
 * This report is a plain text file.
 * Each line corresponds to a different time instant.
 * At each line the following fields are saved, separated by tabulators:
 *   -# Time instant \f$ t \f$.
 *   -# Minimum x coordinate \f$ x_{min} \f$.
 *   -# Minimum y coordinate \f$ y_{min} \f$.
 *   -# Minimum z coordinate \f$ z_{min} \f$ (just in 3D simulations).
 *   -# Maximum x coordinate \f$ x_{max} \f$.
 *   -# Maximum y coordinate \f$ y_{max} \f$.
 *   -# Maximum z coordinate \f$ z_{max} \f$ (just in 3D simulations).
 *   -# Minimum velocity \f$ \vert \mathbf{u}_{min} \vert \f$.
 *   -# Maximum velocity \f$ \vert \mathbf{u}_{max} \vert \f$.
 *
 * The output file will be the first non existent file called `"bounds.%d.dat"`,
 * where `"%d"` is replaced by a unsigned integer.
 */
class Bounds : public Report
{
public:
    /// Constructor
    Bounds();

    /// Destructor
    ~Bounds();

    /** @brief Save the data.
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Get the bounds report file handler.
     * @return The output file handler.
     */
    FILE* fileHandler(){return _file;}
private:
    /** @brief Create the output file.
     * @return false if all gone right, true otherwise.
     */
    bool create();

    /** @brief Close the output file.
     * @return false if all gone right, true otherwise.
     */
    bool close();

    /// Output file
    FILE *_file;

};  // class InputOutput

}}  // namespaces

#endif // REPORT_BOUNDS_H_INCLUDED
