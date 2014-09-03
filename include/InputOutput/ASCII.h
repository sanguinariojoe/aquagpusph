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
 * @brief Particles plain text data files loader/saver.
 * (See Aqua::InputOutput::ASCII for details)
 */

#ifndef ASCII_H_INCLUDED
#define ASCII_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Particles.h>

namespace Aqua{
namespace InputOutput{

/** @class ASCII ASCII.h InputOutput/ASCII.h
 * @brief Plain text particles data files loader/saver.
 *
 * These files are formatted as ASCCI plain text where the particles data are
 * stored by rows, and where the fields are separated by columns.
 *
 * The fields to be saved/loaded are:
 *   -# \f$ \mathbf{r} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{r} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{r} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \mathbf{n} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{n} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{n} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \mathbf{u} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{u} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{u} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ x \f$
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ y \f$
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \rho \f$
 *   -# \f$ \frac{d \rho}{dt} \f$
 *   -# \f$ m \f$
 *   -# moving flag (see Aqua::InputOutput::Fluid::imove)
 *
 * @note Comments are allowed using the symbol `"#"`, such that all the text
 * after this symbol, and in the same line, will be discarded.
 * The fields can be separated by the following symbols:
 *   - `" "`
 *   - `","`
 *   - `";"`
 *   - `"("`
 *   - `")"`
 *   - `"["`
 *   - `"]"`
 *   - `"{"`
 *   - `"}"`
 *   - tabulator
 * @warning Saving the particles data in plain text format may be heavily hard
 * disk demanding, and therefore it is strongly recommended to consider binary
 * formats like Aqua::InputOutput::VTK.
 */
class ASCII : public Particles
{
public:
    /** @brief Constructor
     * @param first First particle managed by this saver/loader.
     * @param n Number of particles managed by this saver/loader.
     * @param ifluid Fluid index.
     */
    ASCII(unsigned int first, unsigned int n, unsigned int ifluid);

    /// Destructor
    ~ASCII();

    /** @brief Save the data.
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Load the data.
     * @return false if all gone right, true otherwise.
     */
    bool load();

private:
    /** @brief Count the number of particles present in the input file.
     * @param f File to be read.
     * @return The number of particles found in the file.
     */
    unsigned int readNParticles(FILE *f);

    /** @brief Conveniently format a read line.
     * @param l Line text.
     */
    void formatLine(char* l);

    /** @brief Count the number of fields in a text line.
     * @param l Line text.
     * @return The number of fields found in the line.
     * @warning It is assumed that the line text has been formatted calling
     * formatLine().
     */
    unsigned int readNFields(char* l);

    /** @brief Create a new file to write.
     * @return The file handler, NULL if errors happened.
     * @see Aqua::InputOutput::Particles::file(const char* basename,
     *                                         unsigned int start_index,
     *                                         unsigned int digits=5)
     */
    FILE* create();
};  // class InputOutput

}}  // namespaces

#endif // ASCII_H_INCLUDED
