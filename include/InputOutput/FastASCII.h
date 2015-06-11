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
 * (See Aqua::InputOutput::FastASCII for details)
 */

#ifndef FastASCII_H_INCLUDED
#define FastASCII_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/ASCII.h>

namespace Aqua{
namespace InputOutput{

/** @class FastASCII FastASCII.h InputOutput/FastASCII.h
 * @brief Plain text particles data files loader/saver.
 *
 * These files are formatted as ASCCI plain text where the particles data are
 * stored by rows, and where the fields are separated by columns.
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
 * @remarks In the case of integer numbers (signed or unsigned) this class does
 * not care about decimal points, just truncating the value, i.e. 1.5 will be
 * interpreted as 1, and -1.5 will be interpreted as -1.
 * @warning Saving the particles data in plain text format may be heavily hard
 * disk demanding, and therefore it is strongly recommended to consider binary
 * formats like Aqua::InputOutput::VTK.
 */
class FastASCII : public ASCII
{
public:
    /** @brief Constructor
     * @param first First particle managed by this saver/loader.
     * @param n Number of particles managed by this saver/loader.
     * @param iset Particles set index.
     */
    FastASCII(unsigned int first, unsigned int n, unsigned int iset);

    /// Destructor
    ~FastASCII();

private:
    /** @brief Extract the field value from a line.
     * @param field Field name.
     * @param line Text line,
     * @param index Index of the particle to read.
     * @param data Data array.
     * @return Remaining text after extracting the field values, NULL if no
     * remaining text lefts to be read, or if the operation has failed.
     */
    char* readField(const char* field,
                    const char* line,
                    unsigned int index,
                    void* data);
};  // class InputOutput

}}  // namespaces

#endif // FastASCII_H_INCLUDED
