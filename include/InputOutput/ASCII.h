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

#ifndef ASCII_H_INCLUDED
#define ASCII_H_INCLUDED

#include <stdio.h>

#include <sphPrerequisites.h>
#include <InputOutput/Particles.h>

namespace Aqua{
namespace InputOutput{

/** \class ASCII ASCII.h InputOutput/ASCII.h
 * ASCII particles data file loader/saver. These files are formatted as ASCCI
 * plain text where particles are stored by rows, where the fields are
 * separated by columns.
 * @note Comments are allowed using the symbol '#'. All the text after this
 * symbol, and in the same line, will be discarded.
 * @note The fields can be separated by the following symbols:
 *   - ' '
 *   - ','
 *   - ';'
 *   - '('
 *   - ')'
 *   - '['
 *   - ']'
 *   - '{'
 *   - '}'
 *   - '\t'
 * @note The expected fields are:
 *   -# \f$\mathbf{r}$\f.\f$x\f$
 *   -# \f$\mathbf{r}$\f.\f$y\f$
 *   -# \f$\mathbf{r}$\f.\f$z\f$ (For 3d cases only)
 *   -# \f$\mathbf{n}$\f.\f$x\f$
 *   -# \f$\mathbf{n}$\f.\f$y\f$
 *   -# \f$\mathbf{n}$\f.\f$z\f$ (For 3d cases only)
 *   -# \f$\mathbf{v}$\f.\f$x\f$
 *   -# \f$\mathbf{v}$\f.\f$y\f$
 *   -# \f$\mathbf{v}$\f.\f$z\f$ (For 3d cases only)
 *   -# \f$\frac{d \mathbf{v}}{d t}$\f.\f$x\f$
 *   -# \f$\frac{d \mathbf{v}}{d t}$\f.\f$y\f$
 *   -# \f$\frac{d \mathbf{v}}{d t}$\f.\f$z\f$ (For 3d cases only)
 *   -# \f$\rho$\f
 *   -# \f$\frac{d \rho}{d t}$\f
 *   -# \f$m$\f
 *   -# moving flag
 *   -#
 * which must be sorted inthe shown way
 */
class ASCII : public Particles
{
public:
	/** Constructor
	 * @param first First particle managed by this saver/loader.
	 * @param n Number of particles managed by this saver/loader.
	 * @param ifluid Fluid index.
	 */
	ASCII(unsigned int first, unsigned int n, unsigned int ifluid);

	/** Destructor
	 */
	~ASCII();

    /** Save the data. The data
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** Useless method
     * @return false.
     */
    bool load();

private:
    /** Count the number of Particles.
     * @param f File to be readed.
     * @return The number of particles found in the file.
     */
    unsigned int readNParticles(FILE *f);

    /** Format a line.
     * @param l String line.
     */
    void formatLine(char* l);

    /** Count the number of fields in a formated text line.
     * @param l String line.
     * @return The number of fields found in the line.
     * @warning It is assumed that formatLine has been called before.
     */
    unsigned int readNFields(char* l);

    /** Create a new file to write
     * @return The file handler, NULL if errors happened.
     */
    FILE* create();
};  // class InputOutput

}}  // namespaces

#endif // ASCII_H_INCLUDED
