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

#ifndef PARTICLES_H_INCLUDED
#define PARTICLES_H_INCLUDED

#include <sphPrerequisites.h>
#include <InputOutput/InputOutput.h>

namespace Aqua{
namespace InputOutput{

/** \class Particles Particles.h InputOutput/Particles.h
 *  Particles file loader/saver base class.
 */
class Particles : public InputOutput
{
public:
	/** Constructor
	 * @param first First particle managed by this saver/loader.
	 * @param n Number of particles managed by this saver/loader.
	 * @param ifluid Fluid index.
	 */
	Particles(unsigned int first, unsigned int n, unsigned int ifluid);

	/** Destructor
	 */
	~Particles();

    /** Save the data
     * @return false if all gone right, true otherwise.
     */
    virtual bool save(){return false;}

    /** Load the data
     * @return false if all gone right, true otherwise.
     */
    virtual bool load(){return false;}

    /** Get the last printed file.
     * @return The last printed file, NULL if a file has not been printed yet.
     */
    const char* file(){return (const char*)_output_file;}

protected:
    /** Get the bounds of the fluid managed by this class
     * @return The bounds (first and last particle indexes).
     */
    uivec2 bounds(){return _bounds;}

    /** Get the fluid index associated with this loader/saver
     * @return The fluid index.
     */
    uivec2 fluidId(){return _ifluid;}

    /** Set a new file.
     * @param filename The new file to work. Optionally a null parameter can
     * be passed in order to clear the stored file.
     */
    void file(const char* filename);

    /** Set the file name as the first ungenerate file.
     * @param basename The base name of the file. In this base name the %d
     * substring will be replaced by the first integer such that the file
     * does not exist.
     * @param startindex First index that will be tested.
     * @return false if all gone right, true otherwise.
     * @note If more than one %d substrings are found, just the first one will
     * be replaced.
     */
    bool file(const char* basename, unsigned int startindex);

private:
    /// Particles managed bounds
    uivec2 _bounds;

    /// Fluid index
    unsigned int _ifluid;

    /// Last file printed
    char* _output_file;

};  // class InputOutput

}}  // namespaces

#endif // PARTICLES_H_INCLUDED
