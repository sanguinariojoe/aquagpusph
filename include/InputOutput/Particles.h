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
 * @brief Particles files manager.
 * (See Aqua::InputOutput::Particles for details)
 */

#ifndef PARTICLES_H_INCLUDED
#define PARTICLES_H_INCLUDED

#include <sphPrerequisites.h>
#include <InputOutput/InputOutput.h>

namespace Aqua{
namespace InputOutput{

/** \class Particles Particles.h InputOutput/Particles.h
 * @brief Particles file loader/saver base class.
 *
 * In AQUAgpusph the input/output managers are divided in 3 different types:
 *   -# The simulation configuration files manager
 *   -# The report file managers
 *   -# The particles output file managers
 *
 * The particles files have 2 main objectives:
 *   -# Particles data loading at the start of simulations.
 *   -# Visualization of the simulation results.
 *
 * @see Aqua::InputOutput::InputOutput
 * @see Aqua::InputOutput::Report
 * @see Aqua::InputOutput::State
 */
class Particles : public InputOutput
{
public:
    /** @brief Constructor
     * @param first First particle managed by this saver/loader.
     * @param n Number of particles managed by this saver/loader.
     * @param ifluid Fluid index.
     */
    Particles(unsigned int first, unsigned int n, unsigned int ifluid);

    /// Destructor
    virtual ~Particles();

    /** @brief Save the data.
     * @return false if all gone right, true otherwise.
     */
    virtual bool save(){return false;}

    /** @brief Load the data.
     * @return false if all gone right, true otherwise.
     */
    virtual bool load(){return false;}

    /** @brief Get the last printed file path.
     * @return The last printed file, NULL if a file has not been printed yet.
     */
    const char* file(){return (const char*)_output_file;}

protected:
    /** @brief Get the particle index bounds of the fluid managed by this class.
     * @return The index bounds (first and last particle).
     */
    uivec2 bounds(){return _bounds;}

    /** @brief Get the fluid index associated with this class
     * @return The fluid index.
     */
    unsigned int fluidId(){return _ifluid;}

    /** @brief Set the file name.
     * @param filename The new file to save/load. Optionally a null parameter
     * can be passed in order to clear the stored file name.
     */
    void file(const char* filename);

    /** Look for the first non existing file name.
     * @param basename The base name of the file. In this base name the `%d`
     * string will be replaced by the first integer such that the file does not
     * exist in the system.
     * @param start_index First index that will be checked.
     * @param digits Number of digits of the replaced integer number. If the
     * number of digits of the integer value are greater than this value this
     * parameter will be ignored, otherwise zeroes will be appended at the left
     * of the decimal representation of the integer.
     * @return false if all gone right, true otherwise.
     * @note If more than one `"%d"` strings are found in @paramname{basename},
     * just the first one will be replaced.
     */
    bool file(const char* basename,
              unsigned int start_index,
              unsigned int digits=5);

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
