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

#include <vector>
#include <sphPrerequisites.h>
#include <ProblemSetup.h>
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
     * @param sim_data Simulation data
     * @param first First particle managed by this saver/loader.
     * @param n Number of particles managed by this saver/loader.
     * @param iset Particles set index.
     */
    Particles(ProblemSetup &sim_data,
              unsigned int first,
              unsigned int n,
              unsigned int iset);

    /// Destructor
    virtual ~Particles();

    /** @brief Get the last printed file path.
     * @return The last printed file, NULL if a file has not been printed yet.
     */
    const std::string file(){return _output_file;}

    /** @brief Wait for the eventual parallel saving threads.
     *
     * Some savers may optionally launch parallel threads to save the data, in
     * an asynchronous way, in order to improve the performance. In such a case,
     * AQUAgpusph shall wait them to finish before proceeding to destroy the
     * data
     */
    virtual void waitForSavers() {return;}
protected:
    /** @brief Get the simulation data structure
     *
     * @return Simulation data
     */
    ProblemSetup& simData() {return _sim_data;}

    /** @brief Get the particle index bounds of the "set of particles" managed
     * by this class.
     * @return The index bounds (first and last particle).
     */
    uivec2 bounds(){return _bounds;}

    /** @brief Get the "particles set" index associated with this class
     * @return The "particles index" index.
     */
    unsigned int setId(){return _iset;}

    /** @brief Register some default arrays:
     *   -# iset
     *   -# id_sorted
     *   -# id_unsorted
     */
    void loadDefault();

    /** @brief Set the file name.
     * @param filename The new file to save/load. Optionally a null parameter
     * can be passed in order to clear the stored file name.
     */
    void file(const std::string filename){_output_file = filename;};

    /** Look for the first non-existing file path
     * @param basename The base name of the file
     * @param start_index First index that will be checked.
     * @param digits Number of digits of the replaced integer number. If the
     * number of digits of the integer value are greater than this value this
     * parameter will be ignored, otherwise zeroes will be appended at the left
     * of the decimal representation of the integer.
     * @return The next non-existing file index.
     * @see Aqua::newFilePath()
     */
    unsigned int file(const std::string basename,
                      unsigned int start_index,
                      unsigned int digits=5);

    /** Download the data from the device, and store it
     * @param fields Fields to download
     * @return host allocated memory
     * @note The returned data must be manually cleared.
     */
    std::vector<void*> download(std::vector<std::string> fields);
private:
    /** Remove the content of the data list.
     * @param data List of memory allocated arrays to be cleared.
     */
    void clearList(std::vector<void*> *data);

    /// Simulation data
    ProblemSetup _sim_data;

    /// Particles managed bounds
    uivec2 _bounds;

    /// Fluid index
    unsigned int _iset;

    /// Last file printed
    std::string _output_file;
};  // class InputOutput

}}  // namespaces

#endif // PARTICLES_H_INCLUDED
