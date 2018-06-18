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
 * @brief Input and output files managing.
 * (See Aqua::InputOutput::FileManager for details)
 */

#ifndef FILEMANAGER_H_INCLUDED
#define FILEMANAGER_H_INCLUDED

#include <string>
#include <vector>

#include <sphPrerequisites.h>
#include <Singleton.h>
#include <CalcServer.h>
#include <InputOutput/State.h>
#include <InputOutput/Log.h>
#include <InputOutput/Particles.h>

namespace Aqua{
/// @namespace Aqua::InputOutput Input/Output data interfaces.
namespace InputOutput{

/** @class FileManager FileManager.h FileManager.h
 * @brief Input/Output files manager.
 * This class acts as a base class, controlling the subclasses which will
 * load/save the files.
 *
 * @see Aqua::InputOutput::State
 * @see Aqua::InputOutput::Log
 * @see Aqua::InputOutput::Energy
 * @see Aqua::InputOutput::Bounds
 * @see Aqua::InputOutput::Particles
 */
class FileManager : public Aqua::Singleton<Aqua::InputOutput::FileManager>
{
public:
    /// Constructor
    FileManager();

    /// Destructor
    ~FileManager();

    /** @brief Set the main XML input file path.
     * 
     * AQUAgpusph simulations are built on top of a XML definition file. Such
     * file can later include another XML definition files, such that modules
     * can be set.
     *
     * @param path XML input file path.
     */
    void inputFile(std::string path);

    /** @brief Get the main XML input file path.
     *
     * AQUAgpusph simulations are built on top of a XML definition file. Such
     * file can later include another XML definition files, such that modules
     * can be set.
     *
     * @return XML input file path.
     */
    std::string inputFile(){return _in_file;}

    /** @brief Get the log file handler.
     *
     * AQUAgpusph is generating, during runtime, an HTML log file, placed in the
     * execution folder, and named log.X.html, where X is replaced by the first
     * unsigned integer which generates a non-existing file.
     *
     * @return Log file handler.
     * @see Aqua::InputOutput::Log
     */
    FILE* logFile();

    /** @brief Load the input files, generating the calculation server.
     * 
     * It require that:
     *    -# Aqua::InputOutput::State load the XML definition files, storing the
     *       data in Aqua::InputOutput::ProblemSetup.
     *    -# Aqua::InputOutput::Particles load the particles fields data,
     *       storing it in Aqua::CalcServer::CalcServer.
     *
     * @return The built Fluid manager, NULL if errors happened.
     */
    CalcServer::CalcServer* load();

    /** @brief Save the output data files.
     *
     * AQUAgpusph is saving both, the XML simulation definition file, and the
     * particles field values in the required formats. This information can be
     * indistinctly used for postprocessing purposes, or as initial condition
     * to resume the simulation.
     *
     * @warning If Python scripts are considered at the simulation, the user is
     * responsible to save the state to can eventually resume the simulation
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Get the last printed file for a specific particles set.
     *
     * @param set Particles set index.
     * @return The last printed file, NULL if a file has not been printed yet.
     */
    std::string file(unsigned int set);

private:
    /// Name of the main XML input file
    std::string _in_file;

    /// The XML simulation definition loader/saver
    State *_state;

    /// The output log file
    Log *_log;

    /// The fluid loaders
    std::vector<Particles*> _loaders;

    /// The fluid savers
    std::vector<Particles*> _savers;

};  // class FileManager

}}  // namespaces

#endif // FILEMANAGER_H_INCLUDED
