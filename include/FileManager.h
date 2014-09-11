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

#include <deque>

#include <sphPrerequisites.h>
#include <Singleton.h>
#include <CalcServer.h>
#include <InputOutput/State.h>
#include <InputOutput/Log.h>
#include <InputOutput/Energy.h>
#include <InputOutput/Bounds.h>
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

    /// Set the input file path.
    /**
     * @param path Path to input file.
     */
    void inputFile(const char* path);

    /// Get input file path.
    /**
     * @return Path to input file.
     */
    const char* inputFile(){return (const char*)_in_file;}

    /// Get the log file handler.
    /**
     * @return Log file handler.
     * @see Aqua::InputOutput::Log
     */
    FILE* logFile();

    /// Get the energy file handler.
    /**
     * @return Energy file handler.
     * @see Aqua::InputOutput::Energy
     */
    FILE* energyFile();

    /// Get the bounds file handler.
    /**
     * @return Bounds file handler.
     * @see Aqua::InputOutput::Bounds
     */
    FILE* boundsFile();

    /// Load the input data files.
    /** It require that:
     *    -# Aqua::InputOutput::State should load the XML definition
     *       files storing the data in Aqua::InputOutput::ProblemSetup.
     *    -# Aqua::InputOutput::Particles should load the particles definition
     *       files storing the data in Aqua::CalcServer::CalcServer.
     *
     * @return The built Fluid manager, NULL if errors happened.
     */
    CalcServer::CalcServer* load();

    /// Save the output data files.
    /** It require that:
     *    -# Aqua::InputOutput::State should save the XML definition
     *       files taking the data from Aqua::InputOutput::ProblemSetup.
     *    -# Aqua::InputOutput::Particles should save the particles definition
     *       files stored in Aqua::InputOutput::Fluid.
     *
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /// Get the last printed file for a specific fluid.
    /**
     * @param ifluid Fluid index.
     * @return The last printed file, NULL if a file has not been printed yet.
     */
    const char* file(unsigned int ifluid);

private:
    /// Name of the main XML input file
    char* _in_file;

    /// The XML simulation definition loader/saver
    State *_state;

    /// The output log file
    Log *_log;

    /// The energy report file
    Energy *_energy;

    /// The bounds report file
    Bounds *_bounds;

    /// The fluid loaders
    std::deque<Particles*> _loaders;

    /// The fluid savers
    std::deque<Particles*> _savers;

};  // class FileManager

}}  // namespaces

#endif // FILEMANAGER_H_INCLUDED
