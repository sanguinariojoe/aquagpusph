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

#include "sphPrerequisites.hpp"
#include "ProblemSetup.hpp"
#include "CalcServer/CalcServer.hpp"
#include "InputOutput/State.hpp"
#include "InputOutput/Particles.hpp"

namespace Aqua {
/// @namespace Aqua::InputOutput Input/Output data interfaces.
namespace InputOutput {

/** @class FileManager FileManager.h FileManager.h
 * @brief Input/Output files manager.
 * This class acts as a base class, controlling the subclasses which will
 * load/save the files.
 *
 * @see Aqua::InputOutput::State
 * @see Aqua::InputOutput::Particles
 * @note The Log file is controlled by Aqua::InputOutput::Logger
 */
class FileManager
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
	 * can be easily created.
	 *
	 * @param path XML input file path.
	 */
	void inputFile(std::string path);

	/** @brief Get the main XML input file path.
	 *
	 * AQUAgpusph simulations are built on top of a XML definition file. Such
	 * file can later include another XML definition files, such that modules
	 * can be easily created.
	 *
	 * @return XML input file path.
	 */
	std::string inputFile() { return _in_file; }

	/** @brief Get the simulation setup, extracted from the XML definition files
	 *
	 * AQUAgpusph simulations are built on top of a XML definition file. Such
	 * file can later include another XML definition files, such that modules
	 * can be easily created.
	 * The simulation data read from such XML files is stored in this
	 * Aqua::InputOutput::ProblemSetup structure
	 *
	 * @return Simualtion data
	 * @warning The returned Aqua::InputOutput::ProblemSetup static object is in
	 * the same scope than this class.
	 */
	ProblemSetup& problemSetup() { return _simulation; }

	/** @brief Load the input files, generating the calculation server.
	 *
	 * Depends on:
	 *    -# Aqua::InputOutput::State to load the XML definition files, storing
	 *       the data in Aqua::InputOutput::ProblemSetup.
	 *    -# Aqua::InputOutput::Particles to load the particles fields data,
	 *       storing it in Aqua::CalcServer::CalcServer.
	 *
	 * @return The built Calculation server, NULL if errors happened.
	 */
	CalcServer::CalcServer* load();

	/** @brief Save the output data files.
	 *
	 * AQUAgpusph is saving both, the XML simulation definition file, and the
	 * particles field values in the required formats. This information can be
	 * indistinctly used for postprocessing purposes, or as initial condition
	 * to resume the simulation.
	 *
	 * @param t Simulation time
	 * @warning If Python scripts are considered at the simulation, the user is
	 * responsible to save the state to can eventually resume the simulation
	 */
	void save(float t);

	/** @brief Wait for the parallel saving threads.
	 *
	 * Some savers may optionally launch parallel threads to save the data, in
	 * an asynchronous way, in order to improve the performance. In such a case,
	 * AQUAgpusph shall wait them to finish before proceeding to destroy the
	 * data
	 * @see Aqua::InputOutput::VTK::waitForSavers()
	 */
	void waitForSavers();

  private:
	/// The XML simulation definition loader/saver
	State _state;

	/// Simulation data read from XML files
	ProblemSetup _simulation;

	/// Name of the main XML input file
	std::string _in_file;

	/// The fluid loaders
	std::vector<Particles*> _loaders;

	/// The fluid savers
	std::vector<Particles*> _savers;
}; // class FileManager

}
} // namespaces

#endif // FILEMANAGER_H_INCLUDED
