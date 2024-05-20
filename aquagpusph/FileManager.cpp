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

#include "FileManager.hpp"
#include "InputOutput/Logger.hpp"
#include "InputOutput/ASCII.hpp"
#include "InputOutput/FastASCII.hpp"
#include "InputOutput/CSV.hpp"
#ifdef HAVE_VTK
#include "InputOutput/VTK.hpp"
#endif // HAVE_VTK

namespace Aqua {
namespace InputOutput {

FileManager::FileManager()
  : _state()
  , _simulation()
  , _in_file("Input.xml")
{
}

FileManager::~FileManager()
{
	for (auto loader : _loaders) {
		delete loader;
	}
	for (auto saver : _savers) {
		delete saver;
	}
}

void
FileManager::inputFile(std::string path)
{
	_in_file = path;
}

CalcServer::CalcServer*
FileManager::load()
{
	// Load the XML definition file
	_state.load(inputFile(), _simulation);

	// Setup the problem setup
	if (_simulation.sets.size() == 0) {
		LOG(L_ERROR, "No sets of particles have been added.\n");
		throw std::runtime_error("No particles sets");
	}

	// Prepare the loaders/savers
	unsigned int i = 0, offset = 0;
	for (auto set : _simulation.sets) {
		if (!set->inputFormat().compare("ASCII")) {
			ASCII* loader = new ASCII(_simulation, i, offset, set->n());
			_loaders.push_back((Particles*)loader);
		} else if (!set->inputFormat().compare("FastASCII")) {
			FastASCII* loader = new FastASCII(_simulation, i, offset, set->n());
			_loaders.push_back((Particles*)loader);
		} else if (!set->inputFormat().compare("CSV")) {
			CSV* loader = new CSV(_simulation, i, offset, set->n());
			_loaders.push_back((Particles*)loader);
		} else if (!set->inputFormat().compare("VTK")) {
#ifdef HAVE_VTK
			VTK* loader = new VTK(_simulation, i, offset, set->n());
			_loaders.push_back((Particles*)loader);
#else
			LOG(L_ERROR, "AQUAgpusph has been compiled without VTK format.\n");
			delete C;
			throw std::runtime_error("VTK support is disabled");
#endif // HAVE_VTK
		} else {
			std::ostringstream msg;
			msg << "Unknow \"" << set->inputFormat() << "\" input file format"
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Unknown input file format");
		}
		// The number of particles can be unknown yet, so better getting it from
		// the loader
		set->n(_loaders.back()->n());
		if (set->n() == 0) {
			std::ostringstream msg;
			msg << "Empty set \"" << i << "\"" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Empty set");
		}
		if (!set->outputFormat().compare("ASCII") ||
		    !set->outputFormat().compare("FastASCII")) {
			ASCII* saver = new ASCII(_simulation, i, offset, set->n());
			_savers.push_back((Particles*)saver);
		} else if (!set->outputFormat().compare("CSV")) {
			CSV* saver = new CSV(_simulation, i, offset, set->n());
			_savers.push_back((Particles*)saver);
		} else if (!set->outputFormat().compare("VTK")) {
#ifdef HAVE_VTK
			VTK* saver = new VTK(_simulation, i, offset, set->n());
			_savers.push_back((Particles*)saver);
#else
			LOG(L_ERROR, "AQUAgpusph has been compiled without VTK format.\n");
			throw std::runtime_error("VTK support is disabled");
#endif // HAVE_VTK
		} else {
			std::ostringstream msg;
			msg << "Unknow \"" << set->outputFormat() << "\" input file format"
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Unknown output file format");
		}
		offset += set->n();
		i++;
	}

	// Build the calculation server
	CalcServer::CalcServer* C = new CalcServer::CalcServer(_simulation);

	// Execute the loaders
	for (auto loader : _loaders) {
		loader->load();
	}

	return C;
}

void
FileManager::save(float t)
{
	// Execute the savers
	for (auto saver : _savers) {
		saver->save(t);
	}

	// Save the XML definition file
	_state.save(_simulation, _savers);
}

void
FileManager::waitForSavers()
{
	LOG(L_INFO, "Waiting for the writers...\n");
	for (auto saver : _savers) {
		saver->waitForSavers();
	}
}

}
} // namespace
