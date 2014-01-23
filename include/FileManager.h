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

#ifndef FILEMANAGER_H_INCLUDED
#define FILEMANAGER_H_INCLUDED

#include <deque>

#include <sphPrerequisites.h>
#include <Singleton.h>
#include <Fluid.h>
#include <InputOutput/State.h>
#include <InputOutput/Log.h>
#include <InputOutput/Energy.h>
#include <InputOutput/Bounds.h>
#include <InputOutput/Particles.h>

namespace Aqua{
/// @namespace InputOutput Input/Output interfaces.
namespace InputOutput{

/** \class FileManager FileManager.h FileManager.h
 *  Input/Output files manager.
 */
class FileManager : public Aqua::Singleton<Aqua::InputOutput::FileManager>
{
public:
	/** Constructor
	 */
	FileManager();

	/** Destructor
	 */
	~FileManager();

	/** Set the input file.
	 * @param path Path to input file.
	 */
	void inputFile(const char* path);

	/** Get input file.
	 * @return Path to input file.
	 */
	const char* inputFile(){return (const char*)_in_file;}

	/** Get the log file handler.
	 * @return Log file handler.
	 */
	FILE* logFile(){return _log->fileHandler();}

	/** Get the energy file handler.
	 * @return Energy file handler.
	 */
	FILE* energyFile();

	/** Get the bounds file handler.
	 * @return Bounds file handler.
	 */
	FILE* boundsFile();

    /** Load the input data files.
     * @return The built Fluid manager, NULL if errors happened.
     */
    Fluid* load();

private:
	/// Name of the input file
	char* _in_file;

    /// The XML simulation definition loader
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
