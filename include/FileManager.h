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
#include <InputOutput/State.h>
#include <InputOutput/Log.h>

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

protected:

private:
	/// Name of the input file
	char* _in_file;

    /// The XML simulation definition loader
    State *_state;

    /// The output log file
    Log *_log;

};  // class FileManager

}}  // namespaces

#endif // FILEMANAGER_H_INCLUDED
