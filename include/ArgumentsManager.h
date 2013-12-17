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

#ifndef ARGUMENTSMANAGER_H_INCLUDED
#define ARGUMENTSMANAGER_H_INCLUDED

#include <sphPrerequisites.h>
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** \class ArgumentsManager ArgumentsManager.h ArgumentsManager.h
 * Input terminal options and parameters manager.
 */
class ArgumentsManager : public Aqua::Singleton<Aqua::InputOutput::ArgumentsManager>
{
public:
	/** Costructor.
	 * @param argc Number of run arguments
	 * @param argv Array of run arguments
	 */
	ArgumentsManager(int argc, char **argv);

	/** Display the program usage.
	 * The program usage must be shown if the user has requested help, or if
	 * wrong/insufficient arguments have been passed.
	 */
	void displayUsage();

	/** Returns if the user has requested the output files reassembly.
	 * @return true if the output files must be reassembled, false otherwise.
	 * @note In some cases the simulations can be started from a previous
	 * saved state, in this case the output will be divided in several files,
	 * which could be convinient to reassembly for the postprocess. The
	 * reassembly process may take some time and memory.
	 * @remarks Only applies to the H5Part formatted output.
	 */
	bool mustReassembly(){return mReassembly;}

	/** Get the number of run arguments.
	 * @return Number of run arguments.
	 */
	int argc(){return mArgc;}

	/** Get the run arguments array.
	 * @return Run arguments.
	 */
	char** argv(){return mArgv;}

private:
	/// true if the output files must be reassembled, false otherwise
	bool mReassembly;

	/// Arguments number
	int mArgc;
	/// Arguments array
	char** mArgv;


};  // class ArgumentsManager

}}  // namespace

#endif // ARGUMENTSMANAGER_H_INCLUDED
