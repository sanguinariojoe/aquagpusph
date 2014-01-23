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

	/** Destructor.
	 */
	~ArgumentsManager();

    /** Parse the runtime arguments.
     * @return false if the excution must continue, true otherwise.
     */
    bool parse();

	/** Display the program usage.
	 * The program usage must be shown if the user has requested help, or if
	 * wrong/insufficient arguments have been passed.
	 */
	void displayUsage();

	/** Get the number of run arguments.
	 * @return Number of run arguments.
	 */
	int argc(){return _argc;}

	/** Get the run arguments array.
	 * @return Run arguments.
	 */
	char** argv(){return _argv;}

private:
	/// Arguments number
	int _argc;
	/// Arguments array
	char** _argv;


};  // class ArgumentsManager

}}  // namespace

#endif // ARGUMENTSMANAGER_H_INCLUDED
