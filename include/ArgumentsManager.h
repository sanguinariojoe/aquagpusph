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
 * @brief Run time input options manager.
 * (See Aqua::InputOutput::ArgumentsManager for details)
 */

#ifndef ARGUMENTSMANAGER_H_INCLUDED
#define ARGUMENTSMANAGER_H_INCLUDED

#include <sphPrerequisites.h>
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class ArgumentsManager ArgumentsManager.h ArgumentsManager.h
 * Input terminal options and parameters manager. The currently valid options
 * are:
 *   - -i, --input=INPUT            XML definition input file (Input.xml as
 *     default value).
 *   - -v, --version                Show the AQUAgpusph version.
 *   - -h, --help                   Show this help page.
 *
 * You can call `"AQUAgpusph --help"` command to see the updated list of
 * available options (in case that this documentation is outdated).
 */
class ArgumentsManager : public Aqua::Singleton<Aqua::InputOutput::ArgumentsManager>
{
public:
    /** @brief Costructor.
     * @param argc Number of run arguments
     * @param argv Array of run arguments
     */
    ArgumentsManager(int argc, char **argv);

    /** @brief Destructor.
     */
    ~ArgumentsManager();

    /** @brief Parse the runtime options.
     * @return false if the excution must continue, true otherwise.
     */
    bool parse();

    /** Display the program usage.
     *
     *  The program usage is shown weather the user has requested help, or
     *  wrong/insufficient options have been used.
     */
    void displayUsage();

    /** @brief Get the number of runtime options.
     * @return Number of runtime options.
     */
    int argc(){return _argc;}

    /** @brief Get the runtime options list.
     * @return Runtime passed options.
     */
    char** argv(){return _argv;}

private:
    /// Number of runtime options
    int _argc;
    /// List of runtime options
    char** _argv;

};  // class ArgumentsManager

}}  // namespace

#endif // ARGUMENTSMANAGER_H_INCLUDED
