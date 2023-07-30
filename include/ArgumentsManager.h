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
#include <FileManager.h>

namespace Aqua {
namespace InputOutput {
namespace CommandLineArgs {

/** @brief Parse the command line arguments.
 *
 * @return false if the excution must continue, true otherwise.
 */
void
parse(int argc, char** argv, FileManager& file_manager);

}
}
} // namespace

#endif // ARGUMENTSMANAGER_H_INCLUDED
