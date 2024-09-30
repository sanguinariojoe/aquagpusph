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

#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>

#include "ArgumentsManager.hpp"
#include "FileManager.hpp"
#include "InputOutput/Logger.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// Short and long runtime options (see
// http://www.gnu.org/software/libc/manual/html_node/Getopt.html#Getopt)
static const char* opts = "l:i:q:d:vh";
static const struct option longOpts[] = {
	{ "log-level", required_argument, NULL, 'l' },
	{ "input", required_argument, NULL, 'i' },
	{ "queues", required_argument, NULL, 'q' },
	{ "dimensions", required_argument, NULL, 'd' },
	{ "version", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, 'h' },
	{ NULL, no_argument, NULL, 0 }
};

namespace Aqua {
namespace InputOutput {
namespace CommandLineArgs {

/** @brief Display the program usage.
 *
 * The program usage is shown in case the user has requested help, or
 * wrong/insufficient command line args have been provided.
 */
void
displayUsage()
{
	std::cout << "Usage:\tAQUAgpusph [Option]..." << std::endl;
	std::cout << "Performs particles based (SPH method) simulations"
	          << std::endl;
	std::cout << "Required arguments for long options are also required for "
	          << "the short ones." << std::endl;
	std::cout << "  -l, --log-level=LEVEL        2 digits number indicating "
	          << "the log level of the main instance and the parallel ones. "
	          << "The levels are: 0=Debug, 1=Info, 2=Warning, 3=error, 4=none."
	          << " (02 by default)"<< std::endl;
	std::cout << "  -i, --input=INPUT            XML definition input file "
	          << "(Input.xml by default)" << std::endl;
	std::cout << "  -q, --queues=QUEUES          Maximum number of OpenCL "
	          << "command queues (15 by default)" << std::endl;
	std::cout << "  -d, --dimensions=DIMS        Simulations dimensions, "
	          << "either 2 or 3 (3 by default)" << std::endl;
	std::cout << "  -v, --version                Show the AQUAgpusph version"
	          << std::endl;
	std::cout << "  -h, --help                   Show this help page"
	          << std::endl;
}

void
parse(int argc, char** argv, FileManager& file_manager)
{
	int index, queues;
	std::string optstr;

#ifdef HAVE_MPI
	const int mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
#else
	const int mpi_rank = 0;
#endif
	int opt = getopt_long(argc, argv, opts, longOpts, &index);
	while (opt != -1) {
		std::ostringstream msg;
		switch (opt) {
			case 'l':
				optstr = trimCopy(optarg);
				if (optstr.length() != 2) {
					msg.str(std::string());
					msg << "Error parsing the log level: '"
					    << optarg << "' should be a 2 digits number"
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::invalid_argument("Invalid log level");
				}
				Aqua::InputOutput::Logger::singleton()->setLevel(
					optstr[(int)(mpi_rank > 0)] - '0');
				break;

			case 'i':
				optstr = trimCopy(optarg);
				file_manager.inputFile(optstr);
				msg.str(std::string());
				msg << "Input file = " << file_manager.inputFile() << std::endl;
				LOG(L_INFO, msg.str());
				break;

			case 'q':
				optstr = trimCopy(optarg);
				try {
					queues = std::stoi(optstr);
				} catch (const std::invalid_argument& e) {
					msg.str(std::string());
					msg << "Error parsing the number of command queues: '"
					    << optarg << "' cannot be converted to an integer"
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw;
				} catch (const std::out_of_range& e) {
					msg.str(std::string());
					msg << "Error parsing the number of command queues: '"
					    << optarg << "' is too large"
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw;
				}
				if (queues <= 0) {
					msg.str(std::string());
					msg << "Error parsing the number of command queues: '"
					    << optarg << "' is not a positive number"
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::invalid_argument("Invalid queues");
				}
				file_manager.problemSetup().nCmdQueues(queues);
				msg.str(std::string());
				msg << file_manager.problemSetup().nCmdQueues()
				    << " command queues can be spawned" << std::endl;
				LOG(L_INFO, msg.str());
				break;

			case 'd':
				optstr = trimCopy(optarg);
				if (optstr == "2") {
					file_manager.problemSetup().dims(2);
					LOG(L_INFO, "2D simulation is selected");
				} else if (optstr == "3") {
					file_manager.problemSetup().dims(3);
					LOG(L_INFO, "3D simulation is selected");
				} else {
					LOG(L_ERROR, std::string("Cannot parse the number of ") +
					             "dimension \"" + optstr + "\"\n");
					std::invalid_argument("Invalid command line argument");
				}
				break;

			case 'v':
				std::cout << "VERSION: " << PACKAGE_VERSION << std::endl
				          << std::endl;
				return;

			case 'h':
				displayUsage();
				break;

			case ':':
			case '?':
				LOG(L_ERROR, "Error parsing the runtime args\n");
				throw std::runtime_error("Invalid args");

			default:
				LOG(L_ERROR, "Invalid command line argument\n");
				displayUsage();
				throw std::invalid_argument("Invalid command line argument");
		}
		opt = getopt_long(argc, argv, opts, longOpts, &index);
	}
}

}
}
} // namespaces
