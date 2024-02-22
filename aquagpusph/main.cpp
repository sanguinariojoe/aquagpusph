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
 * @brief The AQUAgpusph starting point, also known as the host.
 * (See main(int argc, char *argv[]) for details)
 */

/** @mainpage AQUAgpusph (Another QUAlity GPU-SPH).
 * Free CFD code based on SPH (Smoothed Particles Hydrodynamics), licensed
 * under GPLv3 terms, for details read provided LICENSE file.
 *
 * SPH is a meshfree Lagrangian method to solve Navier-Stokes equations).
 *
 * AQUAgpusph is based on the weakly compressible SPH branch (formerly W-SPH or
 * WCSPH), where the pressure field is linked with the density one by a
 * equation of state (EOS).
 *
 * AQUAgpusph has been accelerated with OpenCL allowing you to execute it over
 * CPU, GPU, or eventually every platform developed on future conveniently
 * adapted to the OpenCL standard.
 *
 * For the moment several platforms simultaneous usage is not supported.
 *
 * In order to learn how to use AQUAgpusph an user manual and some examples
 * are provided.
 *
 * AQUAgpusph has been developed by CEHINAV group:
 *   - Jose Luis Cercos Pita
 *   - Antonio Souto Iglesias
 *   - Leo Miguel Gonzalez
 *
 * First release of the code has been published at Wed. 2013-12-12.
 *
 * For more information about this development, please visit the following web
 * page:
 *
 * http://canal.etsin.upm.es/
 *
 * <hr>
 *
 * If you are a developer you should know that AQUAgpusph is divided in 2 main
 * pieces, the Host (see ::main) and the calculation server (see
 * Aqua::CalcServer::CalcServer).
 *
 * The first one is the manager of the simulation, being responsible of the
 * input options and files parsing (see Aqua::InputOutput::ArgumentsManager and
 * Aqua::InputOutput::FileManager), the output files writing (see
 * Aqua::InputOutput::FileManager), and the simulation running (see
 * Aqua::InputOutput::TimeManager).
 *
 * The latter one is the main worker, where the simulation computation is done
 * (see Aqua::CalcServer::CalcServer).
 *
 * You can learn more about this division in the user manual (cahpter 7) and in
 * the developers guide.
 */

#include "sphPrerequisites.hpp"
#include <Python.h>
#include "InputOutput/Logger.hpp"
#include "ArgumentsManager.hpp"
#include "FileManager.hpp"
#include "ProblemSetup.hpp"
#include "CalcServer/CalcServer.hpp"
#include "TimeManager.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

/** @namespace Aqua
 * @brief Main AQUAgpusph namespace.
 */
using namespace Aqua;

/** Aquagpusph starting point.
 *
 *  This method is associated with the host, i.e. The CPU process that controls
 *  the simulation:
 *    -# Parsing the runtime input options and the input files (see
 *       Aqua::InputOutput::ArgumentsManager and
 *       Aqua::InputOutput::FileManager).
 *    -# Setting up the simulation (see Aqua::InputOutput::ProblemSetup).
 *    -# Controlling the simulation time events (see
 *       Aqua::InputOutput::TimeManager).
 *
 * The host is also responsible to create the calculation server (see
 * Aqua::CalcServer::CalcServer).
 *
 * @param argc Number of arguments parsed by terminal.
 * @param argv Array of arguments parsed by terminal.
 */
int
main(int argc, char* argv[])
{
	std::ostringstream msg;
#ifdef HAVE_MPI
	Aqua::MPI::init(&argc, &argv);
#endif

	InputOutput::Logger* logger = new InputOutput::Logger();
	InputOutput::FileManager file_manager;

	std::cout << std::endl;
	std::cout << "\t#########################################################"
	          << std::endl;
	std::cout << "\t#                                                       #"
	          << std::endl;
	std::cout << "\t#    #    ##   #  #   #                           #     #"
	          << std::endl;
	std::cout << "\t#   # #  #  #  #  #  # #                          #     #"
	          << std::endl;
	std::cout << "\t#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #"
	          << std::endl;
	std::cout << "\t#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #"
	          << std::endl;
	std::cout << "\t#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #"
	          << std::endl;
	std::cout << "\t#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #"
	          << std::endl;
	std::cout << "\t#                            # #             #          #"
	          << std::endl;
	std::cout << "\t#                          ##  #             #          #"
	          << std::endl;
	std::cout << "\t#                                                       #"
	          << std::endl;
	std::cout << "\t#########################################################"
	          << std::endl;
	std::cout << "\tAnother QUAlity GPU-SPH, by CEHINAV (UPM) group."
	          << std::endl;
	std::cout << "\t\thttp://canal.etsin.upm.es/" << std::endl;
	std::cout << "\tAuthors:" << std::endl;
	std::cout << "\t\tJose Luis Cercos Pita" << std::endl;
	std::cout << "\t\tLeo Miguel Gonzalez" << std::endl;
	std::cout << "\t\tAntonio Souto-Iglesias" << std::endl;
	std::cout << "\t\tJavier Calderon-Sanchez" << std::endl;
	std::cout << "\tAQUAgpusph Copyright (C) 2012-2018  Jose Luis Cercos-Pita"
	          << std::endl;
	std::cout << "\tThis program comes with ABSOLUTELY NO WARRANTY; for "
	             "details see LICENSE."
	          << std::endl;
	std::cout
	    << "\tThis is free software, and you are welcome to redistribute it"
	    << std::endl;
	std::cout << "\tunder certain conditions; see LICENSE for details."
	          << std::endl;
	std::cout << std::endl;

	InputOutput::CommandLineArgs::parse(argc, argv, file_manager);

	// Now we can load the simulation definition, building the calculation
	// server
	CalcServer::CalcServer* calc_server = NULL;
	try {
		calc_server = file_manager.load();
	} catch (...) {
		if (Py_IsInitialized())
			Py_Finalize();
		return EXIT_FAILURE;
	}

	InputOutput::TimeManager t_manager(file_manager.problemSetup());

	LOG(L_INFO, "Start of simulation...\n");
	logger->printDate();
	logger->initNCurses();

	while (!t_manager.mustStop()) {
		try {
			calc_server->update(t_manager);
			file_manager.save(t_manager.time());
		} catch (const Aqua::CalcServer::user_interruption& e) {
			// The user has interrupted the simulation, just exit normally
#ifndef HAVE_MPI
			// In case of MPI it is not much useful to ask for a new file, since
			// we can hardly enforce the software to wait for the current
			// writers before MPI send a kill signal
			file_manager.save(t_manager.time());
#endif
			break;
		} catch (...) {
			logger->endNCurses();
			sleep(__ERROR_SHOW_TIME__);
			file_manager.waitForSavers();
			logger->printDate();
			msg << "Simulation finished abnormally (t = " << t_manager.time()
			    << " s)" << std::endl
			    << std::endl;
			LOG(L_INFO, msg.str());

			delete logger;
			logger = NULL;
			delete calc_server;
			calc_server = NULL;
			if (Py_IsInitialized())
				Py_Finalize();
			return EXIT_FAILURE;
		}
	}

	logger->endNCurses();
	file_manager.waitForSavers();
	logger->printDate();
	msg << "Simulation finished OK (t = " << t_manager.time() << " s)"
	    << std::endl;
	LOG(L_INFO, msg.str());

	delete logger;
	logger = NULL;
	delete calc_server;
	calc_server = NULL;
	if (Py_IsInitialized())
		Py_Finalize();
#ifdef HAVE_MPI
	Aqua::MPI::barrier(MPI_COMM_WORLD);
	Aqua::MPI::finalize();
#endif

	return EXIT_SUCCESS;
}
