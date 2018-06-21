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

#include <sphPrerequisites.h>
#include <Python.h>
#include <FileManager.h>
#include <ArgumentsManager.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <TimeManager.h>
#include <ScreenManager.h>
#include <unistd.h>

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
 *    -# Generating the simulation setup (see Aqua::InputOutput::ProblemSetup).
 *    -# Controlling the simulation time events (see
 *       Aqua::InputOutput::TimeManager), like the output files updating or the
 *       Simulation end.
 *
 * The host is responsible to create the calculation server (see
 * Aqua::CalcServer::CalcServer) and call it to work as well.
 *
 * @param argc Number of arguments parsed by terminal.
 * @param argv Array of arguments parsed by terminal.
 */
int main(int argc, char *argv[])
{
    std::ostringstream msg;
    InputOutput::FileManager file_manager;
    InputOutput::ScreenManager *S =
        new InputOutput::ScreenManager(&file_manager);
    CalcServer::CalcServer *C = NULL;
    InputOutput::TimeManager *T = NULL;

    printf("\n");
    printf("\t#########################################################\n");
    printf("\t#                                                       #\n");
    printf("\t#    #    ##   #  #   #                           #     #\n");
    printf("\t#   # #  #  #  #  #  # #                          #     #\n");
    printf("\t#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
    printf("\t#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
    printf("\t#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
    printf("\t#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
    printf("\t#                            # #             #          #\n");
    printf("\t#                          ##  #             #          #\n");
    printf("\t#                                                       #\n");
    printf("\t#########################################################\n");
    printf("\tAnother QUAlity GPU-SPH, by CEHINAV.\n");
    printf("\t\thttp://canal.etsin.upm.es/\n");
    printf("\tAuthors:\n");
    printf("\t\tJose Luis Cercos Pita\n");
    printf("\t\tLeo Miguel Gonzalez\n");
    printf("\t\tAntonio Souto-Iglesias\n");
    printf("\tAQUAgpusph Copyright (C) 2012  Jose Luis Cercos-Pita\n");
    printf("\tThis program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n");
    printf("\tThis is free software, and you are welcome to redistribute it\n");
    printf("\tunder certain conditions; see LICENSE for details.\n");
    printf("\n");

    try {
        InputOutput::CommandLineArgs::parse(argc, argv, file_manager);
    } catch (const std::invalid_argument& e) {
        return EXIT_FAILURE;
    }

    // Now we can load the simulation definition
    try {
        C = file_manager.load();
    } catch(...) {
        if(Py_IsInitialized())
            Py_Finalize();
        throw;
    }

    T = new InputOutput::TimeManager(file_manager.problemSetup());

    LOG(L_INFO, "Start of simulation...\n\n");
    S->printDate();

    while(!T->mustStop())
    {
        try {
            C->update();
            file_manager.save();
        } catch (...) {
            sleep(__ERROR_SHOW_TIME__);
            float t = T->time();
            delete S; S = NULL;
            delete T; T = NULL;
            delete C; C = NULL;
            if(Py_IsInitialized())
                Py_Finalize();
            S->printDate();
            msg << "Simulation finished abnormally (Time = " << t
                << " s)" << std::endl << std::endl;
            LOG(L_INFO, msg.str());
            throw;
        }
    }

    float t = T->time();
    delete S; S = NULL;
    delete T; T = NULL;
    delete C; C = NULL;
    if(Py_IsInitialized())
        Py_Finalize();
    S->printDate();
    msg << "Simulation finished OK (Time = " << t
        << " s)" << std::endl << std::endl;
    LOG(L_INFO, msg.str());
    return EXIT_SUCCESS;
}
