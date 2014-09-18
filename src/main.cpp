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

#include <Python.h>
#include <FileManager.h>
#include <ArgumentsManager.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <TimeManager.h>
#include <ScreenManager.h>

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
    InputOutput::ArgumentsManager *A = NULL;
    InputOutput::FileManager *F = NULL;
    InputOutput::ProblemSetup *P = NULL;
    CalcServer::CalcServer *C = NULL;
    InputOutput::TimeManager *T = NULL;
    InputOutput::ScreenManager *S = NULL;
    char msg[256];
    strcpy(msg, "");

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

    A = new InputOutput::ArgumentsManager(argc,argv);
    F = new InputOutput::FileManager();
    P = new InputOutput::ProblemSetup();

    // We must start reading the runtime arguments
    if(A->parse()){
        delete A;
        delete F;
        delete P;
        delete C;
        delete T;
        delete S;
        return EXIT_FAILURE;
    }

    // Now we can load the simulation definition
    C = F->load();
    if(!C){
        delete A;
        delete F;
        delete P;
        delete C;
        delete T;
        delete S;
        if(Py_IsInitialized())
            Py_Finalize();
        return EXIT_FAILURE;
    }

    T = new InputOutput::TimeManager();

    S->addMessageF(1, "Start of simulation...\n\n");
    S->printDate();
    S = new InputOutput::ScreenManager();

    while(!T->mustStop())
    {
        if(C->update()){
            float Time = T->time();
            delete S; S = NULL;
            S->printDate();
            S->addMessageF(1, "Destroying time manager...\n");
            delete T; T = NULL;
            S->addMessageF(1, "Destroying calculation server...\n");
            delete C; C = NULL;
            S->addMessageF(1, "Destroying problem setup...\n");
            delete P; P = NULL;
            S->addMessageF(1, "Destroying arguments manager...\n");
            delete A; A = NULL;
            S->addMessageF(1, "Destroying files manager...\n");
            delete F; F = NULL;
            S->addMessageF(1, "Finishing Python...\n");
            if(Py_IsInitialized())
                Py_Finalize();
            sprintf(msg, "Simulation finished abnormally (Time = %g s)\n\n", Time);
            S->addMessageF(1, msg);
            return EXIT_FAILURE;
        }

        if(F->save()){
            float Time = T->time();
            delete S; S = NULL;
            S->printDate();
            S->addMessageF(1, "Destroying time manager...\n");
            delete T; T = NULL;
            S->addMessageF(1, "Destroying calculation server...\n");
            delete C; C = NULL;
            S->addMessageF(1, "Destroying problem setup...\n");
            delete P; P = NULL;
            S->addMessageF(1, "Destroying arguments manager...\n");
            delete A; A = NULL;
            S->addMessageF(1, "Destroying files manager...\n");
            delete F; F = NULL;
            S->addMessageF(1, "Finishing Python...\n");
            if(Py_IsInitialized())
                Py_Finalize();
            sprintf(msg, "Simulation finished abnormally (Time = %g s)\n\n", Time);
            S->addMessageF(1, msg);
            return EXIT_FAILURE;
        }

        /// @todo let the tool to continue computing
        float Time = T->time();
        delete S; S = NULL;
        S->printDate();
        S->addMessageF(1, "Destroying time manager...\n");
        delete T; T = NULL;
        S->addMessageF(1, "Destroying calculation server...\n");
        delete C; C = NULL;
        S->addMessageF(1, "Destroying problem setup...\n");
        delete P; P = NULL;
        S->addMessageF(1, "Destroying arguments manager...\n");
        delete A; A = NULL;
        S->addMessageF(1, "Destroying files manager...\n");
        delete F; F = NULL;
        S->addMessageF(1, "Finishing Python...\n");
        if(Py_IsInitialized())
            Py_Finalize();
        sprintf(msg, "Simulation finished abnormally (Time = %g s)\n\n", Time);
        S->addMessageF(1, msg);
        return EXIT_FAILURE;
    }

    delete S; S = NULL;
    S->printDate();

    float Time = T->time();
    S->addMessageF(1, "Destroying time manager...\n");
    delete T; T = NULL;
    S->addMessageF(1, "Destroying calculation server...\n");
    delete C; C = NULL;
    S->addMessageF(1, "Destroying problem setup...\n");
    delete P; P = NULL;
    S->addMessageF(1, "Destroying arguments manager...\n");
    delete A; A = NULL;
    S->addMessageF(1, "Destroying files manager...\n");
    delete F; F = NULL;
    S->addMessageF(1, "Finishing Python...\n");
    if(Py_IsInitialized())
        Py_Finalize();
    sprintf(msg, "Simulation finished OK (Time = %g s)\n\n", Time);
    S->addMessageF(1, msg);
    return EXIT_SUCCESS;
}
