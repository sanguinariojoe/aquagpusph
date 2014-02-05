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

/** \mainpage AQUAgpusph (Another QUAlity GPU-SPH).
 * Free CFD code based on SPH (Smoothed Particles Hydrodynamics), licensed
 * under GPLv3 terms, for details read provided license.txt file (SPH is a
 * meshfree Lagrangian method to solve Navier-Stokes equations). \n
 * AQUAgpusph is based on the weakly compressible SPH branch (formerly W-SPH or
 * WCSPH), where the pressure field is linked with the density one by a
 * equation of state (EOS).
 * Ussually the EOS is applied: \n
 * \f$ p = \frac{c_S^2 \rho_0}{\gamma} \left( \left( \frac{\rho}{\rho_0} \right)^{\gamma} - 1 \right) , \f$ \n
 * where \f$ \gamma = 7 \f$ is usually set for water. \n
 * AQUAgpusph has been has been accelerated with OpenCL allowing you to execute
 * it over CPU or GPU based platforms, or eventually, using every platform
 * developed on future conveniently adapted to the OpenCL standard.
 * For the moment several platforms usage simultaneously is not supported. \n
 * In order to learn how to use AQUAgpusph an user manual and some examples are provided. \n
 * AQUAgpusph has been developed by CEHINAV group: \n
 * <ul><li>Cercós Pita, Jose Luis</li>
 * <li>Souto Iglesias, Antonio</li>
 * <li>Miguel Gonzalez, Leo</li></ul>
 * First release of the code have been published at Wed. 2013-12-12. \n
 * For more information about this development, please visit Canal-ETSIN web page: \n
 * http://canal.etsin.upm.es/
 */

#include <FileManager.h>
#include <ArgumentsManager.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <Fluid.h>
#include <TimeManager.h>
#include <ScreenManager.h>

using namespace Aqua;

/** Starting point
 * @param argc Number of arguments parsed by terminal.
 * @param argv Array of arguments parsed by terminal.
 */
int main(int argc, char *argv[])
{
	InputOutput::ArgumentsManager *A = NULL;
	InputOutput::FileManager *F = NULL;
	InputOutput::ProblemSetup *P = NULL;
	InputOutput::Fluid *fluid = NULL;
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
	printf("\t\tCercós Pita, Jose Luis\n");
	printf("\t\tMiguel Gonzalez, Leo\n");
	printf("\t\tSouto Iglesias, Antonio\n");
	printf("\tAQUAgpusph Copyright (C) 2012  Jose Luis Cercos Pita\n");
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
        delete fluid;
        delete C;
        delete T;
        delete S;
        return EXIT_FAILURE;
    }

    // Now we can load the data from the input files
    fluid = F->load();
    if(!fluid){
        delete A;
        delete F;
        delete P;
        delete fluid;
        delete C;
        delete T;
        delete S;
        return EXIT_FAILURE;
    }

	T = new InputOutput::TimeManager();
	C = new CalcServer::CalcServer();

	if(C->setup())
	{
        delete A;
        delete F;
        delete P;
        delete fluid;
        delete C;
        delete T;
        delete S;
		return EXIT_FAILURE;
	}

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
	        S->addMessageF(1, "Destroying fluid host layer...\n");
	        delete fluid; fluid = NULL;
	        S->addMessageF(1, "Destroying calculation server...\n");
	        delete C; C = NULL;
	        S->addMessageF(1, "Destroying problem setup...\n");
	        delete P; P = NULL;
	        S->addMessageF(1, "Destroying arguments manager...\n");
            delete A; A = NULL;
	        S->addMessageF(1, "Destroying files manager...\n");
	        delete F; F = NULL;
	        sprintf(msg, "Simulation finished abnormally (Time = %f s)\n\n", Time);
	        S->addMessageF(1, msg);
		    return EXIT_FAILURE;
		}
		fluid->retrieveData();
		if(F->save())
            return EXIT_FAILURE;
	}

	delete S; S = NULL;
    S->printDate();

	float Time = T->time();
    S->addMessageF(1, "Destroying time manager...\n");
	delete T; T = NULL;
    S->addMessageF(1, "Destroying fluid host layer...\n");
	delete fluid; fluid = NULL;
    S->addMessageF(1, "Destroying calculation server...\n");
	delete C; C = NULL;
    S->addMessageF(1, "Destroying problem setup...\n");
	delete P; P = NULL;
    S->addMessageF(1, "Destroying arguments manager...\n");
	delete A; A = NULL;
    S->addMessageF(1, "Destroying files manager...\n");
	delete F; F = NULL;
	// Exiting
    sprintf(msg, "Simulation finished OK (Time = %f s)\n\n", Time);
    S->addMessageF(1, msg);
	return EXIT_SUCCESS;
}
