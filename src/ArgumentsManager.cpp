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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include <ArgumentsManager.h>
#include <FileManager.h>
#include <ScreenManager.h>

// short options format
static const char *opts = "i:o:nvh?";
// long options format
static const struct option longOpts[] = {
	{ "input", required_argument, NULL, 'i' },
	{ "output-refix", required_argument, NULL, 'o' },
	{ "no-reassembly", no_argument, NULL, 'n' },
	{ "version", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, 'h' },
	{ NULL, no_argument, NULL, 0 }
};
extern int optopt;
extern char *optarg;

namespace Aqua{ namespace InputOutput{

ArgumentsManager::ArgumentsManager(int argc, char **argv)
	: mReassembly(true)
	, mArgc(argc)
	, mArgv(argv)
{
    ScreenManager *S = ScreenManager::singleton();
	FileManager *F = FileManager::singleton();
	int index;
	char msg[256]; strcpy(msg, "");
	int opt = getopt_long( argc, argv, opts, longOpts, &index );
	while( opt != -1 ) {
	    switch( opt ) {
	        case 'i':
	            F->inputFile(optarg);
	            sprintf(msg, "(ArgumentsManager): Input file = %s\n", F->inputFile());
	            S->addMessage(1, msg);
	            break;

	        case 'o':
	            F->outputPrefix(optarg);
	            sprintf(msg, "(ArgumentsManager): Output prefix = %s\n", F->outputPrefix());
	            S->addMessage(1, msg);
	            break;

	        case 'n':
	            sprintf(msg, "(ArgumentsManager): Final output files will not be reassembled\n");
	            S->addMessage(1, msg);
	            mReassembly = 0;
	            break;

	        case 'v':
	            printf("VERSION: ");
	            printf(PACKAGE_VERSION);
	            printf("\n\n");
	            exit(EXIT_SUCCESS);

	        case ':':
	        case '?':
	            printf("\n");
	        case 'h':
	            displayUsage();
	            exit(EXIT_FAILURE);

	        default:
	            sprintf(msg, "(ArgumentsManager): Unhandled exception\n");
	            S->addMessage(3, msg);
	            displayUsage();
	            exit(EXIT_FAILURE);
	    }
	    opt = getopt_long( argc, argv, opts, longOpts, &index );
	}
}

void ArgumentsManager::displayUsage()
{
	printf("Usage:\tAQUAgpusph [Option]...\n");
	printf("   or:\tAQUAgpusph2D [Option]...\n");
	printf("Performs particles based (SPH method) simulation (use AQUAgpusph2D for 2D simulations)\n");
	printf("\n");
	printf("Required arguments for long options are also required for the short ones.\n");
	printf("  -i, --input=INPUT            XML definition input file (Input.xml as default value)\n");
	printf("  -o, --output-prefix=PREFIX   Prefix that will be append to output files\n");
	printf("  -n, --no-reassembly          If simulation has started from a previously saved file,\n");
	printf("                                 the final output will be divided in several files\n");
	printf("                                 which will be reassembled at the end of the\n");
	printf("                                 simulation. This process may take some time, and can\n");
	printf("                                 be disabled using this option.\n");
	printf("  -v, --version                Show the AQUAgpusph version\n");
	printf("  -h, --help                   Show this help page\n");
}

}}  // namespaces
