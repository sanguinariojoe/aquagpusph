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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include <ArgumentsManager.h>
#include <FileManager.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{ namespace CommandLineArgs{

// Short and long runtime options (see
// http://www.gnu.org/software/libc/manual/html_node/Getopt.html#Getopt)
static const char *opts = "i:vh";
static const struct option longOpts[] = {
    { "input", required_argument, NULL, 'i' },
    { "version", no_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
};
extern char *optarg;

/** @brief Display the program usage.
 * 
 * The program usage is shown in case the user has requested help, or
 * wrong/insufficient command line args have been provided.
 */
void displayUsage()
{
    printf("Usage:\tAQUAgpusph [Option]...\n");
    printf("   or:\tAQUAgpusph2D [Option]...\n");
    printf("Performs particles based (SPH method) simulation (use AQUAgpusph2D for 2D simulations)\n");
    printf("\n");
    printf("Required arguments for long options are also required for the short ones.\n");
    printf("  -i, --input=INPUT            XML definition input file (Input.xml by default)\n");
    printf("  -v, --version                Show the AQUAgpusph version\n");
    printf("  -h, --help                   Show this help page\n");
}

void parse(int argc, char **argv, FileManager &file_manager)
{
    int index;

    int opt = getopt_long(argc, argv, opts, longOpts, &index);
    while( opt != -1 ) {
        std::ostringstream msg;
        switch( opt ) {
            case 'i':
                file_manager.inputFile(optarg);
                msg.str(std::string());
                msg << "Input file = " << file_manager.inputFile() << std::endl;
                LOG(L_INFO, msg.str());
                break;

            case 'v':
                printf("VERSION: ");
                printf(PACKAGE_VERSION);
                printf("\n\n");
                return;

            case ':':
            case '?':
                LOG(L_ERROR, "Error parsing the runtime args\n");
                printf("\n");
            case 'h':
                displayUsage();

            default:
                LOG(L_ERROR, "Unhandled exception\n");
                displayUsage();
                throw std::invalid_argument("Invalid command line argument");
        }
        opt = getopt_long(argc, argv, opts, longOpts, &index);
    }
}

}}}  // namespaces
