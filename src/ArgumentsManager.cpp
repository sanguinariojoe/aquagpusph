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

namespace Aqua{ namespace InputOutput{

ArgumentsManager::ArgumentsManager(int argc, char **argv)
    : _argc(argc)
    , _argv(argv)
{
}

ArgumentsManager::~ArgumentsManager()
{
}

bool ArgumentsManager::parse()
{
    ScreenManager *S = ScreenManager::singleton();
    FileManager *F = FileManager::singleton();
    int index;
    char msg[256]; strcpy(msg, "");
    int opt = getopt_long( _argc, _argv, opts, longOpts, &index );
    while( opt != -1 ) {
        switch( opt ) {
            case 'i':
                F->inputFile(optarg);
                sprintf(msg,
                        "Input file = %s\n",
                        F->inputFile());
                S->addMessageF(L_INFO, msg);
                break;

            case 'v':
                printf("VERSION: ");
                printf(PACKAGE_VERSION);
                printf("\n\n");
                return true;

            case ':':
            case '?':
                S->addMessageF(L_ERROR, "Error parsing the options\n");
                printf("\n");
            case 'h':
                displayUsage();
                return true;

            default:
                S->addMessageF(L_ERROR, "Unhandled exception\n");
                displayUsage();
                return true;
        }
        opt = getopt_long( _argc, _argv, opts, longOpts, &index );
    }
    return false;
}

void ArgumentsManager::displayUsage()
{
    printf("Usage:\tAQUAgpusph [Option]...\n");
    printf("   or:\tAQUAgpusph2D [Option]...\n");
    printf("Performs particles based (SPH method) simulation (use AQUAgpusph2D for 2D simulations)\n");
    printf("\n");
    printf("Required arguments for long options are also required for the short ones.\n");
    printf("  -i, --input=INPUT            XML definition input file (Input.xml as default value)\n");
    printf("  -v, --version                Show the AQUAgpusph version\n");
    printf("  -h, --help                   Show this help page\n");
}

}}  // namespaces
