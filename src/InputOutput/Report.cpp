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

#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <InputOutput/Report.h>
#include <ScreenManager.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Report::Report()
    : _output_file(NULL)
{
}

Report::~Report()
{
    if(_output_file)
        delete[] _output_file;
    _output_file = NULL;
}

void Report::file(const char* filename)
{
    size_t len;

    if(_output_file)
        delete[] _output_file;
    _output_file = NULL;

    if(!filename)
        return;

    len = strlen(filename) + 1;
    _output_file = new char[len];
    strcpy(_output_file, filename);
}

bool Report::file(const char* basename, unsigned int startindex)
{
    FILE *f;
    char *newname = NULL, *orig_pos, *dest_pos;
    size_t len;
    unsigned int i = startindex;

    if(!basename)

    if(!strstr(basename, "%d")){
        // We cannot replace nothing in the basename, just test if the file
        // does not exist
        f = fopen(basename, "r");
        if(f){
            // The fail already exist, so we cannot operate
            fclose(f);
            return true;
        }

        file(basename);
        return false;
    }

    while(true){
        if(newname)
            delete[] newname;
        newname = NULL;

        len = strlen(basename) - 1 + numberOfDigits(i);
        newname = new char[len];

        // Copy all the string
        strcpy(newname, basename);
        // Replace the number
        dest_pos = strstr(newname, "%d");
        sprintf(dest_pos, "%u", i);
        // Copy the rest of the original string after the inserted number
        dest_pos += numberOfDigits(i);
        orig_pos = strstr(basename, "%d") + 2;
        strcpy(dest_pos, orig_pos);

        f = fopen(newname, "r");
        if(!f){
            // We found an available slot
            break;
        }
        i++;
    }

    file(newname);
    delete[] newname;
    return false;
}

}}  // namespace
