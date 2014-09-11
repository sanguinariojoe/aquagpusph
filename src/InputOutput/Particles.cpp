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
 * @brief Particles files manager.
 * (See Aqua::InputOutput::Particles for details)
 */

#include <stdlib.h>
#include <string.h>

#include <InputOutput/Particles.h>
#include <ScreenManager.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Particles::Particles(unsigned int first, unsigned int n, unsigned int iset)
    : _iset(iset)
    , _output_file(NULL)
{
    _bounds.x = first;
    _bounds.y = first + n;
}

Particles::~Particles()
{
    if(_output_file)
        delete[] _output_file;
    _output_file = NULL;
}

void Particles::file(const char* filename)
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

bool Particles::file(const char* basename,
                     unsigned int startindex,
                     unsigned int digits)
{
    FILE *f;
    char *newname = NULL, *orig_pos, *dest_pos;
    size_t len;
    unsigned int i = startindex, j;

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

        len = strlen(basename) - 1 + max(numberOfDigits(i), digits);
        newname = new char[len];

        // Copy all the string
        strcpy(newname, basename);
        // Replace the number
        dest_pos = strstr(newname, "%d");
        for(j = 0; j < digits - numberOfDigits(i); j++){
            strcpy(dest_pos, "0");
            dest_pos += 1;
        }
        sprintf(dest_pos, "%u", i);
        // Copy the rest of the original string after the inserted number
        dest_pos += numberOfDigits(i);
        orig_pos = (char*)strstr(basename, "%d") + 2;
        strcpy(dest_pos, orig_pos);

        f = fopen(newname, "r");
        if(!f){
            // We found an available slot
            break;
        }
        fclose(f);

        i++;
    }

    file(newname);
    delete[] newname;
    return false;
}

}}  // namespace
