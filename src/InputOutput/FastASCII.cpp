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
 * @brief Particles plain text data files loader/saver.
 * (See Aqua::InputOutput::FastASCII for details)
 */

#include <stdlib.h>
#include <string.h>

#include <InputOutput/FastASCII.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

#ifndef MAX_LINE_LEN
    #define MAX_LINE_LEN 1024
#endif // MAX_LINE_LEN

namespace Aqua{ namespace InputOutput{

FastASCII::FastASCII(unsigned int first, unsigned int n, unsigned int iset)
    : ASCII(first, n, iset)
{
}

FastASCII::~FastASCII()
{
}

char* FastASCII::readField(const char* field,
                           const char* line,
                           unsigned int index,
                           void* data)
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    Variables* vars = C->variables();
    ArrayVariable *var = (ArrayVariable*)vars->get(field);

    // Extract the variable type data
    unsigned int n = vars->typeToN(var->type());
    size_t type_size = vars->typeToBytes(var->type());
    char *type = new char[strlen(var->type()) + 1];
    strcpy(type, var->type());
    if(strchr(type, '*'))
        strcpy(strchr(type, '*'), "");

    // Point to the chunk of data to become read
    void* ptr = (void*)((char*)data + type_size * index);

    // Start reading the sub-fields
    char* pos = (char*)line;
    for(i = 0; i < n; i++){
        char* end_pos = NULL;
        // Let's use different tools depending on the type to become read
        if(!strcmp(type, "unsigned int") ||
           strstr(type, "uivec")){
            unsigned int val = (unsigned int)strtol(pos, &end_pos, 10);
            if(pos == end_pos){
                char *msg = new char[strlen(var->type()) + 64];
                S->addMessageF(3, "Failure reading a field value\n");
                sprintf(msg, "\tWhile extracting it from \"%s\"\n", pos);
                S->addMessage(0, msg);
            }
            memcpy(ptr, &val, sizeof(unsigned int));
        }
        else if(!strcmp(type, "int") ||
                strstr(type, "ivec")){
            int val = (int)strtol(pos, &end_pos, 10);
            if(pos == end_pos){
                char *msg = new char[strlen(var->type()) + 64];
                S->addMessageF(3, "Failure reading a field value\n");
                sprintf(msg, "\tWhile extracting it from \"%s\"\n", pos);
                S->addMessage(0, msg);
            }
            memcpy(ptr, &val, sizeof(int));
        }
        else{
            float val = (float)strtof(pos, &end_pos);
            if(pos == end_pos){
                char *msg = new char[strlen(var->type()) + 64];
                S->addMessageF(3, "Failure reading a field value\n");
                sprintf(msg, "\tWhile extracting it from \"%s\"\n", pos);
                S->addMessage(0, msg);
            }
            memcpy(ptr, &val, sizeof(float));
        }

        // Go to the next field, we already asserted that there are fields
        // enough, so we don't need to care about that
        ptr = ((char*)ptr) + type_size / n;
        pos = strchr(end_pos, ',');
        if(pos)
            pos++;
    }

    delete[] type;
    return pos;
}

}}  // namespace
