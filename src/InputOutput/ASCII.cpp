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
 * @brief Particles plain text data files loader/saver (with math expressions
 * evaluator).
 * (See Aqua::InputOutput::ASCII for details)
 */

#include <stdlib.h>
#include <string.h>

#include <InputOutput/ASCII.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

#ifndef MAX_LINE_LEN
    #define MAX_LINE_LEN 1024
#endif // MAX_LINE_LEN

namespace Aqua{ namespace InputOutput{

ASCII::ASCII(unsigned int first, unsigned int n, unsigned int iset)
    : Particles(first, n, iset)
{
}

ASCII::~ASCII()
{
}

bool ASCII::load()
{
    FILE *f;
    cl_int err_code;
    char msg[MAX_LINE_LEN + 64], line[MAX_LINE_LEN];
    char *pos = NULL;
    unsigned int i, j, iline, n, N, n_fields, progress;
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    loadDefault();

    sprintf(msg,
            "Loading particles from ASCII file \"%s\"\n",
            P->sets.at(setId())->inputPath());
    S->addMessageF(1, msg);

    f = fopen(P->sets.at(setId())->inputPath(), "r");
    if(!f){
        S->addMessageF(3, "The file is inaccessible.\n");
        return true;
    }

    // Assert that the number of particles is right
    n = bounds().y - bounds().x;
    N = readNParticles(f);
    if(n != N){
        sprintf(msg,
                "Expected %u particles, but the file contains %u ones.\n",
                n,
                N);
        S->addMessageF(3, msg);
        return true;
    }

    // Check the fields to read
    std::deque<char*> fields = P->sets.at(setId())->inputFields();
    if(!fields.size()){
        S->addMessageF(3, "0 fields were set to be read from the file.\n");
        return true;
    }
    bool have_r = false;
    for(i = 0; i < fields.size(); i++){
        if(!strcmp(fields.at(i), "r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        S->addMessageF(3, "\"r\" field was not set to be read from the file.\n");
        return true;
    }
    // Setup an storage
    std::deque<void*> data;
    Variables* vars = C->variables();
    n_fields = 0;
    for(i = 0; i < fields.size(); i++){
        if(!vars->get(fields.at(i))){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was not declared.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        if(!strchr(vars->get(fields.at(i))->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was declared as a scalar.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i));
        n_fields += vars->typeToN(var->type());
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure reading \"%s\" field, which has not length enough.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        void *store = malloc(typesize * n);
        if(!store){
            sprintf(msg,
                    "Failure allocating memory for \"%s\" field.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        data.push_back(store);
    }

    // Read the particles
    rewind(f);
    i = 0;
    iline = 0;
    progress = -1;
    while(fgets(line, MAX_LINE_LEN * sizeof(char), f))
    {
        iline++;

        formatLine(line);
        if(!strlen(line))
            continue;

        unsigned int n_available_fields = readNFields(line);
        if(n_available_fields != n_fields){
            sprintf(msg,
                    "Expected %u fields, but a line contains %u ones.\n",
                    n_fields,
                    n_available_fields);
            S->addMessageF(3, msg);
            sprintf(msg, "\terror found in the line %u.\n", iline);
            S->addMessage(0, msg);
            sprintf(msg, "\t\"%s\".\n", line);
            S->addMessage(0, msg);
            return true;
        }

        pos = line;
        for(j = 0; j < fields.size(); j++){
            pos = readField((const char*)fields.at(j), pos, i, data.at(j));
            if(!pos && (j != fields.size() - 1))
                return true;
        }

        i++;

        if(progress != i * 100 / n){
            progress = i * 100 / n;
            if(!(progress % 10)){
                sprintf(msg, "\t\t%u%%\n", progress);
                S->addMessage(0, msg);
            }
        }
    }

    // Send the data to the server and release it
    for(i = 0; i < fields.size(); i++){
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i));
        size_t typesize = vars->typeToBytes(var->type());
        cl_mem mem = *(cl_mem*)var->get();
        err_code = clEnqueueWriteBuffer(C->command_queue(),
                                        mem,
                                        CL_TRUE,
                                        typesize * bounds().x,
                                        typesize * n,
                                        data.at(i),
                                        0,
                                        NULL,
                                        NULL);
        free(data.at(i)); data.at(i) = NULL;
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending variable \"%s\" to the server.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
        }
    }
    data.clear();

    fclose(f);
    return false;
}

bool ASCII::save()
{
    unsigned int i, j;
    cl_int err_code;
    char msg[256];
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    TimeManager *T = TimeManager::singleton();

    std::deque<char*> fields = P->sets.at(setId())->outputFields();
    if(!fields.size()){
        S->addMessageF(3, "0 fields were set to be saved into the file.\n");
        return true;
    }

    FILE *f = create();
    if(!f)
        return true;

    // Write a head
    fprintf(f, "#########################################################\n");
    fprintf(f, "#                                                       #\n");
    fprintf(f, "#    #    ##   #  #   #                           #     #\n");
    fprintf(f, "#   # #  #  #  #  #  # #                          #     #\n");
    fprintf(f, "#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
    fprintf(f, "#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
    fprintf(f, "#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
    fprintf(f, "#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
    fprintf(f, "#                            # #             #          #\n");
    fprintf(f, "#                          ##  #             #          #\n");
    fprintf(f, "#                                                       #\n");
    fprintf(f, "#########################################################\n");
    fprintf(f, "#\n");
    fprintf(f, "#    File autogenerated by AQUAgpusph\n");
    fprintf(f, "#    t = %g s\n", T->time());
    fprintf(f, "#\n");
    fprintf(f, "#########################################################\n");
    fprintf(f, "\n");
    fflush(f);

    Variables* vars = C->variables();
    for(i = 0; i < fields.size(); i++){
        if(!vars->get(fields.at(i))){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was not declared.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        if(!strchr(vars->get(fields.at(i))->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was declared as an scalar.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i));
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    fields.at(i));
            S->addMessageF(3, msg);
            return true;
        }
    }
    std::deque<void*> data = download(fields);
    if(!data.size()){
        return true;
    }

    for(i = 0; i < bounds().y - bounds().x; i++){
        for(j = 0; j < fields.size(); j++){
            ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(j));
            const char* type_name = var->type();
            if(!strcmp(type_name, "int*")){
                int* v = (int*)data.at(j);
                fprintf(f, "%d,", v[i]);
            }
            else if(!strcmp(type_name, "unsigned int*")){
                unsigned int* v = (unsigned int*)data.at(j);
                fprintf(f, "%u,", v[i]);
            }
            else if(!strcmp(type_name, "float*")){
                float* v = (float*)data.at(j);
                fprintf(f, "%g,", v[i]);
            }
            else if(!strcmp(type_name, "ivec*")){
                #ifdef HAVE_3D
                    ivec* v = (ivec*)data.at(j);
                    fprintf(f, "%d %d %d %d,", v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    ivec* v = (ivec*)data.at(j);
                    fprintf(f, "%d %d,", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "ivec2*")){
                ivec2* v = (ivec2*)data.at(j);
                fprintf(f, "%d %d,", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "ivec3*")){
                ivec3* v = (ivec3*)data.at(j);
                fprintf(f, "%d %d %d,", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "ivec4*")){
                ivec4* v = (ivec4*)data.at(j);
                fprintf(f, "%d %d %d %d,", v[i].x, v[i].y, v[i].z, v[i].w);
            }
            else if(!strcmp(type_name, "uivec*")){
                #ifdef HAVE_3D
                    uivec* v = (uivec*)data.at(j);
                    fprintf(f, "%u %u %u %u,", v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    uivec* v = (uivec*)data.at(j);
                    fprintf(f, "%u %u,", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "uivec2*")){
                uivec2* v = (uivec2*)data.at(j);
                fprintf(f, "%u %u,", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "uivec3*")){
                uivec3* v = (uivec3*)data.at(j);
                fprintf(f, "%u %u %u,", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "uivec4*")){
                uivec4* v = (uivec4*)data.at(j);
                fprintf(f, "%u %u %u %u,", v[i].x, v[i].y, v[i].z, v[i].w);
            }
            else if(!strcmp(type_name, "vec*")){
                #ifdef HAVE_3D
                    vec* v = (vec*)data.at(j);
                    fprintf(f, "%g %g %g %g,", v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    vec* v = (vec*)data.at(j);
                    fprintf(f, "%g %g,", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "vec2*")){
                vec2* v = (vec2*)data.at(j);
                fprintf(f, "%g %g,", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "vec3*")){
                vec3* v = (vec3*)data.at(j);
                fprintf(f, "%g %g %g,", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "vec4*")){
                vec4* v = (vec4*)data.at(j);
                fprintf(f, "%g %g %g %g,", v[i].x, v[i].y, v[i].z, v[i].w);
            }
        }
        fprintf(f, "\n");
        fflush(f);
    }

    for(i = 0; i < fields.size(); i++){
        free(data.at(i)); data.at(i) = NULL;
    }
    data.clear();

    fclose(f);
    return false;
}

unsigned int ASCII::readNParticles(FILE *f)
{
    if(!f)
        return 0;

    char line[MAX_LINE_LEN];
    unsigned int n=0;

    rewind(f);
    while( fgets( line, MAX_LINE_LEN*sizeof(char), f) )
    {
        formatLine(line);
        if(!strlen(line)){
            continue;
        }

        n++;
    }

    return n;
}

void ASCII::formatLine(char* l)
{
    if(!l)
        return;
    if(!strlen(l))
        return;

    unsigned int i, len;

    // Look for a comment and discard it
    if(strchr(l, '#')){
        strcpy(strchr(l, '#'), "");
    }

    // Remove the line break if exist
    if(strchr(l, '\n')){
        strcpy(strchr(l, '\n'), "");
    }

    // Replace all the separators by commas
    const char *separators = " ;()[]{}\t";
    for(i=0; i<strlen(separators); i++){
        while(strchr(l, separators[i])){
            strncpy(strchr(l, separators[i]), ",", 1);
        }
    }

    // Remove all the concatenated separators
    while(char* mempos = strstr(l, ",,")){
        memmove(mempos, mempos + 1, strlen(mempos));
    }

    // Remove the preceding separators
    len = strlen(l);
    while(len){
        if(l[0] != ','){
            break;
        }
        memmove(l, l + 1, len);
        len--;
    }
    // And the trailing ones
    while(len){
        if(l[len - 1] != ','){
            break;
        }
        strcpy(l + len - 1, "");
        len--;
    }
}

unsigned int ASCII::readNFields(char* l)
{
    if(!l){
        return 0;
    }
    if(!strlen(l)){
        return 0;
    }

    unsigned int n = 0;
    char *pos = l;
    while(pos){
        n++;
        pos = strchr(pos, ',');
        if(pos)
            pos++;
    }

    return n;
}

char* ASCII::readField(const char* field,
                       const char* line,
                       unsigned int index,
                       void* data)
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    Variables* vars = C->variables();
    ArrayVariable *var = (ArrayVariable*)vars->get(field);

    unsigned int n = vars->typeToN(var->type());
    size_t type_size = vars->typeToBytes(var->type());

    void* ptr = (void*)((char*)data + type_size * index);
    if(vars->solve(var->type(), line, ptr)){
        return NULL;
    }

    char* pos = (char*)line;
    for(i = 0; i < n; i++){
        pos = strchr(pos, ',');
        if(pos)
            pos++;
    }
    return pos;
}


FILE* ASCII::create(){
    char *basename, msg[1024];
    size_t len;
    FILE *f;
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();

    // Create the file base name
    len = strlen(P->sets.at(setId())->outputPath()) + 8;
    basename = new char[len];
    strcpy(basename, P->sets.at(setId())->outputPath());
    strcat(basename, ".%d.dat");

    if(file(basename, 0)){
        delete[] basename;
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessageF(0, "\tHow do you received this message?.\n");
        return NULL;
    }
    delete[] basename;

    sprintf(msg, "Writing \"%s\" ASCII file...\n", file());
    S->addMessageF(1, msg);

    f = fopen(file(), "w");
    if(!f){
        sprintf(msg,
                "Failure creating the file \"%s\"\n",
                file());
        return NULL;
    }

    return f;
}

}}  // namespace
