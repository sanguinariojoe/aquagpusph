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

#include <InputOutput/ASCII.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <Fluid.h>
#include <AuxiliarMethods.h>
#include <Tokenizer/Tokenizer.h>

#ifndef MAX_LINE_LEN
    #define MAX_LINE_LEN 1024
#endif // MAX_LINE_LEN

#ifndef REQUESTED_FIELDS
    #ifdef HAVE_3D
        #define REQUESTED_FIELDS 17
    #else
        #define REQUESTED_FIELDS 13
    #endif
#endif // REQUESTED_FIELDS

namespace Aqua{ namespace InputOutput{

ASCII::ASCII(unsigned int first, unsigned int n, unsigned int ifluid)
    : Particles(first, n, ifluid)
{
}

ASCII::~ASCII()
{
}

bool ASCII::load()
{
    FILE *f;
    char msg[MAX_LINE_LEN + 64], sentence[MAX_LINE_LEN], *line;
    unsigned int i, iline, n, N, n_fields, progress;
	Tokenizer tok;
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	Fluid *F = Fluid::singleton();

	sprintf(msg,
            "Loading fluid from ASCII file \"%s\"\n",
            P->fluids[fluidId()].in_path);
	S->addMessageF(1, msg);

    f = fopen(P->fluids[fluidId()].in_path, "r");
    if(!f){
        S->addMessage(3, "The file is inaccesible.\n");
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
        S->addMessage(3, msg);
        return true;
    }

    // Read the particles
    rewind(f);
    i = bounds().x;
    iline = 0;
    progress = -1;
    line = new char[MAX_LINE_LEN];
	while( fgets( line, MAX_LINE_LEN*sizeof(char), f) )
	{
	    iline++;

        formatLine(line);
        if(!strlen(line)){
            delete[] line;
            line = new char[MAX_LINE_LEN];
            continue;
        }

        n_fields = readNFields(line);
        if(n_fields != REQUESTED_FIELDS){
            sprintf(msg,
                    "Expected %u fields, but a line contains %u ones.\n",
                    REQUESTED_FIELDS,
                    n_fields);
            S->addMessageF(3, msg);
            sprintf(msg, "\terror found in the line %u.\n", iline);
            S->addMessage(0, msg);
            sprintf(msg, "\t\"%s\".\n", line);
            S->addMessage(0, msg);
            return true;
        }

        tok.registerVariable("id", i);
        F->ifluid[i] = fluidId();
        tok.registerVariable("ifluid", fluidId());
        F->hp[i] = P->SPH_opts.h;
        tok.registerVariable("h", P->SPH_opts.h);

        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->pos[i].x = tok.solve(sentence);
        tok.registerVariable("x", F->pos[i].x);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->pos[i].y = tok.solve(sentence);
        tok.registerVariable("y", F->pos[i].y);

        #ifdef HAVE_3D
            line = strchr(line, ' ');
            strcpy(sentence, line);
            strcpy(strchr(sentence, ' '), "");
            F->pos[i].z = tok.solve(sentence);
            tok.registerVariable("z", F->pos[i].z);
        #endif // HAVE_3D

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->normal[i].x = tok.solve(sentence);
        tok.registerVariable("n.x", F->normal[i].x);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->normal[i].y = tok.solve(sentence);
        tok.registerVariable("n.y", F->normal[i].y);

        #ifdef HAVE_3D
            line = strchr(line, ' ');
            strcpy(sentence, line);
            strcpy(strchr(sentence, ' '), "");
            F->normal[i].z = tok.solve(sentence);
            tok.registerVariable("n.z", F->normal[i].z);
        #endif // HAVE_3D

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->v[i].x = tok.solve(sentence);
        tok.registerVariable("v.x", F->v[i].x);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->v[i].y = tok.solve(sentence);
        tok.registerVariable("v.y", F->v[i].y);

        #ifdef HAVE_3D
            line = strchr(line, ' ');
            strcpy(sentence, line);
            strcpy(strchr(sentence, ' '), "");
            F->v[i].z = tok.solve(sentence);
            tok.registerVariable("v.z", F->v[i].z);
        #endif // HAVE_3D

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->f[i].x = tok.solve(sentence);
        tok.registerVariable("dvdt.x", F->f[i].x);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->f[i].y = tok.solve(sentence);
        tok.registerVariable("dvdt.y", F->f[i].y);

        #ifdef HAVE_3D
            line = strchr(line, ' ');
            strcpy(sentence, line);
            strcpy(strchr(sentence, ' '), "");
            F->f[i].z = tok.solve(sentence);
            tok.registerVariable("dvdt.z", F->f[i].z);
        #endif // HAVE_3D

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->dens[i].x = tok.solve(sentence);
        tok.registerVariable("rho", F->dens[i]);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->drdt[i].x = tok.solve(sentence);
        tok.registerVariable("drhodt", F->drdt[i]);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->mass[i].x = tok.solve(sentence);
        tok.registerVariable("mass", F->mass[i]);

        line = strchr(line, ' ');
        strcpy(sentence, line);
        strcpy(strchr(sentence, ' '), "");
        F->imove[i].x = (int)tok.solve(sentence);
        tok.registerVariable("imove", F->imove[i]);

        i++;

        if(progress != (i - bounds().x) * 100 / n){
            progress = (i - bounds().x) * 100 / n;
            if(!(progress % 10)){
                sprintf(msg, "\t\t%u%%\n", progress);
                S->addMessage(0, msg);
            }
        }

        delete[] line;
        line = new char[MAX_LINE_LEN];
	}

    fclose(f);
    return false;
}

bool ASCII::save()
{
    unsigned int i;

    TimeManager *T = TimeManager::singleton();
	Fluid *F = Fluid::singleton();

    FILE *f = create();
    if(!f)
        return true;

    // Write a head
	fprintf(f, "\t#########################################################\n");
	fprintf(f, "\t#                                                       #\n");
	fprintf(f, "\t#    #    ##   #  #   #                           #     #\n");
	fprintf(f, "\t#   # #  #  #  #  #  # #                          #     #\n");
	fprintf(f, "\t#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
	fprintf(f, "\t#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
	fprintf(f, "\t#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
	fprintf(f, "\t#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
	fprintf(f, "\t#                            # #             #          #\n");
	fprintf(f, "\t#                          ##  #             #          #\n");
	fprintf(f, "\t#                                                       #\n");
	fprintf(f, "\t#########################################################\n");
	fprintf(f, "\t#\n");
	fprintf(f, "\t# File autogenerated by AQUAgpusph\n");
	fprintf(f, "\t# t = %g s\n", T->time());
	fprintf(f, "\n");

    for(i=bounds().x; i<bounds().y; i++){
        fprintf(f, "%g ", F->pos[i].x);
        fprintf(f, "%g ", F->pos[i].y);
        #ifdef HAVE_3D
            fprintf(f, "%g ", F->pos[i].z);
        #endif
        fprintf(f, "%g ", F->normal[i].x);
        fprintf(f, "%g ", F->normal[i].y);
        #ifdef HAVE_3D
            fprintf(f, "%g ", F->normal[i].z);
        #endif
        fprintf(f, "%g ", F->v[i].x);
        fprintf(f, "%g ", F->v[i].y);
        #ifdef HAVE_3D
            fprintf(f, "%g ", F->v[i].z);
        #endif
        fprintf(f, "%g ", F->f[i].x);
        fprintf(f, "%g ", F->f[i].y);
        #ifdef HAVE_3D
            fprintf(f, "%g ", F->f[i].z);
        #endif
        fprintf(f, "%g ", F->dens[i]);
        fprintf(f, "%g ", F->mass[i]);
        fprintf(f, "%d", F->imove[i]);

        fprintf(f, "\n");
    }

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

    unsigned int i;

    // Look for a comment and discard it
    if(strchr(l, '#')){
        strcpy(strchr(l, '#'), "");
    }

    // Remove the line break if exist
    if(strchr(line, '\n')){
        strcpy(strchr(line, '\n'), "");
    }

    // Replace all the separators by spaces
    const char *separators = ",;()[]{}\t";
    for(i=0; i<strlen(separators); i++){
        while(strchr(l, separators[i])){
            strncpy(strchr(l, separators[i]), " ", 1);
        }
    }

    // Remove all the concatenated spaces
    replace_str = strstr(line,"  ");
    while(strstr(line, "  ")){
        strcpy(strstr(line, "  "), strstr(line, "  ") + 1);
    }

    // Remove the preceeding spaces
    len = strlen(line);
    while(len){
        if(line[0] != ' '){
            break;
        }
        strcpy(line, line + 1);
        len--;
    }
    // And the trailing ones
    while(len){
        if(line[len - 1] != ' '){
            break;
        }
        strcpy(line + len - 1, "");
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

    unsigned int n = 1;
    char *pos = l;
    while(pos = strchr(pos, ' ')){
        n++;
    }

    return n;
}

FILE* ASCII::create(){
    char *basename, msg[1024];
    size_t len;
    FILE *f;
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();

    // Create the file base name
    len = strlen(P->fluids[fluidId()].out_path) + 8;
    basename = new char[len];
    strcpy(basename, P->fluids[fluidId()].out_path);
    strcat(basename, "%d.dat");

    if(file(basename, 0)){
        delete[] basename;
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessageF(0, "\tHow do you received this message?.\n");
        return NULL;
    }
    delete[] basename;

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
