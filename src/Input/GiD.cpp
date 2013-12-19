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
#include <math.h>
#include <string.h>

#include <Input/GiD.h>
#include <ScreenManager.h>
#include <Tokenizer/Tokenizer.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

/** Get a conveniently formatted string where the separators ",", ";",  "(",
 * ")" and tabulators are replacerd by spaces. Also the commented line will be
 * removed.
 * @param string String to reformat.
 * @return reformatted string. Empty string will be returned if it is a
 * commented, or empty line.
 */
const char* reFormatLine(const char* string){
	static char line[256]; strcpy(line, string);
	while( (line[0] == ' ') || (line[0] == '\t') )
	    strcpy(line,&line[1]);
	if( (line[0] != '#') && (line[0] != '\n') && (line[0] != EOF) ){
	    // Replace \t by spaces
	    char* replace_str = strstr(line,"\t");
	    while(replace_str){
	        strncpy(replace_str," ",1);
	        replace_str = strstr(line,"\t");
	    }
	    // Erase concatenated spaces
	    replace_str = strstr(line,"  ");
	    while(replace_str){
	        char StrBackup[256];
	        strcpy(StrBackup, &(replace_str[2]));
	        strcpy(&(replace_str[1]),StrBackup);
	        replace_str = strstr(line,"  ");
	    }
	}
	else{
	    strcpy(line,"");
	}
	return line;
}

/** Method that reads the header of a GiD exported file, and returns
 * the number of particles.
 * @note File seek will be moved to the first particle line.
 */
unsigned int readHeader(FILE *input){
	char line[256], msg[256];
	unsigned int n1, n2;
    ScreenManager *S = ScreenManager::singleton();
	while(true){
	    if(!fgets( line, 256*sizeof(char), input)){
	        S->addMessageF(3, "The file seems to be empty.\n");
	        return 0;
	    }
	    const char* line = reFormatLine(line);
	    if(!strlen(line)){      // Comment or empty line
	        continue;
	    }
	    // Two fields separated by spaces must remain
	    if(sscanf(line, "%u %u", &n1, &n2) != 2){
	        S->addMessageF(3, "Bad file header found.\n");
	        sprintf(msg, "\t\"%s\"\n", line);
	        S->addMessage(0, msg);
	        return 0;
	    }
	    // The number of particles must be positive integer
	    if(n2 <= 0){
	        S->addMessageF(3, "Incorrect number of particles readed.\n");
	        sprintf(msg, "\t\"%u\"\n", n2);
	        S->addMessage(0, msg);
	        return 0;
	    }
	    return n2;
	}
}

/** Gets all the data from a line.
 * @param line Good formatted string.
 * @param ifluid Index of the fluid.
 * @param i Index of the particle.
 * @param refd Reference density.
 * @param h Reference kernel height.
 * @param F Fluid host instance.
 * @return false if all gone right, true otherwise.
 */
bool readFields(const char* line, int ifluid, unsigned int i, float refd, float h, Fluid *F){
	unsigned int j, id, n;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	Tokenizer tok;
	const char* pos;
	bool more_data = true;                                               // false when anymore data exist
	char sentence[256]; strcpy(sentence, "");
	#ifndef HAVE_3D
	    //! Test if enoguht spaces exist.
	    pos = strstr(line," ");
	    if(!pos){
	        S->addMessageF(3, "Non readable line.\n");
	        sprintf(msg, "\t\"%s\"\n", line);
	        S->addMessage(0, msg);
	        return true;
	    }
	    for(j=0;j<4;j++){
	        pos = strstr(&pos[1]," ");
	        if(!pos){
                S->addMessageF(3, "Not enoguht data in the line.\n");
                sprintf(msg, "\t6 fields request, %u provided\n", j+1);
                S->addMessage(0, msg);
	            return true;
	        }
	    }

        // Mandatory data
        // ==============
	    id = atoi(line);
	    tok.registerVariable("id", id);
	    pos = strstr(line," ");

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].x))
	        return true;
	    tok.registerVariable("x", F->pos[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].y))
	        return true;
	    tok.registerVariable("y", F->pos[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].x))
	        return true;
	    tok.registerVariable("nx", F->normal[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].y))
	        return true;
	    tok.registerVariable("ny", F->normal[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].x))
	        return true;
	    tok.registerVariable("vx", F->v[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].y))
	        return true;
	    tok.registerVariable("vy", F->v[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    if(!pos){
	        n = strlen(sentence);
	        more_data = false;
	    }
	    else{
	        n = strlen(sentence) - strlen(pos);
	        if(!strcmp(pos, " \n")){
	            more_data = false;
	        }
	    }
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->mass[i]))
	        return true;
	    tok.registerVariable("m", F->mass[i]);

        // Optional data
        // ==============
	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->imove[i] = atoi(sentence);
	        tok.registerVariable("imove", F->imove[i]);
	    }
	    else{
	        F->imove[i] = 1;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence)-1;
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->dens[i]))
	            return true;
	        tok.registerVariable("rho", F->dens[i]);
	    }
	    else{
	        F->dens[i] = refd;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->hp[i]))
	            return true;
	        // Register variable
	        tok.registerVariable("h", F->hp[i]);
	    }
	    else{
	        F->hp[i] = h;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->ifluid[i] = atoi(sentence);
	        tok.registerVariable("ifluid", F->ifluid[i]);
	    }
	    else{
	        F->ifluid[i] = ifluid;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->drdt[i]))
	            return true;
	        tok.registerVariable("drho", F->drdt[i]);
	    }
	    else{
	        F->drdt[i] = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].x))
	            return true;
	        tok.registerVariable("fx", F->f[i].x);
	    }
	    else{
	        F->f[i].x = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].y))
	            return true;
	        tok.registerVariable("fy", F->f[i].y);
	    }
	    else{
	        F->f[i].y = 0.f;
	    }
	    more_data = false;
	#else
	    pos = strstr(line," ");
	    if(!pos){
            S->addMessageF(3, "Non readable line.\n");
            sprintf(msg, "\t\"%s\"\n", line);
            S->addMessage(0, msg);
	        return true;
	    }
	    for(j=0;j<6;j++){
	        pos = strstr(&pos[1]," ");
	        if(!pos){
                S->addMessageF(3, "Not enoguht data in the line.\n");
                sprintf(msg, "\t6 fields request, %u provided\n", j+1);
                S->addMessage(0, msg);
	            return true;
	        }
	    }

        // Mandatory data
        // ==============
	    id = atoi(line);
	    tok.registerVariable("id", id);
	    pos = strstr(line," ");

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].x))
	        return true;
	    tok.registerVariable("x", F->pos[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].y))
	        return true;
	    tok.registerVariable("y", F->pos[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].z))
	        return true;
	    tok.registerVariable("z", F->pos[i].z);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].w))
	        return true;
	    tok.registerVariable("w", F->pos[i].w);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].x))
	        return true;
	    tok.registerVariable("nx", F->normal[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].y))
	        return true;
	    tok.registerVariable("ny", F->normal[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].z))
	        return true;
	    tok.registerVariable("nz", F->normal[i].z);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].w))
	        return true;
	    tok.registerVariable("nw", F->normal[i].w);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].x))
	        return true;
	    tok.registerVariable("vx", F->v[i].x);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].y))
	        return true;
	    tok.registerVariable("vy", F->v[i].y);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].z))
	        return true;
	    tok.registerVariable("vz", F->v[i].z);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].w))
	        return true;
	    tok.registerVariable("vw", F->v[i].w);

	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    if(!pos){
	        n = strlen(sentence);
	        more_data = false;
	    }
	    else{
	        n = strlen(sentence) - strlen(pos);
	        if(!strcmp(pos, " \n")){
	            more_data = false;
	        }
	    }
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->mass[i]))
	        return true;
	    tok.registerVariable("m", F->mass[i]);

        // Optional data
        // ==============
	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->imove[i] = atoi(sentence);
	        tok.registerVariable("imove", F->imove[i]);
	    }
	    else{
	        F->imove[i] = 1;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence)-1;
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->dens[i]))
	            return true;
	        tok.registerVariable("rho", F->dens[i]);
	    }
	    else{
	        F->dens[i] = refd;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->hp[i]))
	            return true;
	        tok.registerVariable("h", F->hp[i]);
	    }
	    else{
	        F->hp[i] = h;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->ifluid[i] = atoi(sentence);
	        tok.registerVariable("ifluid", F->ifluid[i]);
	    }
	    else{
	        F->ifluid[i] = ifluid;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->drdt[i]))
	            return true;
	        tok.registerVariable("drho", F->drdt[i]);
	    }
	    else{
	        F->drdt[i] = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].x))
	            return true;
	        tok.registerVariable("fx", F->f[i].x);
	    }
	    else{
	        F->f[i].x = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].y))
	            return true;
	        tok.registerVariable("fy", F->f[i].y);
	    }
	    else{
	        F->f[i].y = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].z))
	            return true;
	        tok.registerVariable("fz", F->f[i].z);
	    }
	    else{
	        F->f[i].z = 0.f;
	    }

	    if(more_data){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            more_data = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){
	                more_data = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].w))
	            return true;
	        tok.registerVariable("fw", F->f[i].w);
	    }
	    else{
	        F->f[i].w = 0.f;
	    }
	    more_data = false;
	#endif
	return false;
}

bool loadGiD(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F)
{
	unsigned int i=0;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	FILE *input=0;
	//! 1st.- Open file
    sprintf(msg, "Loading fluid from GiD file \"%s\"\n", path);
    S->addMessageF(1, msg);
    sprintf(msg, "\tParticles from %u to %u\n", i0, i0+n-1);
    S->addMessage(0, msg);
	input = fopen(path,"r");
	if(!input){
        sprintf(msg, "Can't open the file.\n");
        S->addMessageF(3, msg);
	    return true;
	}
	//! 2nd.- Read Header.
	unsigned int num_points = readHeader(input);
	if(!num_points)
	    return true;
	if(num_points < n){
        sprintf(msg, "File doesn't contains enought particles for the fluid.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\t%u particles from %u needed.\n", num_points, n);
        S->addMessage(0, msg);
	    return true;
	}
	if(num_points > n){
        sprintf(msg, "File contains more particles than the specified for this fluid.\n");
        S->addMessageF(2, msg);
        sprintf(msg, "\t%u particles will be discarded.\n", num_points - n);
        S->addMessage(0, msg);
	}
	num_points = n;
	//! 3rd.- Read particles
	char line[256];
	unsigned int progress=-1;
	while(fgets( line, 256*sizeof(char), input)){
	    const char* line = reFormatLine(line);
	    if(!strlen(line)){
	        continue;
	    }
	    if(progress != i*100/n){
	        progress = i*100/n;
	        if(!(progress%10)){
                sprintf(msg, "\t\t%u%%\n", progress);
                S->addMessage(0, msg);
	        }
	    }
	    if(readFields(line, ifluid, i+i0, refd, h, F)){
	        return true;
	    }
	    i++;
	    if(i >= n){
	        break;
	    }
	}
	if(i<n){
        sprintf(msg, "The file ends unexpectly.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\tOnly %u particles has been readed from %u particles specified.", i, n);
        S->addMessage(0, msg);
	    return true;
	}
	fclose( input );
    sprintf(msg, "Fluid set!\n");
    S->addMessageF(1, msg);
	return false;
}

}}} // namespace
