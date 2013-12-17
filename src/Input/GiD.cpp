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

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <Input/GiD.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

/** Method that takes a string, and returns another string with all fields
 * splited by spaces. , ; ( ) and tabulators & spaces will be considered
 * separators.
 * @param string String with the data.
 * @return reformatted string. Empty string will be returned if is a commented, or
 * empty line.
 */
const char* reFormatLine(const char* string){
	static char Line[256]; strcpy(Line, string);
	while( (Line[0] == ' ') || (Line[0] == '\t') )
	    strcpy(Line,&Line[1]);
	if( (Line[0] != '#') && (Line[0] != '\n') && (Line[0] != EOF) ){
	    // Replace \t by spaces
	    char* ReplacePoint = strstr(Line,"\t");
	    while(ReplacePoint){
	        strncpy(ReplacePoint," ",1);
	        ReplacePoint = strstr(Line,"\t");
	    }
	    // Erase concatenated spaces
	    ReplacePoint = strstr(Line,"  ");
	    while(ReplacePoint){
	        char StrBackup[256];
	        strcpy(StrBackup, &(ReplacePoint[2]));
	        strcpy(&(ReplacePoint[1]),StrBackup);
	        ReplacePoint = strstr(Line,"  ");
	    }
	}
	else{
	    strcpy(Line,"");
	}
	return Line;
}

/** Method that reads the header of a GiD exported file, and returns
 * the number of particles.
 * @note File seek will be moved to the first particle line.
 */
unsigned int readHeader(FILE *input){
	char Line[256];
	unsigned int n1, n2;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	while(true){
	    if(!fgets( Line, 256*sizeof(char), input)){
	        printf("ERROR (readHeader): File seems empty.\n");
	        return 0;
	    }
	    const char* line = reFormatLine(Line);
	    if(!strlen(line)){      // Comment or empty line
	        continue;
	    }
	    // Two fields separated by spaces must remain
	    if(sscanf(line, "%u %u", &n1, &n2) != 2){
	        S->addMessage(3, "(readHeader): Bad file header.\n");
	        sprintf(msg, "\t\"%s\"\n", Line);
	        S->addMessage(0, msg);
	        return 0;
	    }
	    // The number of particles must be positive integer
	    if(n2 <= 0){
	        S->addMessage(3, "(readHeader): Incorrect number of particles.\n");
	        sprintf(msg, "\t\"%u\"\n", n2);
	        S->addMessage(0, msg);
	        return 0;
	    }
	    return n2;
	}
}

/** Gets all the data from a line.
 * @param Line Good formatted string.
 * @param ifluid Index of the fluid.
 * @param i Index of the particle.
 * @param refd Reference density.
 * @param h Reference kernel height.
 * @param F Fluid host instance.
 * @remarks Good formatted line must have the first field at the start of the line,
 * without any space before it, and the rest of fields separated by spaces.
 * @return 0 if any error happens. \n
 * Error code otherwise.
 */
int readFields(const char* Line, int ifluid, unsigned int i, float refd, float h, Fluid *F){
	unsigned int j, id, n;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	Tokenizer tok;
	const char* pos;
	bool moreData = true;                                               // false when anymore data exist
	char sentence[256]; strcpy(sentence, "");
	#ifndef HAVE_3D
	    //! Test if enoguht spaces exist.
	    pos = strstr(Line," ");
	    if(!pos){
	        S->addMessage(3, "(readFields): Non readable line.\n");
	        sprintf(msg, "\t\"%s\"\n", Line);
	        S->addMessage(0, msg);
	        return 1;
	    }
	    for(j=0;j<4;j++){
	        pos = strstr(&pos[1]," ");
	        if(!pos){
                S->addMessage(3, "(readFields): Not enoguht data.\n");
                sprintf(msg, "\t6 fields request, %u provided\n", j+1);
                S->addMessage(0, msg);
	            return 2;
	        }
	    }
	    //! Read id
	    id = atoi(Line);
	    // Register variable
	    tok.registerVariable("id", id);
	    pos = strstr(Line," ");
	    //! Read x
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("x", F->pos[i].x);
	    //! Read y
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("y", F->pos[i].y);
	    //! Read nx
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("nx", F->normal[i].x);
	    //! Read ny
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("ny", F->normal[i].y);
	    //! Read vx
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("vx", F->v[i].x);
	    //! Read vy
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("vy", F->v[i].y);
	    //! Read mass
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    if(!pos){
	        n = strlen(sentence);
	        moreData = false;
	    }
	    else{
	        n = strlen(sentence) - strlen(pos);
	        if(!strcmp(pos, " \n")){           // Only a space at the end
	            moreData = false;
	        }
	    }
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->mass[i]))
	        return 3;
	    // Register variable
	    tok.registerVariable("m", F->mass[i]);
	    //! Read imove
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->imove[i] = atoi(sentence);
	        // Register variable
	        tok.registerVariable("imove", F->imove[i]);
	    }
	    else{
	        F->imove[i] = 1;
	    }
	    //! Read density
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence)-1;
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->dens[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("rho", F->dens[i]);
	    }
	    else{
	        F->dens[i] = refd;
	    }
	    //! Read kernel height
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->hp[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("h", F->hp[i]);
	    }
	    else{
	        F->hp[i] = h;
	    }
	    //! Read ifluid
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->ifluid[i] = atoi(sentence);
	        // Register variable
	        tok.registerVariable("ifluid", F->ifluid[i]);
	    }
	    else{
	        F->ifluid[i] = ifluid;
	    }
	    //! Read density rate
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->drdt[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("drho", F->drdt[i]);
	    }
	    else{
	        F->drdt[i] = 0.f;
	    }
	    //! Read fx
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].x))
	            return 4;
	        // Register variable
	        tok.registerVariable("fx", F->f[i].x);
	    }
	    else{
	        F->f[i].x = 0.f;
	    }
	    //! Read fy
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].y))
	            return 4;
	        // Register variable
	        tok.registerVariable("fy", F->f[i].y);
	    }
	    else{
	        F->f[i].y = 0.f;
	    }
	    moreData = false;
	#else
	    //! Test if enoguht spaces exist.
	    pos = strstr(Line," ");
	    if(!pos){
            S->addMessage(3, "(readFields): Non readable line.\n");
            sprintf(msg, "\t\"%s\"\n", Line);
            S->addMessage(0, msg);
	        return 1;
	    }
	    for(j=0;j<6;j++){
	        pos = strstr(&pos[1]," ");
	        if(!pos){
                S->addMessage(3, "(readFields): Not enoguht data.\n");
                sprintf(msg, "\t6 fields request, %u provided\n", j+1);
                S->addMessage(0, msg);
	            return 2;
	        }
	    }
	    //! Read id
	    id = atoi(Line);
	    // Register variable
	    tok.registerVariable("id", id);
	    pos = strstr(Line," ");
	    //! Read x
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("x", F->pos[i].x);
	    //! Read y
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("y", F->pos[i].y);
	    //! Read z
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].z))
	        return 3;
	    // Register variable
	    tok.registerVariable("z", F->pos[i].z);
	    //! Read w
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->pos[i].w))
	        return 3;
	    // Register variable
	    tok.registerVariable("w", F->pos[i].w);
	    //! Read nx
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("nx", F->normal[i].x);
	    //! Read ny
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("ny", F->normal[i].y);
	    //! Read nz
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].z))
	        return 3;
	    // Register variable
	    tok.registerVariable("nz", F->normal[i].z);
	    //! Read nw
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->normal[i].w))
	        return 3;
	    // Register variable
	    tok.registerVariable("nw", F->normal[i].w);
	    //! Read vx
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].x))
	        return 3;
	    // Register variable
	    tok.registerVariable("vx", F->v[i].x);
	    //! Read vy
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].y))
	        return 3;
	    // Register variable
	    tok.registerVariable("vy", F->v[i].y);
	    //! Read vz
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].z))
	        return 3;
	    // Register variable
	    tok.registerVariable("vz", F->v[i].z);
	    //! Read vw
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    n = strlen(sentence) - strlen(pos);
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->v[i].w))
	        return 3;
	    // Register variable
	    tok.registerVariable("vw", F->v[i].w);
	    //! Read mass
	    strcpy(sentence, pos);
	    pos = strstr(&pos[1], " ");
	    if(!pos){
	        n = strlen(sentence);
	        moreData = false;
	    }
	    else{
	        n = strlen(sentence) - strlen(pos);
	        if(!strcmp(pos, " \n")){           // Only a space at the end
	            moreData = false;
	        }
	    }
	    strcpy(&sentence[n], "");
	    if(!tok.solve(sentence, &F->mass[i]))
	        return 3;
	    // Register variable
	    tok.registerVariable("m", F->mass[i]);
	    //! Read imove
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->imove[i] = atoi(sentence);
	        // Register variable
	        tok.registerVariable("imove", F->imove[i]);
	    }
	    else{
	        F->imove[i] = 1;
	    }
	    //! Read density
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence)-1;
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->dens[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("rho", F->dens[i]);
	    }
	    else{
	        F->dens[i] = refd;
	    }
	    //! Read kernel height
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->hp[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("h", F->hp[i]);
	    }
	    else{
	        F->hp[i] = h;
	    }
	    //! Read ifluid
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        F->ifluid[i] = atoi(sentence);
	        // Register variable
	        tok.registerVariable("ifluid", F->ifluid[i]);
	    }
	    else{
	        F->ifluid[i] = ifluid;
	    }
	    //! Read density rate
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->drdt[i]))
	            return 4;
	        // Register variable
	        tok.registerVariable("drho", F->drdt[i]);
	    }
	    else{
	        F->drdt[i] = 0.f;
	    }
	    //! Read fx
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].x))
	            return 4;
	        // Register variable
	        tok.registerVariable("fx", F->f[i].x);
	    }
	    else{
	        F->f[i].x = 0.f;
	    }
	    //! Read fy
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].y))
	            return 4;
	        // Register variable
	        tok.registerVariable("fy", F->f[i].y);
	    }
	    else{
	        F->f[i].y = 0.f;
	    }
	    //! Read fz
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].z))
	            return 4;
	        // Register variable
	        tok.registerVariable("fz", F->f[i].z);
	    }
	    else{
	        F->f[i].z = 0.f;
	    }
	    //! Read fw
	    if(moreData){
	        strcpy(sentence, pos);
	        pos = strstr(&pos[1], " ");
	        if(!pos){
	            n = strlen(sentence);
	            moreData = false;
	        }
	        else{
	            n = strlen(sentence) - strlen(pos);
	            if(!strcmp(pos, " \n")){           // Only a space at the end
	                moreData = false;
	            }
	        }
	        strcpy(&sentence[n], "");
	        if(!tok.solve(sentence, &F->f[i].w))
	            return 4;
	        // Register variable
	        tok.registerVariable("fw", F->f[i].w);
	    }
	    else{
	        F->f[i].w = 0.f;
	    }
	    moreData = false;
	#endif
	return 0;
}

int loadGiD(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F)
{
	unsigned int i=0;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	FILE *input=0;
	//! 1st.- Open file
    sprintf(msg, "(loadGiD): Loading fluid from GiD file \"%s\"\n", path);
    S->addMessage(1, msg);
    sprintf(msg, "\tParticles from %u to %u\n", i0, i0+n-1);
    S->addMessage(0, msg);
	input = fopen(path,"r");
	if(!input){
        sprintf(msg, "(loadGiD): Can't open file.\n");
        S->addMessage(3, msg);
	    return 1;
	}
	//! 2nd.- Read Header.
	unsigned int nPoints = readHeader(input);
	if(!nPoints)
	    return 2;
	if(nPoints < n){
        sprintf(msg, "(loadGiD): File doesn't contains enought particles for the fluid.\n");
        S->addMessage(3, msg);
        sprintf(msg, "\t%u particles from %u needed.\n", nPoints, n);
        S->addMessage(0, msg);
	    return 3;
	}
	if(nPoints > n){
        sprintf(msg, "(loadGiD): File contains more particles of the specified for this fluid.\n");
        S->addMessage(2, msg);
        sprintf(msg, "\t%u particles will be discarded.\n", nPoints - n);
        S->addMessage(0, msg);
	}
	nPoints = n;
	//! 3rd.- Read particles
	char Line[256];
	unsigned int Percentage=-1;
	while(fgets( Line, 256*sizeof(char), input)){
	    const char* line = reFormatLine(Line);
	    if(!strlen(line)){      // Comment or empty line
	        continue;
	    }
	    if(Percentage != i*100/n){
	        Percentage = i*100/n;
	        if(!(Percentage%10)){
                sprintf(msg, "\t\t%u%%\n", Percentage);
                S->addMessage(0, msg);
	        }
	    }
	    if(readFields(line, ifluid, i+i0, refd, h, F)){
	        return 4;
	    }
	    i++;
	    if(i >= n){
	        break;
	    }
	}
	if(i<n){
        sprintf(msg, "(loadGiD): File ends unexpectly.\n");
        S->addMessage(3, msg);
        sprintf(msg, "\tOnly %u particles has been readed from %u particles specified.", i, n);
        S->addMessage(0, msg);
	    return 5;
	}
	fclose( input );
    sprintf(msg, "(loadGiD): Fluid set!...\n");
    S->addMessage(1, msg);
	return 0;
}

}}} // namespace
