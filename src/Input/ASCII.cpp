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
#include <Input/ASCII.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

bool loadASCII(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F)
{
	unsigned int i=0;
	FILE *input=0;
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	//! 1st.- Open file
	sprintf(msg, "(loadASCII): Loading fluid from ASCII file \"%s\"\n", path);
	S->addMessage(1, msg);
	sprintf(msg, "\tParticles from %u to %u\n", i0, i0+n-1);
	S->addMessage(0, msg);
	input = fopen(path,"r");
	if(!input){
	    S->addMessage(3, "(loadASCII): Can't open file.\n");
	    return true;
	}
	//! 1st.- Get the number of data lines.
	char Line[256];
	unsigned int nPoints=0;
	int LineStartChar;
	while( fgets( Line, 256*sizeof(char), input) )
	{
	    LineStartChar=0;
	    while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
	        LineStartChar++;
	    if( (Line[LineStartChar] != '#') && (Line[LineStartChar] != '\n') && (Line[LineStartChar] != EOF) )
	        nPoints++;
	}
	if(nPoints < n-1){
	    S->addMessage(3, "(loadASCII): Seems that file doesn't contain enought particles.\n");
	    sprintf(msg, "\t%u particles, %u particles needed.\n", nPoints, n);
        S->addMessage(0, msg);
	    return true;
	}
	else if(nPoints > n){
	    S->addMessage(2, "(loadASCII): File contains more particles than the fluid.\n");
	    sprintf(msg, "\t%u particles will be discarded.\n", nPoints - n);
        S->addMessage(0, msg);
	}
	//! 2nd.- Reads data
	rewind(input);
	unsigned int Percentage=-1;
	while( fgets( Line, 256*sizeof(char), input) )
	{
	    LineStartChar=0;
	    while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
	        LineStartChar++;
	    if( (Line[LineStartChar] != '#') && (Line[LineStartChar] != '\n') && (Line[LineStartChar] != EOF) ){
	        if(i >= n){
	            break;
	        }
	        if(Percentage != i*100/n){
	            Percentage = i*100/n;
	            if(!(Percentage%10)){
                    sprintf(msg, "\t\t%u%%\n", Percentage);
                    S->addMessage(0, msg);
	            }
	        }
	        // Replace , ; ( ) - \t by spaces
	        char *ReplacePoint=0;
	        ReplacePoint = strstr(Line,",");
	        while(ReplacePoint) {
	            strncpy(ReplacePoint," ",1);
	            ReplacePoint = strstr(Line,",");
	        }
	        ReplacePoint = strstr(Line,";");
	        while(ReplacePoint){
	            strncpy(ReplacePoint," ",1);
	            ReplacePoint = strstr(Line,";");
	        }
	        ReplacePoint = strstr(Line,"(");
	        while(ReplacePoint){
	            strncpy(ReplacePoint," ",1);
	            ReplacePoint = strstr(Line,"(");
	        }
	        ReplacePoint = strstr(Line,")");
	        while(ReplacePoint){
	            strncpy(ReplacePoint," ",1);
	            ReplacePoint = strstr(Line,")");
	        }
	        ReplacePoint = strstr(Line,"\t");
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
	        #ifndef HAVE_3D
	            int readed = sscanf(Line, "%g %g %g %g %g %g %g %i %g %g %i %g %g %g",
	               &F->pos[i+i0].x,&F->pos[i+i0].y,&F->normal[i+i0].x,&F->normal[i+i0].y,
	               &F->v[i+i0].x,&F->v[i+i0].y,&F->mass[i+i0],
	               &F->imove[i+i0],&F->dens[i+i0],&F->hp[i+i0],
	               &F->ifluid[i+i0],&F->drdt[i+i0],&F->f[i+i0].x,&F->f[i+i0].y);
	            if( readed < 7)
	            {
	                S->addMessage(3, "(loadASCII): Minimum data has not been provided.\n");
                    sprintf(msg, "\t\"%s\"\n",Line);
                    S->addMessage(0, msg);
                    sprintf(msg, "\tAt least Pos, vel and mass must be set.\n");
                    S->addMessage(0, msg);
	                return true;
	            }
	            // Complete unset data
	            switch(readed){
	                case 7:
	                    F->imove[i+i0]   = 1;
	                case 8:
	                    F->dens[i+i0]    = refd;
	                case 9:
	                    F->hp[i+i0]      = h;
	                case 10:
	                    F->ifluid[i+i0]  = ifluid;
	                case 11:
	                    F->drdt[i+i0]    = 0.f;
	                case 12:
	                    F->f[i+i0].x     = 0.f;
	                case 13:
	                    F->f[i+i0].y     = 0.f;
	                default:
	                    break;
	            }
	        #else
	            int readed = sscanf(Line, "%g %g %g %g %g %g %g %g %g %g %i %g %g %i %g %g %g %g",
	               &F->pos[i+i0].x,&F->pos[i+i0].y,&F->pos[i+i0].z,
	               &F->normal[i+i0].x,&F->normal[i+i0].y,&F->normal[i+i0].z,
	               &F->v[i+i0].x,&F->v[i+i0].y,&F->v[i+i0].z,&F->mass[i+i0],
	               &F->imove[i+i0],&F->dens[i+i0],&F->hp[i+i0],
	               &F->ifluid[i+i0],&F->drdt[i+i0],&F->f[i+i0].x,&F->f[i+i0].y,&F->f[i+i0].z);
	            if( readed < 10)
	            {
	                S->addMessage(3, "(loadASCII): Minimum data has not been provided.\n");
                    sprintf(msg, "\t\"%s\"\n",Line);
                    S->addMessage(0, msg);
                    sprintf(msg, "\tAt least Pos, vel and mass must be set.\n");
                    S->addMessage(0, msg);
	                return true;
	            }
	            // Complete unset data
	            switch(readed){
	                case 10:
	                    F->imove[i+i0]   = 1;
	                case 11:
	                    F->dens[i+i0]    = refd;
	                case 12:
	                    F->hp[i+i0]      = h;
	                case 13:
	                    F->ifluid[i+i0]  = ifluid;
	                case 14:
	                    F->drdt[i+i0]    = 0.f;
	                case 15:
	                    F->f[i+i0].x     = 0.f;
	                case 16:
	                    F->f[i+i0].y     = 0.f;
	                case 17:
	                    F->f[i+i0].z     = 0.f;
	                default:
	                    break;
	            }
	            F->pos[i+i0].w    = 0.f;
	            F->normal[i+i0].w = 0.f;
	            F->v[i+i0].w      = 0.f;
	            F->f[i+i0].w      = 0.f;
	        #endif
	        i++;
	    }
	}
	fclose( input );
	S->addMessage(1, "(loadASCII): Fluid set!...\n");
	return false;
}

}}} // namespaces
