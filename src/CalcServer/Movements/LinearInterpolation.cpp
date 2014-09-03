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
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>
#include <unistd.h>
#include <errno.h>

#include <ScreenManager.h>
#include <TimeManager.h>
#include <CalcServer/Movements/LinearInterpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

LinearInterpolation::LinearInterpolation(const char *data_file)
    : _data_file(0)
    , _time(0.f)
{
    _data.clear();
    _prev_data.clear();
    _next_data.clear();

    if(data_file)
        open(data_file);
}

LinearInterpolation::~LinearInterpolation()
{
    if(_data_file) fclose(_data_file); _data_file=0;
    _data.clear();
    _prev_data.clear();
    _next_data.clear();
}

std::deque<float> LinearInterpolation::update(float t)
{
    unsigned int i;
    if(!_data_file){
        _data.clear();
        return _data;
    }
    // Rewind file if the time is lower than the last stored one
    if(t < _time){
        rewind(_data_file);
        _prev_data.clear();
        _next_data.clear();
    }
    // Ensure that the new time is not bounded by the previously readed ones
    _time = t;
    if(_next_data.size() > 0){
        if(_next_data.at(0) > _time){
            float factor = (_time - _prev_data.at(0)) / (_next_data.at(0) - _prev_data.at(0));
            _data.clear();
            for(i=0;i<_next_data.size();i++){
                float dData = _next_data.at(i) - _prev_data.at(i);
                _data.push_back( _prev_data.at(i) + factor*dData );
            }
            return _data;
        }
    }
    // Read lines while a valid point is located
    _prev_data = _next_data;
    _next_data = readLine();
    if(!_next_data.size()){      // Last point reached
        _data = _prev_data;
        return _data;
    }
    if(_next_data.size() != _prev_data.size()){   // Change the amount of data
        _data.clear();
        return _data;
    }
    while(_next_data.at(0) <= _time){
        _prev_data = _next_data;
        _next_data = readLine();
        if(!_next_data.size()){      // Last point reached
            _data = _prev_data;
            return _data;
        }
        if(_next_data.size() != _prev_data.size()){   // Change the amount of data
            _data.clear();
            return _data;
        }
    }
    // Linear interpolation of the data fields
    float factor = (_time - _prev_data.at(0)) / (_next_data.at(0) - _prev_data.at(0));
    _data.clear();
    for(i=0;i<_next_data.size();i++){
        float dData = _next_data.at(i) - _prev_data.at(i);
        _data.push_back( _prev_data.at(i) + factor*dData );
    }
    return _data;
}

bool LinearInterpolation::open(const char *data_file)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char msg[1024];
    // Clear data
    _prev_data.clear();
    _next_data.clear();
    // Open the file
    if(!data_file){
        S->addMessageF(3, "Invalid file string pointer.\n");
        return false;
    }
    if(!strlen(data_file)){
        S->addMessageF(3, "Empty data file path.\n");
        return false;
    }
    if(_data_file) fclose(_data_file); _data_file=0;
    _data_file = fopen(data_file, "r");
    if(!_data_file){
        sprintf(msg, "Failure while open \"%s\" data file.\n", data_file);
        S->addMessageF(3, msg);
        return false;
    }

    update(_time);
    return true;
}

std::deque<float> LinearInterpolation::readLine()
{
    std::deque<float> data;
    data.clear();
    float var;
    // Read lines unless reach a valid line
    while(true){
        // Read a line
        char Line[256];
        fgets( Line, 256*sizeof(char), _data_file);
        // Look for the start of line
        int LineStartChar=0;
        while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
            LineStartChar++;
        // Ensure that is not a commented or empty line
        if( (Line[LineStartChar] == '#') || (Line[LineStartChar] == '\n') )
            continue;
        // Ensure that end of file has not been reached
        if(Line[LineStartChar] == EOF)
            return data;
        // Erase initial spaces
        strcpy(Line, &Line[LineStartChar]);
        // Edit the line to get good formatted string. Replace , ; ( ) - \t by spaces
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
        // Erase initial spaces
        LineStartChar=0;
        while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
            LineStartChar++;
        strcpy(Line, &Line[LineStartChar]);
        // Erase final newline
        int LineEndChar=strlen(Line)-1;
        if(Line[LineEndChar] == '\n')
            strcpy(&Line[LineEndChar], "\0");
        // Erase eventual final space
        LineEndChar=strlen(Line)-1;
        if(Line[LineEndChar] == ' ')
            strcpy(&Line[LineEndChar], "\0");
        // Add data
        char *ReadPoint=Line-1;
        while(ReadPoint){
            ReadPoint++;
            // Read variable
            if( !sscanf(ReadPoint, "%g", &var) )
                break;
            // Add to the array
            data.push_back(var);
            // Look for next variable
            ReadPoint = strstr(ReadPoint," ");
        }
        break;
    }
    return data;
}

}}} // namespaces
