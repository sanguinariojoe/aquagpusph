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

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <ScreenManager.h>
#include <TimeManager.h>
#include <CalcServer/Movements/C1Interpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

Poly::Poly(std::deque<vec2> p, std::deque<vec2> dp)
    : _x0(0.f)
	, _p(p)
	, _dp(dp)
{
    compute();
}

Poly::~Poly()
{
    _p.clear();
    _dp.clear();
    _k.clear();
}

float Poly::evaluate(float x)
{
    unsigned int i;
    float y = 0.f;
    for(i=0;i<_k.size();i++){
        y += _k.at(i) * pow(x - _x0, (float)(i));
    }
    return y;
}

float Poly::derivate(float x)
{
    unsigned int i;
    float y = 0.f;
    for(i=1;i<_k.size();i++){
        y += i * _k.at(i) * pow(x - _x0, (float)(i-1));
    }
    return y;
}

void Poly::compute()
{
    unsigned int i,j;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    _k.clear();
    // Get the problem dimensions
    unsigned int N = _p.size() + _dp.size();
    if(!N){
	    S->addMessageF(3, "No conditions imposed, so the curve cannot be computed.\n");
	    return;
    }
    // If no point are imposed one DoF will remain ever, so it is an error
    if(!_p.size()){
	    S->addMessageF(3, "No points has been imposed, and therefore the problem is not closed.\n");
	    return;
    }
    _x0 = _p.at(0).x;
    for(i=1;i<_p.size();i++){
        if(_x0 > _p.at(i).x) _x0 = _p.at(i).x;
    }
    Eigen::MatrixXf A(N,N);
    Eigen::VectorXf b(N);
    // Impose the points
    for(i=0;i<_p.size();i++){
        float x = _p.at(i).x - _x0;
        float y = _p.at(i).y;
        b(i) = y;
        for(j=0;j<N;j++){
            A(i,j) = pow(x,(float)(j));
        }
    }
    // Impose the derivatives
    for(i=0;i<_dp.size();i++){
        float x = _dp.at(i).x;
        float y = _dp.at(i).y;
        b(_p.size() + i) = y;
        A(_p.size() + i,0) = 0.f;
        for(j=1;j<N;j++){
            A(_p.size() + i,j) = j * pow(x,(float)(j-1));
        }
    }
    // Solve the linear system
    Eigen::VectorXf k;
    k = A.lu().solve(b);
    for(i=0;i<N;i++){
        _k.push_back(k(i));
    }
}

C1Interpolation::C1Interpolation(const char *data_file)
	: _data_file(NULL)
	, _time(0.f)
	, _poly_time(0.f)
{
	// Open the file (if possible)
	if(data_file)
	    open(data_file);
}

C1Interpolation::~C1Interpolation()
{
    unsigned int i;
	if(_data_file) fclose(_data_file); _data_file=0;
	for(i=0;i<_poly.size();i++){
        delete _poly.at(i); _poly.at(i) = NULL;
	}
    _poly.clear();
	_data.clear();
}

std::deque<float> C1Interpolation::update(float t)
{
	unsigned int i,j;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	if(!_data_file){
        _data.clear();
	    return _data;
	}
	// Rewind the file if the time is lower than the last time
	if(t < _time){
	    rewind(_data_file);
        for(i=0;i<_poly.size();i++){
            delete _poly.at(i); _poly.at(i) = NULL;
        }
        _poly.clear();
	    _poly_time = 0.f;
	}
	_time = t;
	// Test if the computed curve still becomes valid
	if(_poly.size() && t <= _poly_time){
        _data.clear();
	    _data.push_back(t);
        for(i=0;i<_poly.size();i++){
            _data.push_back(_poly.at(i)->evaluate(t));
        }
        return _data;
	}
    // Read a new data block to interpolate
    std::deque<float> *data = NULL;
    unsigned int dataDim = 0;
	while(_poly_time <= _time){
	    std::deque<float> line = readLine();
	    // Test if the last point has been reached
	    if(!line.size()){
	        S->addMessageF(2, "Final time step reached.\n");
	        return _data;
	    }
	    // Test if the file has suddently changed the number of fields
	    if(!data){
            dataDim = line.size();
	        data = new std::deque<float>[dataDim];
	    }
        if(dataDim != line.size()){
	        S->addMessageF(3, "The number of data fields has suddently changed.\n");
            return _data;
        }
        // All seems gone right, so we can append the new data
        for(i=0;i<line.size();i++){
            data[i].push_back(line.at(i));
        }
        _poly_time = line.at(0);
	}

	// -----------------------------
	// Create the new curves
	// -----------------------------
    std::deque<Poly*> poly = _poly;
    _poly.clear();
    std::deque<vec2> p;
    std::deque<vec2> d;
    for(i=1;i<dataDim;i++){
        p.clear();
        d.clear();
        // Start setting the initial derivative and point from the previous curve (if exist)
        if(poly.size()){
            if(poly.size() != dataDim-1){
                S->addMessageF(3, "The number of data fields has changed for the new data block.\n");
                return _data;
            }
            vec2 point;
            point.x = _time;
            point.y = poly.at(i-1)->evaluate(_time);
            vec2 derivative;
            derivative.x = _time;
            derivative.y = poly.at(i-1)->derivate(_time);
            p.push_back(point);
            d.push_back(derivative);
        }
        // Add the new readed points
        for(j=0;j<data[i].size();j++){
            vec2 point;
            point.x = data[0].at(j);
            point.y = data[i].at(j);
            p.push_back(point);
        }
        // Create the new curve
        _poly.push_back(new Poly(p,d));
    }
    // Clean up
    for(i=0;i<poly.size();i++){
        delete poly.at(i); poly.at(i) = NULL;
    }
    poly.clear();
    for(i=0;i<dataDim;i++){
        data[i].clear();
    }
    delete[] data; data = NULL;
    // Compute the requested data
    _data.clear();
    _data.push_back(_time);
    for(i=0;i<_poly.size();i++){
        _data.push_back(_poly.at(i)->evaluate(_time));
    }
    return _data;
}

std::deque<float> C1Interpolation::derivative()
{
    unsigned int i;
    std::deque<float> data;
	if(!_poly.size()){
	    return data;
	}
    data.push_back(_time);
    for(i=0;i<_poly.size();i++){
        data.push_back(_poly.at(i)->derivate(_time));
    }
    return data;
}

bool C1Interpolation::open(const char *data_file)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[1024];
	// Clear previous data
	_poly_time = 0.f;
	// Open the file
	if(!data_file){
	    S->addMessageF(3, "Null string pointer received.\n");
	    return false;
	}
	if(!strlen(data_file)){
	    S->addMessageF(3, "Empty data file path.\n");
	    return false;
	}
	if(_data_file) fclose(_data_file); _data_file=0;
	_data_file = fopen(data_file, "r");
	if(!_data_file){
	    sprintf(msg, "Failure opening \"%s\" data file.\n", data_file);
	    S->addMessageF(3, msg);
	    return false;
	}
	// Go to selected time
	update(_time);
	return true;
}

std::deque<float> C1Interpolation::readLine()
{
	std::deque<float> data;
	data.clear();
	float var;
	// Read lines unless we reach a valid one
	while(true){
	    // Read a line
	    char Line[256];
	    fgets( Line, 256*sizeof(char), _data_file);
	    // Look for the start of it
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
	    // Edit the line to get a conveniently formatted string.
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
	    // Erase initial spaces
	    LineStartChar=0;
	    while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
	        LineStartChar++;
	    strcpy(Line, &Line[LineStartChar]);
	    // Erase final newline char
	    int LineEndChar=strlen(Line)-1;
	    if(Line[LineEndChar] == '\n')
	        strcpy(&Line[LineEndChar], "\0");
	    // Erase trailing spaces
	    LineEndChar=strlen(Line)-1;
	    if(Line[LineEndChar] == ' ')
	        strcpy(&Line[LineEndChar], "\0");
	    // Add the data
	    char *ReadPoint=Line-1;
	    while(ReadPoint){
	        ReadPoint++;
	        // Read the variable
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
