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
// Linear algebra tools
// ----------------------------------------------------------------------------
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

// ----------------------------------------------------------------------------
// Include the Screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/C1Interpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

Poly::Poly(std::deque<vec2> p, std::deque<vec2> dp)
    : mX0(0.f)
	, mP(p)
	, mD(dp)
{
    // Compute the data
    compute();
}

Poly::~Poly()
{
    mP.clear();
    mD.clear();
    mK.clear();
}

float Poly::evaluate(float x)
{
    unsigned int i;
    float y = 0.f;
    for(i=0;i<mK.size();i++){
        y += mK.at(i) * pow(x - mX0, (float)(i));
    }
    return y;
}

float Poly::derivate(float x)
{
    unsigned int i;
    float y = 0.f;
    for(i=1;i<mK.size();i++){
        y += i * mK.at(i) * pow(x - mX0, (float)(i-1));
    }
    return y;
}

void Poly::compute()
{
    unsigned int i,j;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    mK.clear();
    // Get the problem dimensions
    unsigned int N = mP.size() + mD.size();
    if(!N){
	    S->addMessage(3, "(Poly::compute): No conditions imposed, so the curve can't be computed.\n");
	    return;
    }
    // If no any point is imposed one DoF will remain ever, so it is an error
    if(!mP.size()){
	    S->addMessage(3, "(Poly::compute): No points imposed, and therefore the problem is not closed.\n");
	    return;
    }
    mX0 = mP.at(0).x;
    for(i=1;i<mP.size();i++){
        if(mX0 > mP.at(i).x) mX0 = mP.at(i).x;
    }
    Eigen::MatrixXf A(N,N);
    Eigen::VectorXf b(N);
    // Impose the points
    for(i=0;i<mP.size();i++){
        float x = mP.at(i).x - mX0;
        float y = mP.at(i).y;
        b(i) = y;
        for(j=0;j<N;j++){
            A(i,j) = pow(x,(float)(j));
        }
    }
    // Impose the derivatives
    for(i=0;i<mD.size();i++){
        float x = mD.at(i).x;
        float y = mD.at(i).y;
        b(mP.size() + i) = y;
        A(mP.size() + i,0) = 0.f;
        for(j=1;j<N;j++){
            A(mP.size() + i,j) = j * pow(x,(float)(j-1));
        }
    }
    // Solve the linear system
    Eigen::VectorXf k;
    k = A.lu().solve(b);
    for(i=0;i<N;i++){
        mK.push_back(k(i));
    }
}

C1Interpolation::C1Interpolation(const char *dataFile)
	: mDataFile(NULL)
	, _time(0.f)
	, mPolyTime(0.f)
{
	// Open the file (if possible)
	if(dataFile)
	    open(dataFile);
}

C1Interpolation::~C1Interpolation()
{
    unsigned int i;
	if(mDataFile) fclose(mDataFile); mDataFile=0;
	for(i=0;i<mPoly.size();i++){
        delete mPoly.at(i); mPoly.at(i) = NULL;
	}
    mPoly.clear();
	mData.clear();
}

std::deque<float> C1Interpolation::update(float t)
{
	unsigned int i,j;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	if(!mDataFile){
        mData.clear();
	    return mData;
	}
	// Rewind the file if the time is lower than the last time
	if(t < _time){
	    rewind(mDataFile);
        for(i=0;i<mPoly.size();i++){
            delete mPoly.at(i); mPoly.at(i) = NULL;
        }
        mPoly.clear();
	    mPolyTime = 0.f;
	}
	_time = t;
	// Test if the computed curve still being valid
	if(mPoly.size() && t <= mPolyTime){
        mData.clear();
	    mData.push_back(t);
        for(i=0;i<mPoly.size();i++){
            mData.push_back(mPoly.at(i)->evaluate(t));
        }
        return mData;
	}
    // Read a new data block
    std::deque<float> *data = NULL, empty; // We don't know the number of fields yet
    unsigned int dataDim = 0;
	while(mPolyTime <= _time){
	    std::deque<float> line = readLine();
	    // Test if the last point has been reached
	    if(!line.size()){
	        S->addMessage(2, "(C1Interpolation::update): Final time step reached.\n");
	        return mData;
	    }
	    // Test if the file has suddently changed the number of fields
	    if(!data){
            dataDim = line.size();
	        data = new std::deque<float>[dataDim];
	    }
        if(dataDim != line.size()){
	        S->addMessage(3, "(C1Interpolation::update): The number of data fields has suddently changed.\n");
            return mData;
        }
        // All seems gone right, so we can append the new data
        for(i=0;i<line.size();i++){
            data[i].push_back(line.at(i));
        }
        mPolyTime = line.at(0);
	}

	// -----------------------------
	// Create the new curves
	// -----------------------------
    std::deque<Poly*> poly = mPoly;
    mPoly.clear();
    std::deque<vec2> p;
    std::deque<vec2> d;
    for(i=1;i<dataDim;i++){
        p.clear();
        d.clear();
        // Start setting the initial derivative and point from the previous curve (if exist)
        if(poly.size()){
            if(poly.size() != dataDim-1){
                S->addMessage(3, "(C1Interpolation::update): The number of data fields has changed for the new data block.\n");
                return mData;
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
        // Add the readed points
        for(j=0;j<data[i].size();j++){
            vec2 point;
            point.x = data[0].at(j);
            point.y = data[i].at(j);
            p.push_back(point);
        }
        // Create the new curve
        mPoly.push_back(new Poly(p,d));
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
    mData.clear();
    mData.push_back(_time);
    for(i=0;i<mPoly.size();i++){
        mData.push_back(mPoly.at(i)->evaluate(_time));
    }
    return mData;
}

std::deque<float> C1Interpolation::derivative()
{
    unsigned int i;
    std::deque<float> data;
	if(!mPoly.size()){
	    return data;
	}
    data.push_back(_time);
    for(i=0;i<mPoly.size();i++){
        data.push_back(mPoly.at(i)->derivate(_time));
    }
    return data;
}

bool C1Interpolation::open(const char *dataFile)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[1024];
	// Clear previous data
	mPolyTime = 0.f;
	// Open the file
	if(!dataFile){
	    S->addMessage(3, "(C1Interpolation::open): Null string pointer received.\n");
	    return false;
	}
	if(!strlen(dataFile)){
	    S->addMessage(3, "(C1Interpolation::open): Empty data file path.\n");
	    return false;
	}
	if(mDataFile) fclose(mDataFile); mDataFile=0;
	mDataFile = fopen(dataFile, "r");
	if(!mDataFile){
	    sprintf(msg, "(C1Interpolation::open): Can't open \"%s\" data file.\n", dataFile);
	    S->addMessage(3, msg);
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
	// Read lines unless reach a valid line
	while(true){
	    // Read a line
	    char Line[256];
	    fgets( Line, 256*sizeof(char), mDataFile);
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
	        if( !sscanf(ReadPoint, "%f", &var) )
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
