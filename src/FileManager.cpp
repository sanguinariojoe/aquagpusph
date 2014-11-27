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
 * @brief Input and output files managing.
 * (See Aqua::InputOutput::FileManager for details)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <FileManager.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <InputOutput/ASCII.h>
#ifdef HAVE_VTK
    #include <InputOutput/VTK.h>
#endif // HAVE_VTK

namespace Aqua{ namespace InputOutput{

FileManager::FileManager()
    : _in_file(NULL)
    , _log(NULL)
{
    inputFile("Input.xml");
    _log = new Log();
    _state = new State();
}

FileManager::~FileManager()
{
    unsigned int i;

    if(_in_file)
        delete[] _in_file;
    _in_file = NULL;
    if(_state)
        delete _state;
    _state = NULL;
    if(_log)
        delete _log;
    _log = NULL;

    for(i=0; i<_loaders.size(); i++){
        delete _loaders.at(i);
        _loaders.at(i) = NULL;
    }
    _loaders.clear();

    for(i=0; i<_savers.size(); i++){
        delete _savers.at(i);
        _savers.at(i) = NULL;
    }
    _savers.clear();
}

void FileManager::inputFile(const char* path)
{
    size_t len;

    if(_in_file)
        delete[] _in_file;
    _in_file = NULL;

    if(!path)
        return;

    len = strlen(path) + 1;
    _in_file = new char[len];
    strcpy(_in_file, path);
}

FILE* FileManager::logFile()
{
    if(_log)
        return _log->fileHandler();
}

CalcServer::CalcServer* FileManager::load()
{
    unsigned int i, n=0;
    char msg[1024];
    ProblemSetup *P = ProblemSetup::singleton();
    ScreenManager *S = ScreenManager::singleton();

    // Load the XML definition file
    if(_state->load()){
        return NULL;
    }

    // Setup the problem setup
    if(P->perform()) {
        return NULL;
    }

    // Build the calculation server
    CalcServer::CalcServer *C = new CalcServer::CalcServer();
    if(!C)
        return NULL;
    if(C->setup()){
        delete C; C = NULL;
        return NULL;
    }

    // Now we can build the loaders/savers
    for(i = 0; i < P->sets.size(); i++){
        if(!strcmp(P->sets.at(i)->inputFormat(), "ASCII")){
            ASCII *loader = new ASCII(n, P->sets.at(i)->n(), i);
            if(!loader){
                delete C; C = NULL;
                return NULL;
            }
            _loaders.push_back((Particles*)loader);
        }
        else if(!strcmp(P->sets.at(i)->inputFormat(), "VTK")){
            #ifdef HAVE_VTK
                VTK *loader = new VTK(n, P->sets.at(i)->n(), i);
                if(!loader){
                    delete C; C = NULL;
                    return NULL;
                }
                _loaders.push_back((Particles*)loader);
            #else
                S->addMessageF(3, "AQUAgpusph has been compiled without VTK format.\n");
                return NULL;
            #endif // HAVE_VTK
        }
        else{
            sprintf(msg,
                    "Unknow \"%s\" input file format.\n",
                    P->sets.at(i)->inputFormat());
            S->addMessageF(3, msg);
            delete C; C = NULL;
            return NULL;
        }
        if(!strcmp(P->sets.at(i)->outputFormat(), "ASCII")){
            ASCII *saver = new ASCII(n, P->sets.at(i)->n(), i);
            if(!saver){
                delete C; C = NULL;
                return NULL;
            }
            _savers.push_back((Particles*)saver);
        }
        else if(!strcmp(P->sets.at(i)->outputFormat(), "VTK")){
            #ifdef HAVE_VTK
                VTK *saver = new VTK(n, P->sets.at(i)->n(), i);
                if(!saver){
                    delete C; C = NULL;
                    return NULL;
                }
                _savers.push_back((Particles*)saver);
            #else
                S->addMessageF(3, "AQUAgpusph has been compiled without VTK format.\n");
                return NULL;
            #endif // HAVE_VTK
        }
        else{
            sprintf(msg,
                    "Unknow \"%s\" output file format.\n",
                    P->sets.at(i)->outputFormat());
            S->addMessageF(3, msg);
            delete C; C = NULL;
            return NULL;
        }
        n += P->sets.at(i)->n();
    }

    // Execute the loaders
    for(i = 0; i < _loaders.size(); i++){
        if(_loaders.at(i)->load())
            return NULL;
    }

    return C;
}

bool FileManager::save()
{
    unsigned int i;

    // Execute the savers
    for(i = 0; i < _savers.size(); i++){
        if(_savers.at(i)->save())
            return true;
    }

    // Save the XML definition file
    if(_state->save()){
        return true;
    }

    return false;
}

const char* FileManager::file(unsigned int i){
    return _savers.at(i)->file();
}


}}  // namespace
