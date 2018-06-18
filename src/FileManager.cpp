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
#include <InputOutput/FastASCII.h>
#ifdef HAVE_VTK
    #include <InputOutput/VTK.h>
#endif // HAVE_VTK

namespace Aqua{ namespace InputOutput{

FileManager::FileManager()
    : _state(NULL)
    , _log(NULL)
{
    inputFile("Input.xml");
    _state = new State();
    _log = new Log();
}

FileManager::~FileManager()
{
    unsigned int i;

    if(_state) delete _state; _state = NULL;
    if(_log) delete _log; _log = NULL;

    for(auto i = _loaders.begin(); i != _loaders.end(); i++) {
        delete *i;
    }
    for(auto i = _savers.begin(); i != _savers.end(); i++) {
        delete *i;
    }
}

void FileManager::inputFile(std::string path)
{
    _in_file = path;
}

FILE* FileManager::logFile()
{
    if(_log)
        return _log->fileHandler();
    return NULL;
}

CalcServer::CalcServer* FileManager::load()
{
    unsigned int i, n=0;
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
        delete C;
        return NULL;
    }

    // Now we can build the loaders/savers
    for(i = 0; i < P->sets.size(); i++){
        if(!strcmp(P->sets.at(i)->inputFormat(), "ASCII")){
            ASCII *loader = new ASCII(n, P->sets.at(i)->n(), i);
            if(!loader){
                delete C;
                return NULL;
            }
            _loaders.push_back((Particles*)loader);
        }
        else if(!strcmp(P->sets.at(i)->inputFormat(), "FastASCII")){
            FastASCII *loader = new FastASCII(n, P->sets.at(i)->n(), i);
            if(!loader){
                delete C;
                return NULL;
            }
            _loaders.push_back((Particles*)loader);
        }
        else if(!strcmp(P->sets.at(i)->inputFormat(), "VTK")){
            #ifdef HAVE_VTK
                VTK *loader = new VTK(n, P->sets.at(i)->n(), i);
                if(!loader){
                    delete C;
                    return NULL;
                }
                _loaders.push_back((Particles*)loader);
            #else
                S->addMessageF(L_ERROR, "AQUAgpusph has been compiled without VTK format.\n");
                return NULL;
            #endif // HAVE_VTK
        }
        else{
            std::ostringstream msg;
            msg << "Unknow \"" << P->sets.at(i)->inputFormat()
                << "\" input file format" << std::endl;
            S->addMessageF(L_ERROR, msg.str());
            delete C;
            return NULL;
        }
        if(!strcmp(P->sets.at(i)->outputFormat(), "ASCII")){
            ASCII *saver = new ASCII(n, P->sets.at(i)->n(), i);
            if(!saver){
                delete C;
                return NULL;
            }
            _savers.push_back((Particles*)saver);
        }
        else if(!strcmp(P->sets.at(i)->outputFormat(), "VTK")){
            #ifdef HAVE_VTK
                VTK *saver = new VTK(n, P->sets.at(i)->n(), i);
                if(!saver){
                    delete C;
                    return NULL;
                }
                _savers.push_back((Particles*)saver);
            #else
                S->addMessageF(L_ERROR, "AQUAgpusph has been compiled without VTK format.\n");
                return NULL;
            #endif // HAVE_VTK
        }
        else{
            std::ostringstream msg;
            msg << "Unknow \"" << P->sets.at(i)->outputFormat()
                << "\" input file format" << std::endl;
            S->addMessageF(L_ERROR, msg.str());
            delete C;
            return NULL;
        }
        n += P->sets.at(i)->n();
    }

    // Execute the loaders
    for(auto loader = _loaders.begin(); loader != _loaders.end(); loader++) {
        if((*loader)->load()) {
            delete C;
            return NULL;
        }
    }

    return C;
}

bool FileManager::save()
{
    unsigned int i;

    // Execute the savers
    for(auto saver = _savers.begin(); saver != _savers.end(); saver++) {
        if((*saver)->save())
            return true;
    }

    // Save the XML definition file
    if(_state->save()){
        return true;
    }

    return false;
}

std::string FileManager::file(unsigned int i){
    return _savers.at(i)->file();
}


}}  // namespace
