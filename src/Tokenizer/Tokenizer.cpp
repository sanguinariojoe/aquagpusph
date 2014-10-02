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
 * @brief Math expression evaluator.
 * (See Aqua::Tokenizer for details)
 */

#include <matheval.h>
#include <stdlib.h>

#include <Tokenizer/Tokenizer.h>
#include <ScreenManager.h>

using namespace std;

namespace Aqua{

Tokenizer::Tokenizer()
{
    // Register default variables
    registerVariable("pi", M_PI);
    registerVariable("e", M_E);
}

Tokenizer::~Tokenizer()
{
    _variables.clear();
}

bool Tokenizer::registerVariable(const char* name, float value)
{
    // Look for the variable in order to know if it already exist
    if(isVariable(name)){
        // The variable already exist, let we take it and modify
        _variables[name] = value;
        return true;
    }
    _variables.insert(make_pair(name, value));
    return false;
}

bool Tokenizer::unregisterVariable(const char* name)
{
    if(!isVariable(name)){
        return false;
    }
    map<string, float>::iterator var = _variables.find(name);
    _variables.erase(var);
    return true;
}

void Tokenizer::clearVariables()
{
    _variables.clear();
    defaultVariables();
}

bool Tokenizer::isVariable(const char* name)
{
    map<string, float>::iterator var = _variables.find(name);
    if(var != _variables.end())
        return true;
    return false;
}

float Tokenizer::variable(const char* name)
{
    if(!isVariable(name)){
        return 0.f;
    }
    return _variables[name];
}

float Tokenizer::solve(const char* eq, bool *error)
{
    Aqua::InputOutput::ScreenManager *S;
    char msg[1024], *test;
    void *f;
    char **names;
    double *values;
    int n, i;
    float result;
    S = Aqua::InputOutput::ScreenManager::singleton();

    if(error)
        *error = false;

    // First test if the equation is directly a number
    result = strtof(eq, &test);
    if(test != eq){
        // We have been able to convert it, but we must test that there are
        // no remaining data to compute
        if(!strlen(test)){
            return result;
        }
    }

    f = evaluator_create((char*)eq);
    if(!f){
        S->addMessageF(3, "Invalid math expression to evaluate:\n");
        sprintf(msg, "\t\"%s\"\n", eq);
        S->addMessage(0, msg);
        if(error)
            *error = true;
        return 0.f;
    }

    evaluator_get_variables(f, &names, &n);
    values = new double[n];
    if(!values){
        S->addMessageF(3, "Failure allocating memory for the variables.\n");
        evaluator_destroy(f);
        if(error)
            *error = true;
        return 0.f;
    }
    for(i=0;i<n;i++){
        if(!isVariable(names[i])){
            S->addMessageF(3, "Impossible to evaluate a variable\n");
            sprintf(msg,
                    "\t\"%s\" variable has not been registered\n",
                    names[i]);
            S->addMessage(0, msg);
            evaluator_destroy(f);
            delete[] values;
            if(error)
                *error = true;
            return 0.f;
        }
        values[i] = variable(names[i]);
    }

    result = evaluator_evaluate(f, n, names, values);

    evaluator_destroy(f);
    delete[] values;
    return result;
}

void Tokenizer::defaultVariables()
{
    // Pi and e are registered variables out of the box
}

}  // namespace
