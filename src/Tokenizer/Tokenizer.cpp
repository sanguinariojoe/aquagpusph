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

// #include <matheval.h>
#include <stdlib.h>

#include <Tokenizer/Tokenizer.h>
#include <ScreenManager.h>

using namespace std;

namespace Aqua{

Tokenizer::Tokenizer()
{
    struct lconv *lc;
    char *s, msg[256];
    Aqua::InputOutput::ScreenManager *S = Aqua::InputOutput::ScreenManager::singleton();

    // Set the decimal-point character (which is depending on the locale)
    s = setlocale(LC_NUMERIC, NULL);
    if(strcmp(s, "C")){
        sprintf(msg, "\"%s\" numeric locale found\n", s);
        S->addMessageF(1, msg);
        S->addMessage(0, "\tIt is replaced by \"C\"\n");
        setlocale(LC_NUMERIC, "C");
    }
    lc = localeconv();
    s = lc->decimal_point;
    if(strcmp(s, ".")){
        sprintf(msg, "\"%s\" decimal point character found\n", s);
        S->addMessageF(2, msg);
        S->addMessage(0, "\tIt is replaced by \".\"\n");
        lc->decimal_point = ".";
    }
    s = lc->thousands_sep;
    if(strcmp(s, "")){
        sprintf(msg, "\"%s\" thousands separator character found\n", s);
        S->addMessageF(2, msg);
        S->addMessage(0, "\tIt is removed\n");
        lc->thousands_sep = "";
    }
}

Tokenizer::~Tokenizer()
{
    p.ClearConst();
}

bool Tokenizer::registerVariable(const char* name, float value)
{
    bool overwritten = false;
    // Look for the variable in order to know if it already exist
    if(isVariable(name)){
        // The variable already exist
        overwritten = true;
    }
    p.DefineConst(name, (mu::value_type)value);
    return  overwritten;
}

void Tokenizer::clearVariables()
{
    p.ClearVar();
    defaultVariables();
}

bool Tokenizer::isVariable(const char* name)
{
    mu::valmap_type cmap = p.GetConst();
    if (cmap.size())
    {
        mu::valmap_type::const_iterator item = cmap.begin();
        for (; item!=cmap.end(); ++item){
            if(!strcmp(name, item->first.c_str())){
                return true;
            }
        }
    }
    return false;
}

float Tokenizer::variable(const char* name)
{
    if(!isVariable(name)){
        return 0.f;
    }
    mu::valmap_type cmap = p.GetConst();
    return (float)cmap[name];
}


float Tokenizer::solve(const char* eq, bool *error)
{
    Aqua::InputOutput::ScreenManager *S;
    char msg[1024], *test;
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

    // No way, let's evaluate the new expression
    p.SetExpr(eq);
    try
    {
        result = (float)p.Eval();
    }
    catch(mu::Parser::exception_type &e)
    {
        sprintf(msg, "Error evaluating \"%s\"\n", e.GetExpr().c_str());
        S->addMessageF(3, msg);
        sprintf(msg, "\t%s\n", e.GetMsg().c_str());
        S->addMessage(0, msg);
        sprintf(msg, "\tToken %s in position %d\n", e.GetToken().c_str(),
                                                    e.GetPos());
        S->addMessage(0, msg);
    }

    return result;
}

void Tokenizer::defaultVariables()
{
    // Pi and e are registered variables out of the box
}

}  // namespace
