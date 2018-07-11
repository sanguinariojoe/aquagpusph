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

#include <Tokenizer/Tokenizer.h>
#include <InputOutput/Logger.h>

using namespace std;

namespace Aqua{

Tokenizer::Tokenizer()
{
    struct lconv *lc;
    char *s;

    // Set the decimal-point character (which is depending on the locale)
    s = setlocale(LC_NUMERIC, NULL);
    if(strcmp(s, "C")){
        std::ostringstream msg;
        msg << "\"" << s << "\" numeric locale found" << std::endl;
        LOG(L_INFO, msg.str());
        LOG0(L_DEBUG, "\tIt is replaced by \"C\"\n");
        setlocale(LC_NUMERIC, "C");
    }
    lc = localeconv();
    s = lc->decimal_point;
    if(strcmp(s, ".")){
        std::ostringstream msg;
        msg << "\"" << s << "\" decimal point character found" << std::endl;
        LOG(L_WARNING, msg.str());
        LOG0(L_DEBUG, "\tIt is replaced by \".\"\n");
        lc->decimal_point = ".";
    }
    s = lc->thousands_sep;
    if(strcmp(s, "")){
        std::ostringstream msg;
        msg << "\"" << s << "\" thousands separator character found" << std::endl;
        LOG(L_WARNING, msg.str());
        LOG0(L_DEBUG, "\tIt is removed\n");
        lc->thousands_sep = "";
    }
}

Tokenizer::~Tokenizer()
{
    p.ClearConst();
}

bool Tokenizer::registerVariable(const std::string name, float value)
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

bool Tokenizer::isVariable(const std::string name)
{
    mu::valmap_type cmap = p.GetConst();
    if (cmap.size())
    {
        mu::valmap_type::const_iterator item = cmap.begin();
        for (; item!=cmap.end(); ++item){
            if(!name.compare(item->first.c_str())){
                return true;
            }
        }
    }
    return false;
}

float Tokenizer::variable(const std::string name)
{
    if(!isVariable(name)){
        return 0.f;
    }
    mu::valmap_type cmap = p.GetConst();
    return (float)cmap[name.c_str()];
}


float Tokenizer::solve(const std::string eq, bool *error)
{
    float result;

    if(error)
        *error = false;

    // First try a straight number conversion
    try {
        std::string::size_type sz;
        result = std::stof(eq, &sz);
        if (sz == eq.size()) {
            // There is remaining content
            return result;
        }
    }
    catch(...){
        // No possible conversion (by a variety of errors), just proceed with
        // the parser
    }

    // No way, let's evaluate it as a expression
    p.SetExpr(eq);
    try
    {
        result = (float)p.Eval();
    }
    catch(mu::Parser::exception_type &e)
    {
        std::ostringstream msg;
        msg << "Error evaluating \"" << e.GetExpr() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\t" << e.GetMsg() << std::endl;
        LOG0(L_DEBUG, msg.str());
        msg.str("");
        msg << "\tToken " << e.GetToken()
            << " in position " << e.GetPos() << std::endl;
        LOG0(L_DEBUG, msg.str());
    }

    return result;
}

void Tokenizer::defaultVariables()
{
    // Pi and e are registered variables out of the box
}

}  // namespace
