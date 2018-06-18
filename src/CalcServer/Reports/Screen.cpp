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
 * @brief On screen runtime output.
 * (See Aqua::CalcServer::Reports::Screen for details)
 */

#include <CalcServer/Reports/Screen.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Screen::Screen(const char* tool_name,
               const char* fields,
               const char* color,
               bool bold)
    : Report(tool_name, fields)
    , _color(NULL)
    , _bold(bold)
{
    _color = new char[strlen(color) + 1];
    strcpy(_color, color);
}

Screen::~Screen()
{
    if(_color) delete[] _color; _color=NULL;
}

bool Screen::setup()
{
    unsigned int i;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the report \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    // Set the color in lowercase
    for(i = 0; i < strlen(_color); i++){
        _color[i] = tolower(_color[i]);
    }

    if(Report::setup()){
        return true;
    }

    return false;
}

bool Screen::_execute()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    S->writeReport(data(), _color, _bold);
    return false;
}

}}} // namespace
