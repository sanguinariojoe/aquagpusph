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

#include <algorithm>
#include <InputOutput/Logger.h>
#include <CalcServer/Reports/Screen.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Screen::Screen(const std::string tool_name,
               const std::string fields,
               const std::string color,
               bool bold)
    : Report(tool_name, fields)
    , _color(color)
    , _bold(bold)
{
}

Screen::~Screen()
{
}

void Screen::setup()
{
    unsigned int i;

    std::ostringstream msg;
    msg << "Loading the report \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    // Set the color in lowercase
    std::transform(_color.begin(), _color.end(), _color.begin(), ::tolower);

    Report::setup();
}

void Screen::_execute()
{
    InputOutput::Logger::singleton()->writeReport(data(), _color, _bold);
}

}}} // namespace
