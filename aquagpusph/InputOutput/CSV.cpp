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
 * @brief Particles plain text data files loader/saver.
 * (See Aqua::InputOutput::FastASCII for details)
 */

#include <string>
#include "CSV.hpp"
#include "aquagpusph/ProblemSetup.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {
namespace InputOutput {

CSV::CSV(ProblemSetup& sim_data,
                     unsigned int iset,
                     size_t first,
                     size_t n,
                     const std::string file_ext,
                     const char sep)
  : FastASCII(sim_data, iset, first, n, file_ext)
  , _sep(sep)
  , _has_header(false)
{
}

CSV::~CSV() {}

void
CSV::load()
{
	_has_header = false;
	ASCII::load();
}

void
CSV::print_file()
{
	ASCII::print_file(_sep, _sep);
}

void
CSV::print_header(std::ofstream& f) const
{
	const char str = _sep == '"' ? '\'' : '"'; 
	auto fields = simData().sets.at(setId())->outputFields();
	for (unsigned int i = 0; i < fields.size(); i++) {
		f << str << fields[i] << str;
		if (i < fields.size() - 1)
			f << _sep;
	}
	f << std::endl;
	f.flush();
}

void
CSV::formatLine(std::string& l)
{
	trim(l);

	if (_has_header)
		return;

	// Split the string and check all the fields. If a field is not a number,
	// we discard the line as a header
	for (auto field : split(l, _sep)) {
		try {
			std::stod(field);
		} catch(...) {
			l = "";
			break;
		}
	}

	_has_header = true;
}

} // InputOutput::
} // Aqua::
