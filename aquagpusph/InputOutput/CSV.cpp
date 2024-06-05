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
#include "aquagpusph/CalcServer/CalcServer.hpp"
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
	auto vars = CalcServer::CalcServer::singleton()->variables();
	const char str = _sep == '"' ? '\'' : '"';
	const char* extensions[16] = { "_x",  "_y",  "_z",  "_w",
	                               "_yx", "_yy", "_yz", "_yw",
	                               "_zx", "_zy", "_zz", "_zw",
	                               "_wx", "_wy", "_wz", "_ww" };
	auto fields = simData().sets.at(setId())->outputFields();
	for (unsigned int i = 0; i < fields.size(); i++) {
		const auto var = vars->get(fields[i]);
		const auto n_subfields = vars->typeToN(var->type());
		if (n_subfields == 1)
			f << str << fields[i] << str;
		else {
			for (unsigned int j = 0; j < n_subfields; j++) {
				f << str << fields[i] << extensions[j] << str;
				if (j < n_subfields - 1)
					f << _sep;
			}
		}
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

size_t
CSV::readNParticles(std::ifstream& f)
{
	_has_header = false;
	size_t s = ASCII::readNParticles(f);
	_has_header = false;
	return s;
}

} // InputOutput::
} // Aqua::
