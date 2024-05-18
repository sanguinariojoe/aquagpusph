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
#include "FastASCII.hpp"
#include "Logger.hpp"
#include "aquagpusph/ProblemSetup.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {
namespace InputOutput {

FastASCII::FastASCII(ProblemSetup& sim_data,
                     unsigned int iset,
                     size_t first,
                     size_t n,
                     const std::string file_ext)
  : ASCII(sim_data, iset, first, n, file_ext)
{
}

FastASCII::~FastASCII() {}

std::string
FastASCII::readField(const std::string field,
                     const std::string line,
                     size_t index,
                     void* data)
{
	unsigned int i;
	Variables* vars = CalcServer::CalcServer::singleton()->variables();
	ArrayVariable* var = (ArrayVariable*)vars->get(field);

	// Extract the variable type data
	unsigned int n = vars->typeToN(var->type());
	size_t type_size = vars->typeToBytes(var->type());
	std::string type = trimCopy(var->type());
	if (type.back() == '*') {
		type.pop_back();
	}

	// Point to the chunk of data to become read
	void* ptr = (void*)((char*)data + type_size * index);

	std::string remaining = line;
	for (i = 0; i < n; i++) {
		std::string::size_type end_pos;
		try {
			if (!type.compare("unsigned int") ||
			    startswith(type, "uivec")) {
				uicl val =
				    (uicl)std::stoul(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(uicl));
			} else if (!type.compare("unsigned long") ||
			           startswith(type, "ulvec")) {
				ulcl val =
				    (ulcl)std::stoull(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(ulcl));
			} else if (!type.compare("int") ||
			           startswith(type, "ivec")) {
				icl val = (icl)std::stoi(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(icl));
			} else if (!type.compare("long") ||
			           startswith(type, "lvec")) {
				lcl val = (lcl)std::stoi(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(lcl));
			} else if (!type.compare("float") ||
			           startswith(type, "vec")) {
				fcl val = (fcl)std::stof(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(fcl));
			} else {
				dcl val = std::stod(remaining, &end_pos);
				memcpy(ptr, &val, sizeof(dcl));
			}
		} catch (const std::invalid_argument& e) {
			std::ostringstream msg;
			msg << "Cannot extract a number from \"" << remaining
			    << "\" string." << std::endl;
			LOG(L_ERROR, msg.str());
			throw;
		} catch (const std::out_of_range& e) {
			std::ostringstream msg;
			msg << "The number extracted from \"" << remaining
			    << "\" string overflows \"" << type << "\" type." << std::endl;
			LOG(L_ERROR, msg.str());
			throw;
		}

		// Go to the next field, we already asserted that there are fields
		// enough, so we don't need to care about that
		ptr = ((char*)ptr) + type_size / n;
		remaining = remaining.substr(end_pos);
		if (remaining.find(',') != std::string::npos)
			remaining = remaining.substr(remaining.find(',') + 1);
	}
	return remaining;
}

}
} // namespace
