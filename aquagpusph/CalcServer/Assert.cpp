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
 * @brief Check that a condition holds true, or trhow a fatal error otherwise.
 * (See Aqua::CalcServer::Assert for details)
 */

#include <math.h>

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "CalcServer.hpp"
#include "Assert.hpp"

namespace Aqua {
namespace CalcServer {

Assert::Assert(const std::string name, const std::string condition, bool once)
  : ScalarExpression(name, condition, "int", once)
{
}

Assert::~Assert() {}

void
Assert::_solve()
{
	int result;
	ScalarExpression::_solve();
	// Check the result
	memcpy(&result, getValue(), sizeof(int));
	if (result == 0) {
		std::stringstream msg;
		msg << "Assertion error. The expression \"" << getExpression()
		    << "\" on tool \"" << name() << "\" is false" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Assertion error");
	}
}

}
} // namespaces
