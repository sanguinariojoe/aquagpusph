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
 * @brief Check a condition to enable/disable all the tools in its scope
 * (See Aqua::CalcServer::Conditional and Aqua::CalcServer::End for details)
 */

#include <sstream>

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "CalcServer.hpp"
#include "Conditional.hpp"

namespace Aqua {
namespace CalcServer {

Conditional::Conditional(const std::string name,
                         const std::string condition,
                         bool once)
  : ScalarExpression(name, condition, "int", once)
  , _ending_tool(NULL)
  , _result(true)
{
}

Conditional::~Conditional() {}

void
Conditional::setup()
{
	std::vector<Tool*> tools = CalcServer::singleton()->tools();

	// Get the next tool in case the condition is not fulfilled
	int i = id_in_pipeline();
	if (i < 0) {
		std::ostringstream msg;
		msg << "Tool \"" << name() << "\" cannot be found in the pipeline"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid tool");
	}

	int scope = 1;
	while ((unsigned int)i < tools.size() - 1) {
		i++;
		scope += tools.at(i)->scope_modifier();
		if (!scope)
			break;
	}
	if (scope > 0) {
		std::ostringstream msg;
		msg << "Unclosed scope for Tool \"" << name() << "\". Add and End tool"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Unbalanced scope");
	}
	// We cannot use next_tool() attribute since that tool has not been setup
	// yet
	if ((unsigned int)i == tools.size() - 1)
		_ending_tool = NULL;
	else
		_ending_tool = tools.at(i + 1);

	ScalarExpression::setup();
}

Tool*
Conditional::next_tool()
{
	cl_event event = getEvent();
	cl_int err_code = clWaitForEvents(1, &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure waiting for \"") + name() +
	                       "\" conditional event.");
	if (_result)
		return Tool::next_tool();
	return _ending_tool;
}

void
Conditional::_solve()
{
	ScalarExpression::_solve();
	// Check the result
	memcpy(&_result, getValue(), sizeof(int));
}

While::While(const std::string name, const std::string condition, bool once)
  : Conditional(name, condition, once)
{
}

While::~While() {}

If::If(const std::string name, const std::string condition, bool once)
  : Conditional(name, condition, once)
{
}

If::~If() {}

Tool*
If::next_tool()
{
	Tool* next_tool = Conditional::next_tool();
	// There are several possibilities here:
	// 1.- We have evaluated the condition, and it has been true, so we want
	//     to skip the next evaluation (_result = false), when the associated
	//     End gives back us the control
	// 2.- We have evaluated the condition, and it has been false, we are
	//     jumping out of our scope. Thus next time we has control we want to
	//     evaluate the condition again (_result = true)
	// 3.- End has gave back control to us, so _result = false (see point 1).
	//     this situation is exactly the same than point 2.
	// Therefore, all the cases are covered just simply swaping _result value
	_result = !_result;
	return next_tool;
}

void
If::_solve()
{
	// Execute the tool just if _result is true. Otherwise is an End tool which
	// has gave back the control to us
	// See If::next_tool()
	if (_result)
		Conditional::_solve();
}

End::End(const std::string name, bool once)
  : Tool(name, once)
{
}

End::~End() {}

void
End::setup()
{
	std::vector<Tool*> tools = CalcServer::singleton()->tools();

	// Locate the opening scope tool (which will be always the next tool on the
	// pipeline)
	int i = id_in_pipeline();
	if (i < 0) {
		std::ostringstream msg;
		msg << "Tool \"" << name() << "\" cannot be found in the pipeline"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid tool");
	}

	int scope = 1;
	while (i > 0) {
		i--;
		scope -= tools.at(i)->scope_modifier();
		if (!scope)
			break;
	}
	if (scope > 0) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" cannot be associated to any scope opening tool (If/While)"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Unbalanced scope");
	}
	next_tool(tools.at(i));
}

}
} // namespaces
