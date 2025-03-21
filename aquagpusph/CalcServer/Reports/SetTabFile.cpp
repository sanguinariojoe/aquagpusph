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
 * @brief Runtime output file.
 * (See Aqua::CalcServer::Reports::SetTabFile for details)
 */

#include <sstream>
#include <iomanip>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "SetTabFile.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

SetTabFile::SetTabFile(const std::string tool_name,
                       const std::string fields,
                       size_t first,
                       size_t n,
                       const std::string output_file,
                       unsigned int ipf,
                       float fps)
  : Report(tool_name, fields, ipf, fps)
  , _output_file("")
{
	_bounds.x = first;
	_bounds.y = first + n;
	try {
		unsigned int i = 0;
		_output_file = newFilePath(output_file, i, 1);
	} catch (std::invalid_argument& e) {
		std::ostringstream msg;
		_output_file = setStrConstantsCopy(output_file);
		msg << "Overwriting '" << _output_file << "'" << std::endl;
		LOG(L_WARNING, msg.str());
	}
}

SetTabFile::~SetTabFile()
{
	if (_f.is_open())
		_f.close();
	clearList(&_data);
}

void
SetTabFile::setup()
{
	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	_f.open(_output_file.c_str(), std::ios::out);
	CalcServer* C = CalcServer::singleton();
	if (C->device_addr_bits() == 64) {
		constexpr auto max_precision{
			std::numeric_limits<double>::digits10 + 1};
		_f << std::setprecision(max_precision);
	} else {
		constexpr auto max_precision{
			std::numeric_limits<float>::digits10 + 1};
		_f << std::setprecision(max_precision);
	}

	Report::setup();

	// Write the header
	_f << "# Time ";
	std::vector<InputOutput::Variable*> vars = variables();
	for (ulcl i = _bounds.x; i < _bounds.y; i++) {
		for (auto var : vars) {
			_f << var->name() << "_" << i << " ";
		}
	}
	_f << std::endl;
	_f.flush();
}

#define __ASCII_WRITE_SCALAR_BLOCK(TYPE)                                       \
	TYPE* v = (TYPE* )ptr;                                                     \
	_f << v[i];

#define __ASCII_WRITE_MATRIX_BLOCK(TYPE)                                       \
	TYPE* v = (TYPE *)ptr;                                                     \
	_f << '(';                                                                 \
	for (unsigned int k = 0; k < nc - 1; k++) {                                \
		_f << v[i].s[k] << ',';                                                \
	}                                                                          \
	_f << v[i].s[nc - 1] << ')';                                               \

#define __ASCII_WRITE_VEC_BLOCK(TYPE)                                          \
	if (nc == 2) {                                                             \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 2 )                                \
	} else if (nc == 3) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 3 )                                \
	} else if (nc == 4) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 4 )                                \
	} else if (nc == 8) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 8 )                                \
	} else if (nc == 16) {                                                     \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 16 )                               \
	}

void
SetTabFile::print()
{
	std::vector<InputOutput::Variable*> vars = variables();
	for (size_t i = 0; i < bounds().y - bounds().x; i++) {
		for (size_t j = 0; j < vars.size(); j++) {
			Aqua::InputOutput::ArrayVariable* var =
				(Aqua::InputOutput::ArrayVariable*)vars.at(j);
			const std::string type_name = var->type();
			const void* ptr = _data.at(j);
			const unsigned int nc = Aqua::InputOutput::Variables::typeToN(
				type_name);
			if (startswith(var->type(), "int")) {
				__ASCII_WRITE_SCALAR_BLOCK(icl)
			} else if (startswith(var->type(), "long")) {
				__ASCII_WRITE_SCALAR_BLOCK(lcl)
			} else if (startswith(var->type(), "unsigned int")) {
				__ASCII_WRITE_SCALAR_BLOCK(uicl)
			} else if (startswith(var->type(), "unsigned long")) {
				__ASCII_WRITE_SCALAR_BLOCK(ulcl)
			} else if (startswith(var->type(), "float")) {
				__ASCII_WRITE_SCALAR_BLOCK(fcl)
			} else if (startswith(var->type(), "double")) {
				__ASCII_WRITE_SCALAR_BLOCK(dcl)
			} else if (startswith(var->type(), "ivec")) {
				__ASCII_WRITE_VEC_BLOCK(ivec)
			} else if (startswith(var->type(), "lvec")) {
				__ASCII_WRITE_VEC_BLOCK(lvec)
			} else if (startswith(var->type(), "uivec")) {
				__ASCII_WRITE_VEC_BLOCK(uivec)
			} else if (startswith(var->type(), "ulvec")) {
				__ASCII_WRITE_VEC_BLOCK(ulvec)
			} else if (startswith(var->type(), "vec")) {
				__ASCII_WRITE_VEC_BLOCK(vec)
			} else if (startswith(var->type(), "dvec")) {
				__ASCII_WRITE_VEC_BLOCK(dvec)
			} else if (startswith(var->type(), "matrix")) {
				__ASCII_WRITE_VEC_BLOCK(vec)
			}
			_f << ' ';
		}
	}
	_f << std::endl;
	_f.flush();
	cl_int err_code;
	err_code = clSetUserEventStatus(getUserEvent(), CL_COMPLETE);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure setting as complete the tool \"") +
		name() + "\".");
	err_code = clReleaseEvent(getUserEvent());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure releasing the user event in tool \"") +
		name() + "\".");
}

/** @brief Callback called when all the variables required by
 * Aqua::CalcServer::SetTabFile have been downloaded.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::SetTabFile::print()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Screen
 * tool (or the inherited one)
 */
void CL_CALLBACK
settabfile_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (SetTabFile*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name() << "\" due to dependency errors."
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	tool->print();
}

cl_event
SetTabFile::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
	if (!mustUpdate()) {
		return NULL;
	}

	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	// Print the time instant
	_f << C->variables()->get("t")->asString() << " ";

	// Get the data to be printed
	std::vector<InputOutput::Variable*> vars = variables();
	for (auto var : vars) {
		if (var->type().find('*') == std::string::npos) {
			std::stringstream msg;
			msg << "The report \"" << name()
			    << "\" may not save scalar variables (\"" << var->name()
			    << "\")." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
		size_t typesize = C->variables()->typeToBytes(var->type());
		size_t len = var->size() / typesize;
		if (len < bounds().y) {
			std::stringstream msg;
			msg << "The report \"" << name() << "\" may not save field \""
			    << var->name() << "\" because is not long enough." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
	}
	std::vector<cl_event> new_events = download(vars);

	// Register the callback to work
	cl_event user_event = setCallback(new_events, settabfile_cb);
	for (auto event : new_events) {
		err_code = clReleaseEvent(event);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure releasing the downloading events in tool \"") +
		        name() + "\".");
	}
	return user_event;
}

std::vector<cl_event>
SetTabFile::download(std::vector<InputOutput::Variable*> vars)
{
	std::vector<cl_event> events; // vector storage is continuous memory
	size_t typesize, len;
	CalcServer* C = CalcServer::singleton();

	if (_data.size() != vars.size()) {
		for (auto var : vars) {
			typesize = C->variables()->typeToBytes(var->type());
			len = var->size() / typesize;
			if (len < bounds().y) {
				std::stringstream msg;
				msg << "The report \"" << name() << "\" may not save field \""
				    << var->name() << "\" because is not long enough."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				clearList(&_data);
				throw std::runtime_error("Invalid variable type");
			}
			void* store = malloc(typesize * (bounds().y - bounds().x));
			if (!store) {
				std::stringstream msg;
				msg << "Failure allocating "
				    << typesize * (bounds().y - bounds().x)
				    << " bytes for the field \"" << var->name() << "\"."
				    << std::endl;
				clearList(&_data);
				throw std::bad_alloc();
			}
			_data.push_back(store);
		}
	}

	for (size_t i = 0; i < vars.size(); i++) {
		cl_event event;
		typesize = C->variables()->typeToBytes(vars[i]->type());
		try {
			event = C->getUnsortedMem(vars[i]->name().c_str(),
			                          typesize * bounds().x,
			                          typesize * (bounds().y - bounds().x),
			                          _data[i]);
		} catch (...) {
			throw;
		}

		events.push_back(event);
	}

	return events;
}

void
SetTabFile::clearList(std::vector<void*>* data)
{
	for (auto d : *data) {
		if (d)
			free(d);
	}
	data->clear();
}

}
}
} // namespace
