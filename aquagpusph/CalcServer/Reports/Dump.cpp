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
 * @brief Runtime array dumping to an output file.
 * (See Aqua::CalcServer::Reports::Dump for details)
 */

#include <sstream>
#include <fstream>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "Dump.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

Dump::Dump(const std::string tool_name,
           const std::string fields,
           const std::string output_file,
           bool binary)
  : Report(tool_name, fields, 1, 0.f)
  , _output_file(output_file)
  , _binary(binary)
  , _n(0)
{
}

Dump::~Dump()
{
	for (auto d : _data) {
		if (d)
			free(d);
	}
}

void
Dump::setup()
{
	std::ostringstream msg;
	msg << "Loading the dump report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Report::setup();
}

#define __ASCII_WRITE_SCALAR_BLOCK(TYPE)                                       \
	TYPE* v = (TYPE* )ptr;                                                     \
	f << v[i];

#define __ASCII_WRITE_MATRIX_BLOCK(TYPE)                                       \
	TYPE* v = (TYPE *)ptr;                                                     \
	for (auto k = 0; k < nc - 1; k++) {                                        \
		f << v[i].s[k] << ' ';                                                 \
	}                                                                          \
	f << v[i].s[nc - 1];                                                       \

#define __ASCII_WRITE_VEC_BLOCK(TYPE)                                          \
	if (nc == 2) {                                                             \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 2 )                                \
	} else if (nc == 3) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 3 )                                \
	} else if (nc == 4) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 4 )                                \
	} else if (nc == 8) {                                                      \
		__ASCII_WRITE_MATRIX_BLOCK( TYPE ## 8 )                                \
	}

void
Dump::print()
{
	unsigned int findex=0;
	const auto fpath = newFilePath(_output_file, findex);
	std::ios_base::openmode mode = _binary ?
		std::ios::out | std::ios::binary:
		std::ios::out;
	std::ofstream f(fpath, mode);
	std::vector<InputOutput::Variable*> vars = variables();

	if (_binary) {
		for (size_t i = 0; i < vars.size(); i++) {
			InputOutput::ArrayVariable* var =
			    (InputOutput::ArrayVariable*)vars[i];
			size_t datasize = var->size();
			f.write(reinterpret_cast<const char*>(_data[i]), datasize);
		}
	} else {
		for (size_t i = 0; i < _n; i++) {
			for (size_t j = 0; j < vars.size(); j++) {
				InputOutput::ArrayVariable* var =
					(InputOutput::ArrayVariable*)vars.at(j);
				const std::string type_name = var->type();
				const void* ptr = _data.at(j);
				const auto nc = Aqua::InputOutput::Variables::typeToN(
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
					__ASCII_WRITE_MATRIX_BLOCK(matrix)
				}
				if (j < vars.size() - 1)
					f << ' ';
			}
			f << std::endl;
			f.flush();
		}
	}
	f.close();
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
 * Aqua::CalcServer::Dump have been downloaded.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::Dump::print()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Screen
 * tool (or the inherited one)
 */
void CL_CALLBACK
dumpfile_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (Dump*)user_data;
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
Dump::_execute(const std::vector<cl_event> events)
{
	// We are always dumping
	// if (!mustUpdate()) {
	// 	return NULL;
	// }

	size_t n = 0;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	// Get the data to be printed
	std::vector<InputOutput::Variable*> vars = variables();
	if (_data.size() < vars.size())
		_data.resize(vars.size(), NULL);
	for (auto var : vars) {
		if (var->type().find('*') == std::string::npos) {
			std::stringstream msg;
			msg << "The dump report \"" << name()
			    << "\" may not save scalar variables (\"" << var->name()
			    << "\")." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
		size_t typesize = C->variables()->typeToBytes(var->type());
		size_t len = var->size() / typesize;
		if (!len) {
			std::stringstream msg;
			msg << "The dump report \"" << name()
			    << "\" cannot print the array \"" << var->name()
				<< "\", because it has null length." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable length");
		}
		if (!n)
			n = len;
		if (len > n) {
			std::stringstream msg;
			msg << "The dump report \"" << name()
			    << "\" cannot save field \""
			    << var->name() << "\" because is too long" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable length");
		}
	}
	if (n > _n) {
		for (auto d : _data) {
			if (d)
				free(d);
			d = NULL;
		}
	}
	_n = n;
	if (!_n) {
		std::stringstream msg;
		msg << "The dump report \"" << name()
			<< "\" has nothing to do" << std::endl;
		LOG(L_WARNING, msg.str());
		return NULL;
	}
	std::vector<cl_event> new_events = download(vars);

	// Register the callback to work
	cl_event user_event = setCallback(new_events, dumpfile_cb);
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
Dump::download(std::vector<InputOutput::Variable*> vars)
{
	std::vector<cl_event> events;
	size_t typesize, len;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	for (unsigned int i = 0; i < _data.size(); i++) {
		InputOutput::ArrayVariable* var = (InputOutput::ArrayVariable*)vars[i];
		size_t datasize = var->size();
		// Allocate memory if needed
		if (!_data[i]) {
			void* store = malloc(datasize);
			if (!store) {
				std::stringstream msg;
				msg << "Failure allocating " << datasize
					<< " bytes for the field \"" << var->name() << "\"."
					<< std::endl;
				throw std::bad_alloc();
			}
			_data[i] = store;
		}

		cl_event event, wait_event = var->getWritingEvent();
		err_code = clEnqueueReadBuffer(C->command_queue(),
		                               *((cl_mem*)var->get()),
		                               CL_FALSE,
		                               0,
		                               datasize,
		                               _data[i],
		                               1,
		                               &wait_event,
		                               &event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure receiving the variable \"") +
		                   var->name() + "\"");
		events.push_back(event);
	}

	return events;
}

}
}
} // namespace
