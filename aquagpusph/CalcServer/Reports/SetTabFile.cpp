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

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "SetTabFile.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

SetTabFile::SetTabFile(const std::string tool_name,
                       const std::string fields,
                       unsigned int first,
                       unsigned int n,
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
	} catch (std::invalid_argument e) {
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
	unsigned int i;

	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	_f.open(_output_file.c_str(), std::ios::out);

	Report::setup();

	// Write the header
	_f << "# Time ";
	std::vector<InputOutput::Variable*> vars = variables();
	for (i = _bounds.x; i < _bounds.y; i++) {
		for (auto var : vars) {
			_f << var->name() << "_" << i;
		}
	}
	_f << std::endl;
	_f.flush();
}

void
SetTabFile::print()
{
	std::vector<InputOutput::Variable*> vars = variables();
	for (unsigned int i = 0; i < bounds().y - bounds().x; i++) {
		for (unsigned int j = 0; j < vars.size(); j++) {
			InputOutput::ArrayVariable* var =
			    (InputOutput::ArrayVariable*)vars.at(j);
			const std::string type_name = var->type();
			if (!type_name.compare("int*")) {
				int* v = (int*)_data.at(j);
				_f << v[i] << " ";
			} else if (!type_name.compare("unsigned int*")) {
				unsigned int* v = (unsigned int*)_data.at(j);
				_f << v[i] << " ";
			} else if (!type_name.compare("float*")) {
				float* v = (float*)_data.at(j);
				_f << v[i] << " ";
			} else if (!type_name.compare("ivec*")) {
#ifdef HAVE_3D
				ivec* v = (ivec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
#else
				ivec* v = (ivec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
#endif // HAVE_3D
			} else if (!type_name.compare("ivec2*")) {
				ivec2* v = (ivec2*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
			} else if (!type_name.compare("ivec3*")) {
				ivec3* v = (ivec3*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ") ";
			} else if (!type_name.compare("ivec4*")) {
				ivec4* v = (ivec4*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
			} else if (!type_name.compare("uivec*")) {
#ifdef HAVE_3D
				uivec* v = (uivec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
#else
				uivec* v = (uivec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
#endif // HAVE_3D
			} else if (!type_name.compare("uivec2*")) {
				uivec2* v = (uivec2*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
			} else if (!type_name.compare("uivec3*")) {
				uivec3* v = (uivec3*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ") ";
			} else if (!type_name.compare("uivec4*")) {
				uivec4* v = (uivec4*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
			} else if (!type_name.compare("vec*")) {
#ifdef HAVE_3D
				vec* v = (vec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
#else
				vec* v = (vec*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
#endif // HAVE_3D
			} else if (!type_name.compare("vec2*")) {
				vec2* v = (vec2*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << ") ";
			} else if (!type_name.compare("vec3*")) {
				vec3* v = (vec3*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ") ";
			} else if (!type_name.compare("vec4*")) {
				vec4* v = (vec4*)_data.at(j);
				_f << "(" << v[i].x << "," << v[i].y << "," << v[i].z << ","
				   << v[i].w << ") ";
			}
		}
	}
	_f << std::endl;
	_f.flush();
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
SetTabFile::_execute(const std::vector<cl_event> events)
{
	if (!mustUpdate()) {
		return NULL;
	}

	unsigned int i, j;
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
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure releasing the downloading events.\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}
	return user_event;
}

std::vector<cl_event>
SetTabFile::download(std::vector<InputOutput::Variable*> vars)
{
	std::vector<cl_event> events; // vector storage is continuous memory
	size_t typesize, len;
	cl_int err_code;
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

	for (unsigned int i = 0; i < vars.size(); i++) {
		cl_event event;
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
