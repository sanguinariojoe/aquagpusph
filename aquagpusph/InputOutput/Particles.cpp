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
 * @brief Particles files manager.
 * (See Aqua::InputOutput::Particles for details)
 */

#include <string>
#include <sstream>

#include "Particles.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {
namespace InputOutput {

Particles::Particles(ProblemSetup& sim_data,
                     unsigned int iset,
                     size_t first,
                     size_t n)
  : _sim_data(sim_data)
  , _iset(iset)
  , _time(0.f)
  , _user_event(NULL)
{
	_bounds.x = first;
	_bounds.y = first + n;
}

Particles::~Particles()
{
	if (_user_event)
		clReleaseEvent(_user_event);
	for (auto const& mem : _data)
		free(mem.second);
}

/** @brief Callback called when all the fields have been downloaded.
 *
 * This function is just redirecting the work to
 * Aqua::InputOutput::Particles::print_file()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or more dependencies failed.
 * @param user_data A casted pointer to the Aqua::InputOutput::Particles
 * tool (or the inherited one)
 */
void CL_CALLBACK
particles_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (Particles*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		LOG(L_WARNING, "Skipping the file saving due to errors.\n");
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	tool->print_file();
}

void
Particles::save(float t)
{
	std::vector<std::string> fields =
	    simData().sets.at(setId())->outputFields();
	if (!fields.size()) {
		LOG(L_ERROR, "0 fields were set to be saved into the file.\n");
		throw std::runtime_error("No fields have been marked to be saved");
	}

	// Just one instance at a time
	waitForSavers();
	_time = t;

	cl_int err_code;
	auto C = CalcServer::CalcServer::singleton();
	auto trigger = download(fields);
	auto event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure creating the user event \n.");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	_user_event = event;

	// So it is time to register our callback on our trigger
	err_code = clRetainEvent(event);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure retaining the user event\n.");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clSetEventCallback(trigger, CL_COMPLETE, particles_cb, this);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure registering the printing callback.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
}

void
Particles::loadDefault()
{
	unsigned int i;
	cl_int err_code;
	ArrayVariable* var;
	cl_mem mem;
	auto C = CalcServer::CalcServer::singleton();

	size_t n = bounds().y - bounds().x;
	uicl* iset = new uicl[n];
	// Let's assume we have 64bits. Otherwise we will recast this during the
	// assingments
	ulcl* id = new ulcl[n];
	size_t id_size = sizeof(ulcl);

	if (!iset || !id) {
		LOG(L_ERROR, "Failure allocating memory.\n");
	}

	for (i = 0; i < n; i++) {
		iset[i] = setId();
		if (C->device_addr_bits() == 64)
			id[i] = bounds().x + i;
		else {
			((uicl*)id)[i] = (uicl)(bounds().x + i);
			id_size = sizeof(uicl);
		}
	}

	Variables* vars = C->variables();
	var = (ArrayVariable*)vars->get("iset");
	mem = *(cl_mem*)var->get();
	err_code = clEnqueueWriteBuffer(C->command_queue(),
	                                mem,
	                                CL_TRUE,
	                                sizeof(uicl) * bounds().x,
	                                sizeof(uicl) * n,
	                                iset,
	                                0,
	                                NULL,
	                                NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending variable \"iset\" to the server.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	var = (ArrayVariable*)vars->get("id");
	mem = *(cl_mem*)var->get();
	err_code = clEnqueueWriteBuffer(C->command_queue(),
	                                mem,
	                                CL_TRUE,
	                                id_size * bounds().x,
	                                id_size * n,
	                                id,
	                                0,
	                                NULL,
	                                NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending variable \"id\" to the server.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	var = (ArrayVariable*)vars->get("id_sorted");
	mem = *(cl_mem*)var->get();
	err_code = clEnqueueWriteBuffer(C->command_queue(),
	                                mem,
	                                CL_TRUE,
	                                id_size * bounds().x,
	                                id_size * n,
	                                id,
	                                0,
	                                NULL,
	                                NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending variable \"id_sorted\" to the server.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	var = (ArrayVariable*)vars->get("id_unsorted");
	mem = *(cl_mem*)var->get();
	err_code = clEnqueueWriteBuffer(C->command_queue(),
	                                mem,
	                                CL_TRUE,
	                                id_size * bounds().x,
	                                id_size * n,
	                                id,
	                                0,
	                                NULL,
	                                NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR,
		    "Failure sending variable \"id_unsorted\" to the server.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}

	delete[] iset;
	iset = NULL;
	delete[] id;
	id = NULL;
}

unsigned int
Particles::file(const std::string basename,
                unsigned int startindex,
                unsigned int digits)
{
	unsigned int i = startindex;
	try {
		file(newFilePath(basename, i, digits));
	} catch (std::invalid_argument& e) {
		std::ostringstream msg;
		msg << "It is forbidden to overwrite particles output file '"
		    << setStrConstantsCopy(basename) << "'" << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
	return i;
}

cl_event
Particles::download(std::vector<std::string> fields)
{
	std::vector<cl_event> events;
	size_t typesize, len;
	cl_int err_code;
	auto C = CalcServer::CalcServer::singleton();
	Variables* vars = C->variables();

	{
		std::ostringstream msg;
		msg << "Reading particles set " << setId() << "..." << std::endl;
		LOG(L_INFO, msg.str());
	}

	// We are doing all this process on a new parallel command queue
	C->command_queue(CalcServer::CalcServer::cmd_queue::cmd_queue_new);

	for (auto field : fields) {
		if (!vars->get(field)) {
			std::ostringstream msg;
			msg << "Can't download undeclared variable \"" << field << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable");
		}
		if (!vars->get(field)->isArray()) {
			std::ostringstream msg;
			msg << "Variable \"" << field << "\" is a scalar." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
		ArrayVariable* var = (ArrayVariable*)vars->get(field);
		typesize = vars->typeToBytes(var->type());
		len = var->size() / typesize;
		if (len < bounds().y) {
			std::ostringstream msg;
			msg << "Variable \"" << field << "\" is not long enough."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			msg.str("");
			msg << "length = " << bounds().y << "is required, but just " << len
			    << " components are available." << std::endl;
			LOG0(L_DEBUG, msg.str());
			throw std::runtime_error("Invalid variable length");
		}
		if (_data.find(field) == _data.end()) {
			void* store = malloc(typesize * (bounds().y - bounds().x));
			if (!store) {
				std::ostringstream msg;
				msg << "Failure allocating "
				    << typesize * (bounds().y - bounds().x)
				    << "bytes for variable \"" << field << "\"." << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::bad_alloc();
			}
			_data.insert(std::make_pair(field, store));
		}

		auto event = C->getUnsortedMem(var->name().c_str(),
		                               typesize * bounds().x,
		                               typesize * (bounds().y - bounds().x),
		                               _data.at(field));
		events.push_back(event);
	}

	// Join all the events together
	cl_event trigger;
	err_code = clEnqueueMarkerWithWaitList(
	    C->command_queue(), events.size(), events.data(), &trigger);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure setting the trigger for data saving.\n");
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// And we do not want everybody waiting for us, so we create another
	// parallel command queue to continue the work, just in case
	C->command_queue(CalcServer::CalcServer::cmd_queue::cmd_queue_new);

	return trigger;
}

}
} // namespace
