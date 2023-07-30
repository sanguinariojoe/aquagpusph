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
 * @brief link-list based neighbours location algorithm.
 * (See Aqua::CalcServer::LinkList for details)
 * @note Hardcoded versions of the files CalcServer/LinkList.cl.in and
 * CalcServer/LinkList.hcl.in are internally included as a text array.
 */

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/LinkList.h>
#include <algorithm>

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/LinkList.hcl"
#include "CalcServer/LinkList.cl"
#endif
std::string LINKLIST_INC = xxd2string(LinkList_hcl_in, LinkList_hcl_in_len);
std::string LINKLIST_SRC = xxd2string(LinkList_cl_in, LinkList_cl_in_len);

LinkList::LinkList(const std::string tool_name,
                   const std::string input,
                   bool once)
  : Tool(tool_name, once)
  , _input_name(input)
  , _cell_length(0.f)
  , _min_pos(NULL)
  , _max_pos(NULL)
  , _ihoc(NULL)
  , _ihoc_lws(0)
  , _ihoc_gws(0)
  , _icell(NULL)
  , _icell_lws(0)
  , _icell_gws(0)
  , _ll(NULL)
  , _ll_lws(0)
  , _ll_gws(0)
{
	std::stringstream min_pos_name;
	min_pos_name << tool_name << "->Min. Pos.";
	std::string min_pos_op = "c.x = (a.x < b.x) ? a.x : b.x;\nc.y = (a.y < "
	                         "b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z < "
	                         "b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
	_min_pos = new Reduction(
	    min_pos_name.str(), input, "r_min", min_pos_op, "VEC_INFINITY");
	std::stringstream max_pos_name;
	max_pos_name << tool_name << "->Max. Pos.";
	std::string max_pos_op = "c.x = (a.x > b.x) ? a.x : b.x;\nc.y = (a.y > "
	                         "b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z > "
	                         "b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
	_max_pos = new Reduction(
	    max_pos_name.str(), input, "r_max", max_pos_op, "-VEC_INFINITY");
	std::stringstream sort_name;
	sort_name << tool_name << "->Radix-Sort";
	_sort = new RadixSort(sort_name.str());
}

LinkList::~LinkList()
{
	if (_min_pos)
		delete _min_pos;
	_min_pos = NULL;
	if (_max_pos)
		delete _max_pos;
	_max_pos = NULL;
	if (_sort)
		delete _sort;
	_sort = NULL;
	if (_ihoc)
		clReleaseKernel(_ihoc);
	_ihoc = NULL;
	if (_icell)
		clReleaseKernel(_icell);
	_icell = NULL;
	if (_ll)
		clReleaseKernel(_ll);
	_ll = NULL;
	for (auto arg : _ihoc_args) {
		free(arg);
	}
	_ihoc_args.clear();
	for (auto arg : _icell_args) {
		free(arg);
	}
	_icell_args.clear();
	for (auto arg : _ll_args) {
		free(arg);
	}
	_ll_args.clear();
}

void
LinkList::setup()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();

	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();

	// Setup the reduction tools
	_min_pos->setup();
	_max_pos->setup();

	// Compute the cells length
	InputOutput::Variable* s = vars->get("support");
	InputOutput::Variable* h = vars->get("h");
	_cell_length = *(float*)s->get() * *(float*)h->get();

	// Setup the kernels
	setupOpenCL();

	// Setup the radix-sort
	_sort->setup();

	// _input_name at front and icell at back on forward purpose!
	std::vector<std::string> deps = { _input_name, "ihoc", "icell" };
	setDependencies(deps);
}

cl_event
LinkList::_execute(const std::vector<cl_event> events_prior)
{
	cl_int err_code;
	cl_event event, event_wait;
	std::vector<cl_event> events;
	CalcServer* C = CalcServer::singleton();

	// Reduction steps to find maximum and minimum position
	_min_pos->execute();
	_max_pos->execute();

	// We should refresh the events adding the new one (we can just keep the
	// outdated ones, which are already retained). The new events existence are
	// granted while we don't set new events to those variables, which would not
	// gonna happen until we returns the new event herein.
	// The new event comes from the first dependency (see setup())
	std::copy(
	    events_prior.begin(), events_prior.end(), std::back_inserter(events));
	events.push_back(getDependencies().front()->getEvent());

	// Compute the number of cells, and eventually allocate memory for ihoc
	nCells();
	allocate();

	// Check the validity of the variables
	setVariables();

	// Compute the cell of each particle
	cl_uint num_events_in_wait_list = events.size();
	const cl_event* event_wait_list = events.size() ? events.data() : NULL;
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _icell,
	                                  1,
	                                  NULL,
	                                  &_icell_gws,
	                                  &_icell_lws,
	                                  num_events_in_wait_list,
	                                  event_wait_list,
	                                  &event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure executing \"iCell\" from tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	getDependencies().front()->setEvent(event);
	getDependencies().back()->setEvent(event);
	err_code = clReleaseEvent(event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing transactional \"iCell\" event from tool \""
		    << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// Sort the particles from the cells
	_sort->execute();

	// Now our transactional event is the one coming from sorting algorithm
	// Such a new event can be taken from the last dependency (see setup())
	// This transactional event SHALL NOT BE RELEASED. It is automagically
	// destroyed when no more variables use it
	event_wait = getDependencies().back()->getEvent();

	// Compute the head of cells
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _ihoc,
	                                  1,
	                                  NULL,
	                                  &_ihoc_gws,
	                                  &_ihoc_lws,
	                                  1,
	                                  &event_wait,
	                                  &event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure executing \"iHoc\" from tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	event_wait = event; // This new transactional event should be released

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _ll,
	                                  1,
	                                  NULL,
	                                  &_ll_gws,
	                                  &_ll_lws,
	                                  1,
	                                  &event_wait,
	                                  &event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure executing \"linkList\" from tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clReleaseEvent(event_wait);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing transactional \"linkList\" event from tool \""
		    << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	return event;
}

void
LinkList::setupOpenCL()
{
	unsigned int i;
	uivec4 n_cells;
	unsigned int n_radix, N;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Create a header for the source code where the operation will be placed
	std::ostringstream source;
	source << LINKLIST_INC << LINKLIST_SRC;

	std::vector<cl_kernel> kernels =
	    compile(source.str(), { "iHoc", "iCell", "linkList" });
	_ihoc = kernels.at(0);
	_icell = kernels.at(1);
	_ll = kernels.at(2);

	err_code = clGetKernelWorkGroupInfo(_ihoc,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_ihoc_lws,
	                                    NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size (\"iHoc\").\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	if (_ihoc_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"iHoc\".\n");
		std::stringstream msg;
		msg << "\t" << _ihoc_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	n_cells = *(uivec4*)vars->get("n_cells")->get();
	_ihoc_gws = roundUp(n_cells.w, _ihoc_lws);
	const char* _ihoc_vars[3] = { "ihoc", "N", "n_cells" };
	for (i = 0; i < 3; i++) {
		err_code = clSetKernelArg(_ihoc,
		                          i,
		                          vars->get(_ihoc_vars[i])->typesize(),
		                          vars->get(_ihoc_vars[i])->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure sending \"" << _ihoc_vars[i]
			    << "\" argument to \"iHoc\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		_ihoc_args.push_back(malloc(vars->get(_ihoc_vars[i])->typesize()));
		memcpy(_ihoc_args.at(i),
		       vars->get(_ihoc_vars[i])->get(),
		       vars->get(_ihoc_vars[i])->typesize());
	}

	err_code = clGetKernelWorkGroupInfo(_icell,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_icell_lws,
	                                    NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size (\"iCell\").\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	if (_icell_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"iCell\".\n");
		std::stringstream msg;
		msg << "\t" << _icell_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	n_radix = *(unsigned int*)vars->get("n_radix")->get();
	_icell_gws = roundUp(n_radix, _icell_lws);
	const char* _icell_vars[8] = {
		"icell", _input_name.c_str(), "N", "n_radix",
		"r_min", "support",           "h", "n_cells"
	};
	for (i = 0; i < 8; i++) {
		err_code = clSetKernelArg(_icell,
		                          i,
		                          vars->get(_icell_vars[i])->typesize(),
		                          vars->get(_icell_vars[i])->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure sending \"" << _icell_vars[i]
			    << "\" argument to \"iCell\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		_icell_args.push_back(malloc(vars->get(_icell_vars[i])->typesize()));
		memcpy(_icell_args.at(i),
		       vars->get(_icell_vars[i])->get(),
		       vars->get(_icell_vars[i])->typesize());
	}

	err_code = clGetKernelWorkGroupInfo(_ll,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_ll_lws,
	                                    NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size (\"linkList\").\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	if (_ll_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"linkList\".\n");
		std::stringstream msg;
		msg << "\t" << _ll_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	N = *(unsigned int*)vars->get("N")->get();
	_ll_gws = roundUp(N, _ll_lws);
	const char* _ll_vars[3] = { "icell", "ihoc", "N" };
	for (i = 0; i < 3; i++) {
		err_code = clSetKernelArg(_ll,
		                          i,
		                          vars->get(_ll_vars[i])->typesize(),
		                          vars->get(_ll_vars[i])->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure sending \"" << _ll_vars[i]
			    << "\" argument to \"iCell\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		_ll_args.push_back(malloc(vars->get(_ll_vars[i])->typesize()));
		memcpy(_ll_args.at(i),
		       vars->get(_ll_vars[i])->get(),
		       vars->get(_ll_vars[i])->typesize());
	}
}

void
LinkList::nCells()
{
	vec pos_min, pos_max;
	InputOutput::Variables* vars = CalcServer::singleton()->variables();

	if (!_cell_length) {
		std::stringstream msg;
		msg << "Zero cell length detected in the tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid number of cells");
	}

	pos_min = *(vec*)vars->get("r_min")->get();
	pos_max = *(vec*)vars->get("r_max")->get();

	_n_cells.x = (unsigned int)((pos_max.x - pos_min.x) / _cell_length) + 6;
	_n_cells.y = (unsigned int)((pos_max.y - pos_min.y) / _cell_length) + 6;
#ifdef HAVE_3D
	_n_cells.z = (unsigned int)((pos_max.z - pos_min.z) / _cell_length) + 6;
#else
	_n_cells.z = 1;
#endif
	_n_cells.w = _n_cells.x * _n_cells.y * _n_cells.z;
}

void
LinkList::allocate()
{
	uivec4 n_cells;
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	if (vars->get("n_cells")->type().compare("uivec4")) {
		std::stringstream msg;
		msg << "\"n_cells\" has and invalid type for \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"n_cells\" type is \""
		    << vars->get("n_cells")->type()
		    << "\", while \"uivec4\" was expected" << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid n_cells type");
	}

	n_cells = *(uivec4*)vars->get("n_cells")->get();

	if (_n_cells.w <= n_cells.w) {
		n_cells.x = _n_cells.x;
		n_cells.y = _n_cells.y;
		n_cells.z = _n_cells.z;
		vars->get("n_cells")->set(&n_cells);
		return;
	}

	// We have no alternative, we must sync here
	InputOutput::Variable* ihoc_var = vars->get("ihoc");
	event = ihoc_var->getEvent();
	clWaitForEvents(1, &event);
	cl_mem mem = *(cl_mem*)ihoc_var->get();
	if (mem)
		clReleaseMemObject(mem);
	mem = NULL;

	mem = clCreateBuffer(C->context(),
	                     CL_MEM_READ_WRITE,
	                     _n_cells.w * sizeof(unsigned int),
	                     NULL,
	                     &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure allocating device memory in the tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL allocation error");
	}

	vars->get("n_cells")->set(&_n_cells);
	ihoc_var->set(&mem);
	_ihoc_gws = roundUp(_n_cells.w, _ihoc_lws);
}

void
LinkList::setVariables()
{
	unsigned int i;
	cl_int err_code;
	InputOutput::Variables* vars = CalcServer::singleton()->variables();

	const char* _ihoc_vars[3] = { "ihoc", "N", "n_cells" };
	for (i = 0; i < 3; i++) {
		InputOutput::Variable* var = vars->get(_ihoc_vars[i]);
		if (!memcmp(var->get(), _ihoc_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_ihoc, i, var->typesize(), var->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the variable \"" << _ihoc_vars[i]
			    << "\" to the tool \"" << name() << "\" (\"iHoc\")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		memcpy(_ihoc_args.at(i), var->get(), var->typesize());
	}

	const char* _icell_vars[8] = {
		"icell", _input_name.c_str(), "N", "n_radix",
		"r_min", "support",           "h", "n_cells"
	};
	for (i = 0; i < 8; i++) {
		InputOutput::Variable* var = vars->get(_icell_vars[i]);
		if (!memcmp(var->get(), _icell_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_icell, i, var->typesize(), var->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the variable \"" << _ihoc_vars[i]
			    << "\" to the tool \"" << name() << "\" (\"iCell\")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		memcpy(_icell_args.at(i), var->get(), var->typesize());
	}

	const char* _ll_vars[3] = { "icell", "ihoc", "N" };
	for (i = 0; i < 3; i++) {
		InputOutput::Variable* var = vars->get(_ll_vars[i]);
		if (!memcmp(var->get(), _ll_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_ll, i, var->typesize(), var->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the variable \"" << _ihoc_vars[i]
			    << "\" to the tool \"" << name() << "\" (\"linkList\")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		memcpy(_ll_args.at(i), var->get(), var->typesize());
	}
}

}
} // namespace
