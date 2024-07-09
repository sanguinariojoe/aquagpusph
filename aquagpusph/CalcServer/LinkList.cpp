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

#include <sstream>
#include <algorithm>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "CalcServer.hpp"
#include "LinkList.hpp"
#include "RadixSort.hpp"
#include "Sort.hpp"
#include "SetScalar.hpp"
#include "Kernel.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/LinkList.hcl"
#include "aquagpusph/CalcServer/LinkList.cl"
#endif
std::string LINKLIST_INC = xxd2string(LinkList_hcl_in, LinkList_hcl_in_len);
std::string LINKLIST_SRC = xxd2string(LinkList_cl_in, LinkList_cl_in_len);

LinkList::LinkList(const std::string tool_name,
                   const std::string input,
                   const std::string input_min,
                   const std::string input_max,
                   const std::string ihoc,
                   const std::string icell,
                   const std::string n_cells,
                   const std::string permutations,
                   const std::string inv_permutations,
                   bool recompute_grid,
                   const std::string sorter,
                   bool once)
  : Tool(tool_name, once)
  , _input_name(input)
  , _input_min_name(input_min)
  , _input_max_name(input_max)
  , _ihoc_name(ihoc)
  , _icell_name(icell)
  , _n_cells_name(n_cells)
  , _cell_length(0.f)
  , _min_pos(NULL)
  , _max_pos(NULL)
  , _recompute_minmax(recompute_grid)
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
	if (recompute_grid) {
		_min_pos = new Reduction(tool_name + "->Min. Pos.",
		                         input,
		                         "r_min",
		                         "c = min(a, b);",
		                         "VEC_INFINITY");
		_max_pos = new Reduction(tool_name + "->Max. Pos.",
		                         input,
		                         "r_max",
		                         "c = max(a, b);",
		                         "-VEC_INFINITY");
	}
	if (toLowerCopy(sorter) == "radix-sort") {
		_sort = new RadixSort(tool_name + "->Radix-Sort",
		                      icell, permutations, inv_permutations);
	} else if (toLowerCopy(sorter) == "bitonic") {
		_sort = new Sort(tool_name + "->Bitonic-Sort",
		                 icell, permutations, inv_permutations);
	} else {
		LOG(L_ERROR, std::string("Invalid sorter '") + sorter + "'\n");
		throw std::invalid_argument("Invalid sorter on link-list");
	}
	Profiler::substages({ new ScalarProfile("n_cells", this),
	                      new EventProfile("icell", this),
	                      new EventProfile("ihoc", this),
	                      new EventProfile("link-list", this) });
}

LinkList::~LinkList()
{
	if (_min_pos)
		delete _min_pos;
	if (_max_pos)
		delete _max_pos;
	if (_sort)
		delete _sort;
	if (_ihoc)
		clReleaseKernel(_ihoc);
	if (_icell)
		clReleaseKernel(_icell);
	if (_ll)
		clReleaseKernel(_ll);
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
	if (_recompute_minmax) {
		_min_pos->setup();
		_max_pos->setup();
	}

	// Compute the cells length
	InputOutput::Variable* s = vars->get("support");
	InputOutput::Variable* h = vars->get("h");
	_cell_length = *(float*)s->get() * *(float*)h->get();

	// Setup the kernels
	setupOpenCL();

	// Setup the radix-sort
	_sort->setup();

	std::vector<std::string> inputs = {
		_input_name, "N", "support", "h"
	};
	std::vector<std::string> outputs = { _ihoc_name,
		                                 _icell_name,
		                                 _n_cells_name };
	if (_recompute_minmax) {
		outputs.push_back(_input_min_name);
		outputs.push_back(_input_max_name);
	} else {
		inputs.push_back(_input_min_name);
		inputs.push_back(_input_max_name);
	}
	setDependencies(inputs, outputs);

	// This tool shall mark ihoc variable as reallocatable
	InputOutput::ArrayVariable* ihoc =
	    (InputOutput::ArrayVariable*)getOutputDependencies()[0];
	ihoc->reallocatable(true);
}

template<>
void
LinkList::nCells<vec2>()
{
	auto vars = CalcServer::singleton()->variables();

	if (!_cell_length) {
		std::stringstream msg;
		msg << "Zero cell length detected in the tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid number of cells");
	}

	const vec2 pos_min = *(vec2*)vars->get("r_min")->get_async();
	const vec2 pos_max = *(vec2*)vars->get("r_max")->get_async();

	_n_cells.x = narrow_cast<ulcl>((pos_max.x - pos_min.x) / _cell_length) + 6;
	_n_cells.y = narrow_cast<ulcl>((pos_max.y - pos_min.y) / _cell_length) + 6;
	_n_cells.z = 1;
	const ulcl nw = _n_cells.x * _n_cells.y * _n_cells.z;
	_n_cells.w = nw;
}

template<>
void
LinkList::nCells<vec4>()
{
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	if (!_cell_length) {
		std::stringstream msg;
		msg << "Zero cell length detected in the tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid number of cells");
	}

	const vec4 pos_min = *(vec4*)vars->get("r_min")->get(true);
	const vec4 pos_max = *(vec4*)vars->get("r_max")->get(true);

	_n_cells.x = narrow_cast<ulcl>((pos_max.x - pos_min.x) / _cell_length) + 6;
	_n_cells.y = narrow_cast<ulcl>((pos_max.y - pos_min.y) / _cell_length) + 6;
	_n_cells.z = narrow_cast<ulcl>((pos_max.z - pos_min.z) / _cell_length) + 6;
	const ulcl nw = _n_cells.x * _n_cells.y * _n_cells.z;
	_n_cells.w = nw;
}

void
LinkList::allocate()
{
	cl_int err_code;
	auto C = CalcServer::singleton();

	_ihoc_gws = roundUp<size_t>(_n_cells.w, _ihoc_lws);

	// Check if the size of ihoc is already bigger than the required one
	auto ihoc_var = getOutputDependencies()[0];
	const size_t tsize = InputOutput::Variables::typeToBytes(ihoc_var->type());
	if (ihoc_var->size() / tsize >= _n_cells.w)
		return;

	// We have no alternative, we must sync here
	cl_mem mem = *(cl_mem*)ihoc_var->get();
	ihoc_var->sync();
	if (mem) {
		err_code = clReleaseMemObject(mem);
		CHECK_OCL_OR_THROW(err_code,
			std::string("Failure releasing ") +
			std::to_string(ihoc_var->size()) +
			" bytes on the device memory for tool \"" + name() + "\".");
	}
	mem = NULL;

	const size_t ihoc_t_size = ihoc_var->typesize();
	mem = clCreateBuffer(C->context(),
	                     CL_MEM_READ_WRITE,
	                     _n_cells.w * sizeof(ihoc_t_size),
	                     NULL,
	                     &err_code);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure allocating ") +
		std::to_string(_n_cells.w * sizeof(ihoc_t_size)) +
		" bytes on the device memory for tool \"" + name() + "\".");
	ihoc_var->set_async(&mem);
}

void
LinkList::setVariables()
{
	unsigned int i;
	cl_int err_code;
	InputOutput::Variables* vars = CalcServer::singleton()->variables();

	std::vector<std::string> _ihoc_vars = { _ihoc_name, "N", _n_cells_name };
	for (i = 0; i < _ihoc_vars.size(); i++) {
		InputOutput::Variable* var = vars->get(_ihoc_vars[i]);
		if (!memcmp(var->get_async(), _ihoc_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_ihoc, i, var->typesize(), var->get_async());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure setting the variable \"") +
		                       _ihoc_vars[i] + "\" to the tool \"" + name() +
		                       "\" (\"iHoc\").");
		memcpy(_ihoc_args.at(i), var->get_async(), var->typesize());
	}

	std::vector<std::string> _icell_vars = {
		_icell_name, _input_name, "N", _input_min_name,
		"support", "h", _n_cells_name
	};
	for (i = 0; i < _icell_vars.size(); i++) {
		InputOutput::Variable* var = vars->get(_icell_vars[i]);
		if (!memcmp(var->get_async(), _icell_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_icell, i, var->typesize(), var->get_async());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure setting the variable \"") +
		                       _ihoc_vars[i] + "\" to the tool \"" + name() +
		                       "\" (\"iCell\").");
		memcpy(_icell_args.at(i), var->get_async(), var->typesize());
	}

	std::vector<std::string> _ll_vars = { _icell_name, _ihoc_name, "N" };
	for (i = 0; i < _ll_vars.size(); i++) {
		InputOutput::Variable* var = vars->get(_ll_vars[i]);
		if (!memcmp(var->get_async(), _ll_args.at(i), var->typesize())) {
			continue;
		}
		err_code = clSetKernelArg(_ll, i, var->typesize(), var->get_async());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure setting the variable \"") +
		                       _ihoc_vars[i] + "\" to the tool \"" + name() +
		                       "\" (\"linkList\").");
		memcpy(_ll_args.at(i), var->get_async(), var->typesize());
	}
}

cl_event
LinkList::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	cl_event event;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	auto icell = getOutputDependencies()[1];
	auto n_cells = getOutputDependencies()[2];
	InputOutput::Variable *r_min, *r_max;

	if (_recompute_minmax) {
		// Reduction steps to find maximum and minimum position
		_min_pos->execute();
		_max_pos->execute();
		r_min = getOutputDependencies()[3];
		r_max = getOutputDependencies()[4];
	} else {
		r_min = getInputDependencies()[5];
		r_max = getInputDependencies()[6];
	}

	// Compute the number of cells and allocate ihoc accordingly
	const cl_event ncells_events[2] = { r_min->getWritingEvent(),
		                                r_max->getWritingEvent() };
	err_code = clWaitForEvents(2, ncells_events);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure waiting for the reductions on tool \"") + name() +
	        "\".");

	dynamic_cast<ScalarProfile*>(Profiler::substages()[0])->start();
	if (C->have_3d())
		nCells<vec4>();
	else
		nCells<vec2>();
	auto n_cells_var = getOutputDependencies()[2];
	if (!vars->isSameType(n_cells_var->type(), "svec4")) {
		std::stringstream msg;
		msg << "\"n_cells\" has and invalid type for \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"n_cells\" type is \"" << n_cells_var->type()
		    << "\", while \"" << vars->typeAlias("uivec4") << "\" was expected"
		    << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid n_cells type");
	}

	if (C->device_addr_bits() == 64)
		n_cells_var->set_async(&_n_cells);
	else {
		uivec4 n_cells = nCells32();
		n_cells_var->set_async(&n_cells);
	}

	if (_recompute_minmax) {
		std::stringstream report;
		report << name() << ":" << std::endl;
		report << r_min->name() << "=" << r_min->asString(true) << " ";
		report << r_max->name() << "=" << r_max->asString(true) << " ";
		report << n_cells->name() << "=" << n_cells->asString(true);
		InputOutput::Logger::singleton()->writeReport(report.str());
	}

	allocate();

	// Set the new allocated variables
	setVariables();
	dynamic_cast<ScalarProfile*>(Profiler::substages()[0])->end();

	// Compute the cell of each particle
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _icell,
	                                  1,
	                                  NULL,
	                                  &_icell_gws,
	                                  &_icell_lws,
	                                  events.size(),
	                                  events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing \"iCell\" from tool \"") +
	                       name() + "\".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
		name() + "\" when executing \"iCell\".");

	{
		auto profiler = dynamic_cast<EventProfile*>(Profiler::substages()[1]);
		profiler->start(event);
		profiler->end(event);
	}

	// Time to sort the icell array
	icell->setWritingEvent(event);
	n_cells->addReadingEvent(event);
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure releasing transactional \"iCell\" event from tool \"") +
	        name() + "\".");

	// Sort the particles from the cells
	_sort->execute();

	// Now our transactional event is the one coming from the sorting algorithm
	// Such a new event can be taken from the icell dependency
	// This transactional event SHALL NOT BE RELEASED. It is automagically
	// destroyed when no more variables use it
	auto event_wait = icell->getWritingEvent();

	// Reinit the head of cells
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _ihoc,
	                                  1,
	                                  NULL,
	                                  &_ihoc_gws,
	                                  &_ihoc_lws,
	                                  1,
	                                  &event_wait,
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing \"iHoc\" from tool \"") +
	                       name() + "\".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
		name() + "\" when executing \"iHoc\".");

	event_wait = event;
	{
		auto profiler = dynamic_cast<EventProfile*>(Profiler::substages()[2]);
		profiler->start(event);
		profiler->end(event);
	}

	// And compute them
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _ll,
	                                  1,
	                                  NULL,
	                                  &_ll_gws,
	                                  &_ll_lws,
	                                  1,
	                                  &event_wait,
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"linkList\" from tool \"") + name() +
	        "\".");
	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure releasing transactional \"linkList\" event in tool \"") +
	        name() + "\".");
	{
		auto profiler = dynamic_cast<EventProfile*>(Profiler::substages()[3]);
		profiler->start(event);
		profiler->end(event);
	}

	return event;
}

void
LinkList::setupOpenCL()
{
	unsigned int i;
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
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure querying the work group size (\"iHoc\") in tool \"") +
	        name() + "\".");
	if (_ihoc_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"iHoc\".\n");
		std::stringstream msg;
		msg << "\t" << _ihoc_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	size_t n_cells;
	if (C->device_addr_bits() == 64)
		n_cells = ((ulvec4*)vars->get("n_cells")->get())->w;
	else
		n_cells = ((uivec4*)vars->get("n_cells")->get())->w;
	_ihoc_gws = roundUp<size_t>(n_cells, _ihoc_lws);

	std::vector<std::string> _ihoc_vars = {
		_ihoc_name, "N", _n_cells_name
	};
	for (i = 0; i < _ihoc_vars.size(); i++) {
		err_code = clSetKernelArg(_ihoc,
		                          i,
		                          vars->get(_ihoc_vars[i])->typesize(),
		                          vars->get(_ihoc_vars[i])->get());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure sending \"") + _ihoc_vars[i] +
		                       "\" argument to \"iHoc\" in tool \"" + name() +
		                       "\".");
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
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure querying the work group size (\"iCell\") in tool \"") +
	        name() + "\".");
	if (_icell_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"iCell\".\n");
		std::stringstream msg;
		msg << "\t" << _icell_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	size_t N;
	if (C->device_addr_bits() == 64)
		N = *(ulcl*)vars->get("N")->get();
	else
		N = *(uicl*)vars->get("N")->get();
	_icell_gws = roundUp<size_t>(N, _icell_lws);

	std::vector<std::string> _icell_vars = {
		_icell_name, _input_name, "N",
		_input_min_name, "support", "h", _n_cells_name
	};
	for (i = 0; i < _icell_vars.size(); i++) {
		err_code = clSetKernelArg(_icell,
		                          i,
		                          vars->get(_icell_vars[i])->typesize(),
		                          vars->get(_icell_vars[i])->get());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure sending \"") + _icell_vars[i] +
		                       "\" argument to \"iCell\" in tool \"" + name() +
		                       "\".");
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
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure querying the work group size (\"linkList\") in tool \"") +
	        name() + "\".");
	if (_ll_lws < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory for \"linkList\".\n");
		std::stringstream msg;
		msg << "\t" << _ll_lws
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	_ll_gws = roundUp<size_t>(N, _ll_lws);

	std::vector<std::string> _ll_vars = {
		_icell_name, _ihoc_name, "N"
	};
	for (i = 0; i < _ll_vars.size(); i++) {
		err_code = clSetKernelArg(_ll,
		                          i,
		                          vars->get(_ll_vars[i])->typesize(),
		                          vars->get(_ll_vars[i])->get());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure sending \"") + _ll_vars[i] +
		                       "\" argument to \"iCell\" in tool \"" + name() +
		                       "\".");
		_ll_args.push_back(malloc(vars->get(_ll_vars[i])->typesize()));
		memcpy(_ll_args.at(i),
		       vars->get(_ll_vars[i])->get(),
		       vars->get(_ll_vars[i])->typesize());
	}
}

}
} // namespace
