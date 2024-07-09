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
 * @brief Methods to perform a radix sort using the GPU (or any device
 * supported by OpenCL).
 * (See Aqua::CalcServer::RadixSort for details)
 * @note Hardcoded versions of the files CalcServer/RadixSort.cl.in and
 * CalcServer/RadixSort.hcl.in are internally included as a text array.
 */

#include <sstream>
#include <limits>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "RadixSort.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/RadixSort.hcl"
#include "aquagpusph/CalcServer/RadixSort.cl"
#endif
std::string RADIXSORT_INC = xxd2string(RadixSort_hcl_in, RadixSort_hcl_in_len);
std::string RADIXSORT_SRC = xxd2string(RadixSort_cl_in, RadixSort_cl_in_len);

RadixSort::RadixSort(const std::string tool_name,
                     const std::string variable,
                     const std::string permutations,
                     const std::string inv_permutations,
                     bool once)
  : Tool(tool_name, once)
  , _var_name(variable)
  , _perms_name(permutations)
  , _inv_perms_name(inv_permutations)
  , _var(NULL)
  , _perms(NULL)
  , _inv_perms(NULL)
  , _n(0)
  , _n_padded(0)
  , _init_kernel(NULL)
  , _histograms_kernel(NULL)
  , _scan_kernel(NULL)
  , _paste_kernel(NULL)
  , _sort_kernel(NULL)
  , _inv_perms_kernel(NULL)
  , _in_vals(NULL)
  , _out_vals(NULL)
  , _in_permut(NULL)
  , _out_permut(NULL)
  , _histograms(NULL)
  , _global_sums(NULL)
  , _temp_mem(NULL)
  , _items(_ITEMS)
  , _groups(_GROUPS)
  , _bits(_STEPBITS)
  , _radix(_RADIX)
  , _histo_split(_HISTOSPLIT)
{
	Profiler::substages({ new EventProfile("init", this),
	                      new EventProfile("radix-sort", this),
	                      new EventProfile("keys", this),
	                      new EventProfile("values", this),
	                      new EventProfile("inverse keys", this) });
}

RadixSort::~RadixSort()
{
	if (_init_kernel)
		clReleaseKernel(_init_kernel);
	if (_histograms_kernel)
		clReleaseKernel(_histograms_kernel);
	if (_scan_kernel)
		clReleaseKernel(_scan_kernel);
	if (_paste_kernel)
		clReleaseKernel(_paste_kernel);
	if (_sort_kernel)
		clReleaseKernel(_sort_kernel);
	if (_inv_perms_kernel)
		clReleaseKernel(_inv_perms_kernel);
	if (_in_vals)
		clReleaseMemObject(_in_vals);
	if (_out_vals)
		clReleaseMemObject(_out_vals);
	if (_in_permut)
		clReleaseMemObject(_in_permut);
	if (_out_permut)
		clReleaseMemObject(_out_permut);
	if (_histograms)
		clReleaseMemObject(_histograms);
	if (_global_sums)
		clReleaseMemObject(_global_sums);
	if (_temp_mem)
		clReleaseMemObject(_temp_mem);
}

void
RadixSort::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();

	// Get the variables
	variables();

	// Setup the working tools
	setupOpenCL();
}

cl_event
RadixSort::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
	cl_int err_code;
	size_t max_val;
	cl_event event;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Fisrt we copy everything on our transactional memory objects, which are
	// conveniently padded
	// We also initialize the permutations as the particles id, i.e. each
	// particle is converted on itself.
	auto event_init = init();
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
		name() + "\" at the initialization.");
	auto init_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(0));
	init_profiler->start(event_init);
	init_profiler->end(event_init);

	// Get maximum key bits and the associated needed passes
	if (C->device_addr_bits() == 64)
		max_val = std::numeric_limits<ulcl>::max() >> 1;
	else
		max_val = std::numeric_limits<uicl>::max() >> 1;
	if (_var_name == "icell") {
		size_t n_cells;
		if (C->device_addr_bits() == 64)
			n_cells = ((ulvec4*)vars->get("n_cells")->get())->w;
		else
			n_cells = ((uivec4*)vars->get("n_cells")->get())->w;
		max_val = n_cells;
	}
	max_val = nextPowerOf2(max_val);
	for (_key_bits = 1; (max_val & 1) == 0; max_val >>= 1, _key_bits++)
		;
	_key_bits = roundUp<size_t>(_key_bits, _bits);
	if (_key_bits > C->device_addr_bits()) {
		LOG(L_ERROR,
		    std::string("Resultant keys=") + std::to_string(_key_bits) +
		    "bits overflows the " + std::to_string(C->device_addr_bits()) +
		    "bits available for unsigned int types\n");
		throw std::runtime_error("Unsigned int overflow");
	}
	_n_pass = _key_bits / _bits;
	size_t uint_size = C->device_addr_bits() / 8;

	// Time to sort everything up
	auto radixsort_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(1));
	cl_event event_wait = NULL;
	for (_pass = 0; _pass < _n_pass; _pass++) {
		event_wait = histograms(event_init, event_wait);
		if (_pass == 0)
			radixsort_profiler->start(event_wait);
		event_wait = scan(event_wait);
		event_wait = reorder(event_init, event_wait);
		err_code = clFlush(C->command_queue());
		CHECK_OCL_OR_THROW(err_code,
			std::string("Failure flushing the command queue in tool \"") +
			name() + "\" at pass " + std::to_string(_pass) + ".");
	}

	err_code = clReleaseEvent(event_init);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing permutations event in tool \"") +
	        name() + "\".");
	radixsort_profiler->end(event_wait);

	// Copy back the results
	auto var_events = _var->getReadingEvents();
	// We positively know that the last reading event is posterior to the
	// last writing event, because of the clEnqueueCopyBuffer() above
	// We can anyway wait for it
	var_events.push_back(_var->getWritingEvent());
	var_events.push_back(event_wait);
	err_code = clEnqueueCopyBuffer(C->command_queue(),
	                               _in_vals,
	                               *(cl_mem*)_var->get(),
	                               0,
	                               0,
	                               _n * uint_size,
	                               var_events.size(),
	                               var_events.data(),
	                               &event);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure copying the sort keys in tool \"") + name() +
			"\".");
	_var->setWritingEvent(event);
	auto vals_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(3));
	vals_profiler->start(event);
	vals_profiler->end(event);
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing copying event in tool \"") + name() +
	        "\".");

	auto perms_events = _perms->getReadingEvents();
	perms_events.push_back(_perms->getWritingEvent());
	perms_events.push_back(event_wait);
	err_code = clEnqueueCopyBuffer(C->command_queue(),
	                               _in_permut,
	                               *(cl_mem*)_perms->get(),
	                               0,
	                               0,
	                               _n * uint_size,
	                               perms_events.size(),
	                               perms_events.data(),
	                               &event);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure copying the permutations in tool \"") + name() +
			"\".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
		name() + "\" when copying.");

	_perms->setWritingEvent(event);
	auto keys_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(2));
	keys_profiler->start(event);
	keys_profiler->end(event);

	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing permutations event in tool \"") +
	        name() + "\".");

	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing transactional event in tool \"") +
	        name() + "\".");

	event_wait = inversePermutations();
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
		name() + "\" when computing the inverse permutations.");
	auto inv_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(4));
	inv_profiler->start(event_wait);

	// The events associated to _var and _perms was already set during the
	// execution of this function. However, if we return the event coming
	// from inversePermutations(), their events will get overwritten by
	// Aqua::CalcServer::Tool::execute()
	// Thus we can join all the events together with a marker
	try {
		event = C->marker(C->command_queue(), { _var->getWritingEvent(),
		                                        _perms->getWritingEvent(),
		                                        event_wait });
	} catch (std::runtime_error& e) {
		LOG(L_ERROR, std::string("While joining the output events in tool ") +
		             name() + ".\n");
	}
	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure releasing the inverse permutations event in tool \"") +
	        name() + "\".");
	inv_profiler->end(event);

	return event;
}

cl_event
RadixSort::init()
{
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();

	err_code =
	    clSetKernelArg(_init_kernel, 1, sizeof(cl_mem), (void*)&_in_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"init\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_init_kernel, 2, sizeof(cl_mem), (void*)&_in_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"init\" in tool \"") +
	        name() + "\".");

	auto var_event = _var->getWritingEvent();
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _init_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  NULL,
	                                  1,
	                                  &var_event,
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"init\" within in tool \"") + name() +
	        "\".");

	_var->addReadingEvent(event);

	return event;
}

cl_event
RadixSort::histograms(cl_event keys_event, cl_event histograms_event)
{
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();
	size_t local_work_size = _items;
	size_t global_work_size = _groups * _items;

	err_code =
	    clSetKernelArg(_histograms_kernel, 0, sizeof(cl_mem), (void*)&_in_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"histogram\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_histograms_kernel, 2, sizeof(uicl), (void*)&_pass);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"histogram\" in tool \"") +
	        name() + "\".");

	std::vector<cl_event> events = { keys_event };
	if (histograms_event)
		events.push_back(histograms_event);

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _histograms_kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &local_work_size,
	                                  events.size(),
	                                  events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"histogram\" within in tool \"") +
	        name() + "\".");

	if (histograms_event) {
		err_code = clReleaseEvent(histograms_event);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure releasing \"histogram\" event in tool \"") +
		        name() + "\".");
	}

	return event;
}

cl_event
RadixSort::scan(cl_event event_histo)
{
	cl_int err_code;
	cl_event event, event_wait = event_histo;
	CalcServer* C = CalcServer::singleton();
	size_t global_work_size = _radix * _groups * _items / 2;
	size_t local_work_size = global_work_size / _histo_split;

	// 1st scan
	// ========
	err_code =
	    clSetKernelArg(_scan_kernel, 0, sizeof(cl_mem), (void*)&_histograms);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"scan\" (1st call) ") +
	        "in tool \"" + name() + "\".");
	err_code =
	    clSetKernelArg(_scan_kernel, 2, sizeof(cl_mem), (void*)&_global_sums);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"scan\" (1st call) ") +
	        "in tool \"" + name() + "\".");

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _scan_kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &local_work_size,
	                                  1,
	                                  &event_wait,
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"scan\" (1st call) in tool \"") +
	        name() + "\".");
	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing \"scan\" (1st call) event in tool \"") +
	        name() + "\".");
	event_wait = event;

	// 2nd scan
	// ========
	global_work_size = _histo_split / 2;
	local_work_size = global_work_size;
	err_code =
	    clSetKernelArg(_scan_kernel, 0, sizeof(cl_mem), (void*)&_global_sums);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"scan\" (2nd call) ") +
	        "in tool \"" + name() + "\".");
	err_code =
	    clSetKernelArg(_scan_kernel, 2, sizeof(cl_mem), (void*)&_temp_mem);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"scan\" (2nd call) ") +
	        "in tool \"" + name() + "\".");

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _scan_kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &local_work_size,
	                                  1,
	                                  &event_wait,
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"scan\" (2nd call) in tool \"") +
	        name() + "\".");
	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing \"scan\" (2nd call) event in tool \"") +
	        name() + "\".");
	event_wait = event;

	// Histograms paste
	// ================
	global_work_size = _radix * _groups * _items / 2;
	local_work_size = global_work_size / _histo_split;

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _paste_kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &local_work_size,
	                                  1,
	                                  &event_wait,
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing \"paste\" in tool \"") +
	                       name() + "\".");
	err_code = clReleaseEvent(event_wait);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing \"paste\" event in tool \"") + name() +
	        "\".");

	return event;
}

cl_event
RadixSort::reorder(cl_event perms_event, cl_event histograms_event)
{
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();
	size_t local_work_size = _items;
	size_t global_work_size = _groups * _items;

	err_code =
	    clSetKernelArg(_sort_kernel, 0, sizeof(cl_mem), (void*)&_in_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_sort_kernel, 1, sizeof(cl_mem), (void*)&_out_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(_sort_kernel, 3, sizeof(uicl), (void*)&_pass);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 3 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_sort_kernel, 4, sizeof(cl_mem), (void*)&_in_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 4 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_sort_kernel, 5, sizeof(cl_mem), (void*)&_out_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 5 to \"sort\" in tool \"") +
	        name() + "\".");

	std::vector<cl_event> events = { perms_event, histograms_event };
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _sort_kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &local_work_size,
	                                  events.size(),
	                                  events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing \"sort\" in tool \"") +
	                       name() + "\".");
	err_code = clReleaseEvent(histograms_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing \"sort\" event in tool \"") + name() +
	        "\".");

	// Swap the memory identifiers for the next pass
	cl_mem d_temp;

	d_temp = _in_vals;
	_in_vals = _out_vals;
	_out_vals = d_temp;

	d_temp = _in_permut;
	_in_permut = _out_permut;
	_out_permut = d_temp;

	return event;
}

cl_event
RadixSort::inversePermutations()
{
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();

	err_code = clSetKernelArg(
	    _inv_perms_kernel, 0, sizeof(cl_mem), (void*)&_in_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"inversePermutation\" ") +
	        "in tool \"" + name() + "\".");
	err_code = clSetKernelArg(
	    _inv_perms_kernel, 1, sizeof(cl_mem), _inv_perms->get_async());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"inversePermutation\" ") +
	        "in tool \"" + name() + "\".");

	auto perms_events = _inv_perms->getReadingEvents();
	perms_events.push_back(_inv_perms->getWritingEvent());
	// We just wrote on _perms before calling this function, so no need to add
	// its reading events
	perms_events.push_back(_perms->getWritingEvent());
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _inv_perms_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  NULL,
	                                  perms_events.size(),
	                                  perms_events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"inversePermutation\" in tool \"") +
	        name() + "\".");

	return event;
}

void
RadixSort::variables()
{
	size_t n;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	std::string size_t_name = (C->device_addr_bits() == 64) ?
		"unsigned long*" :
		"unsigned int*";
	// Check and get the variables
	if (!vars->get(_var_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared variable \"" << _var_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_var_name)->type() != size_t_name) {
		std::ostringstream msg;
		msg << "Tool \"" << name() << "\" cannot process variable \""
		    << _var_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"" << size_t_name << "\" type was expected, but \""
		    << vars->get(_var_name)->type() << "\" has been received."
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_var = (InputOutput::ArrayVariable*)vars->get(_var_name);

	if (!vars->get(_perms_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared permutations variable \""
		    << _perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_perms_name)->type() != size_t_name) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" cannot process permutations variable \"" << _perms_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"" << size_t_name << "\" type was expected, but \""
		    << vars->get(_var_name)->type() << "\" has been received."
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_perms = (InputOutput::ArrayVariable*)vars->get(_perms_name);

	if (!vars->get(_inv_perms_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared inverse permutations variable "
		       "\""
		    << _inv_perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_inv_perms_name)->type() != size_t_name) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" cannot process inverse permutations variable \""
		    << _inv_perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"" << size_t_name << "\" type was expected, but \""
		    << vars->get(_var_name)->type() << "\" has been received."
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_inv_perms = (InputOutput::ArrayVariable*)vars->get(_inv_perms_name);

	// Check the lengths
	_n = _var->size() / vars->typeToBytes(_var->type());
	_n_padded = nextPowerOf2(roundUp<size_t>(_n, _ITEMS * _GROUPS));

	n = _perms->size() / vars->typeToBytes(_perms->type());
	if (n != _n) {
		std::ostringstream msg;
		msg << "Lengths mismatch in tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _var->name() << "\" has length, n=" << _n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _perms->name() << "\" has length, n=" << n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable length");
	}
	n = _inv_perms->size() / vars->typeToBytes(_inv_perms->type());
	if (n != _n) {
		std::ostringstream msg;
		msg << "Lengths mismatch in tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _var->name() << "\" has length, n=" << _n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _inv_perms->name() << "\" has length, n=" << n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable length");
	}

	setOutputDependencies({ _var, _perms, _inv_perms });
}

void
RadixSort::setupOpenCL()
{
	std::ostringstream source;
	source << RADIXSORT_INC << RADIXSORT_SRC;
	std::ostringstream flags;
	flags << " -D_BITS=" << _bits << " -D_RADIX=" << _radix;
	std::vector<cl_kernel> kernels = compile(
	    source.str(),
	    { "init", "histogram", "scan", "paste", "sort", "inversePermutation" },
	    flags.str());
	_init_kernel = kernels.at(0);
	_histograms_kernel = kernels.at(1);
	_scan_kernel = kernels.at(2);
	_paste_kernel = kernels.at(3);
	_sort_kernel = kernels.at(4);
	_inv_perms_kernel = kernels.at(5);

	// Check and correct _items, _groups and _histo_split
	setupDims();

	// Setup the memory objects
	setupMems();
	setupArgs();

	std::ostringstream msg;
	msg << "\titems: " << _items << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tgroups: " << _groups << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tsplits: " << _histo_split << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tn: " << _n << " -> " << _n_padded << std::endl;
	LOG0(L_DEBUG, msg.str());
}

void
RadixSort::setupDims()
{
	cl_int err_code;
	size_t max_local_work_size, sort_local_work_size, scan_local_work_size;
	CalcServer* C = CalcServer::singleton();

	// For the _histograms_kernel and _sort_kernel _items can be used as the
	// upper bound
	err_code = clGetKernelWorkGroupInfo(_histograms_kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &max_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure getting CL_KERNEL_WORK_GROUP_SIZE from ") +
	        "\"histogram\" in tool " + name() + "\".");
	err_code = clGetKernelWorkGroupInfo(_sort_kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &sort_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure getting CL_KERNEL_WORK_GROUP_SIZE from ") +
	        "\"sort\" in tool " + name() + "\".");
	if (sort_local_work_size < max_local_work_size)
		max_local_work_size = sort_local_work_size;
	if (max_local_work_size < _items)
		_items = max_local_work_size;
	if (!isPowerOf2(_items))
		_items = nextPowerOf2(_items) / 2;

	// The _scan_kernel can be used to set an upper bound to the number of
	// histogram splits
	err_code = clGetKernelWorkGroupInfo(_scan_kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &scan_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure getting CL_KERNEL_WORK_GROUP_SIZE from ") +
	        "\"scan\" in tool " + name() + "\".");
	if (scan_local_work_size < _histo_split / 2)
		_histo_split = 2 * scan_local_work_size;
	if (!isPowerOf2(_histo_split))
		_histo_split = nextPowerOf2(_histo_split) / 2;

	// Finally using _scan_kernel and _paste_kernel we can adjust _groups,
	// _items and _radix
	err_code = clGetKernelWorkGroupInfo(_paste_kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &max_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure getting CL_KERNEL_WORK_GROUP_SIZE from ") +
	        "\"paste\" in tool " + name() + "\".");
	max_local_work_size = min(max_local_work_size, scan_local_work_size);
	while (max_local_work_size < _radix * _groups * _items / 2 / _histo_split) {
		// We can't increase _histo_split, so we may start decreasing the number
		// of items
		_items /= 2;
		if (_items < __CL_MIN_LOCALSIZE__) {
			_items = __CL_MIN_LOCALSIZE__;
			break;
		}
	}
	while (max_local_work_size < _radix * _groups * _items / 2 / _histo_split) {
		// We have reached the minimum possible number of items, so we can
		// continue decreasing the number of groups
		_groups /= 2;
		if (!_groups) {
			_groups = 1;
			break;
		}
	}
	if (max_local_work_size < _radix * _groups * _items / 2 / _histo_split) {
		// We can try to reduce the radix, but it is a bad business
		LOG(L_ERROR,
		    "Failure setting a number of items and groups compatible with "
		    "\"scan\" and \"paste\".\n");
		LOG0(L_DEBUG,
		     "\tYou can try to recompile the code decreasing "
		     "__CL_MIN_LOCALSIZE__\n");
		throw std::runtime_error("OpenCL error");
	}

	_local_work_size = getLocalWorkSize(C->command_queue());
	_global_work_size = getGlobalWorkSize(_n_padded, _local_work_size);
}

void
RadixSort::setupMems()
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	if (_in_vals)
		clReleaseMemObject(_in_vals);
	_in_vals = NULL;
	if (_out_vals)
		clReleaseMemObject(_out_vals);
	_out_vals = NULL;
	if (_in_permut)
		clReleaseMemObject(_in_permut);
	_in_permut = NULL;
	if (_out_permut)
		clReleaseMemObject(_out_permut);
	_out_permut = NULL;
	if (_histograms)
		clReleaseMemObject(_histograms);
	_histograms = NULL;
	if (_global_sums)
		clReleaseMemObject(_global_sums);
	_global_sums = NULL;
	if (_temp_mem)
		clReleaseMemObject(_temp_mem);
	_temp_mem = NULL;
	allocatedMemory(0);

	// Get the memory identifiers
	size_t type_size = (C->device_addr_bits() == 64) ?
		sizeof(ulcl) : sizeof(uicl);

	_in_vals = clCreateBuffer(C->context(),
	                          CL_MEM_READ_WRITE,
	                          _n_padded * type_size,
	                          NULL,
	                          &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_out_vals = clCreateBuffer(C->context(),
	                           CL_MEM_READ_WRITE,
	                           _n_padded * type_size,
	                           NULL,
	                           &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_in_permut = clCreateBuffer(C->context(),
	                            CL_MEM_READ_WRITE,
	                            _n_padded * type_size,
	                            NULL,
	                            &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_out_permut = clCreateBuffer(C->context(),
	                             CL_MEM_READ_WRITE,
	                             _n_padded * type_size,
	                             NULL,
	                             &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_histograms =
	    clCreateBuffer(C->context(),
	                   CL_MEM_READ_WRITE,
	                   (_radix * _groups * _items) * type_size,
	                   NULL,
	                   &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_global_sums = clCreateBuffer(C->context(),
	                              CL_MEM_READ_WRITE,
	                              _histo_split * type_size,
	                              NULL,
	                              &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");
	_temp_mem = clCreateBuffer(
	    C->context(), CL_MEM_READ_WRITE, sizeof(unsigned int), NULL, &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating device memory in tool \"") + name() +
	        "\".");

	allocatedMemory((4 * _n_padded + (_radix * _groups * _items) +
	                 _histo_split + 1) * type_size);
}

template<typename T>
void
RadixSort::setupTypedArgs()
{
	cl_int err_code;

	T n = narrow_cast<T>(_n);
	T n_padded = narrow_cast<T>(_n_padded);

	err_code =
	    clSetKernelArg(_init_kernel, 0, sizeof(cl_mem), _var->get_async());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"init\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(_init_kernel, 3, sizeof(T), (void*)&n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 3 to \"init\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_init_kernel, 4, sizeof(T), (void*)&n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 4 to \"init\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(
	    _histograms_kernel, 1, sizeof(cl_mem), (void*)&_histograms);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"histogram\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(
	    _histograms_kernel, 3, sizeof(T) * _radix * _items, NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 3 to \"histogram\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(
	    _histograms_kernel, 4, sizeof(T), (void*)&n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 4 to \"histogram\" in tool \"") +
	        name() + "\".");

	unsigned int maxmemcache =
	    max(_histo_split, _items * _groups * _radix / _histo_split);
	err_code =
	    clSetKernelArg(_scan_kernel, 1, sizeof(T) * maxmemcache, NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"scan\" in tool \"") +
	        name() + "\".");

	err_code =
	    clSetKernelArg(_paste_kernel, 0, sizeof(cl_mem), (void*)&_histograms);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"paste\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_paste_kernel, 1, sizeof(cl_mem), (void*)&_global_sums);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"paste\" in tool \"") +
	        name() + "\".");

	err_code =
	    clSetKernelArg(_sort_kernel, 2, sizeof(cl_mem), (void*)&_histograms);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(
	    _sort_kernel, 6, sizeof(T) * _radix * _items, NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 6 to \"sort\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_sort_kernel, 7, sizeof(T), (void*)&n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 7 to \"sort\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(
		_inv_perms_kernel, 1, sizeof(cl_mem), _inv_perms->get());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"inversePermutation\" ") +
	        "in tool \"" + name() + "\".");
	err_code =
	    clSetKernelArg(_inv_perms_kernel, 2, sizeof(T), (void*)&n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"inversePermutation\" ") +
	        "in tool \"" + name() + "\".");
}

void
RadixSort::setupArgs()
{
	if (CalcServer::singleton()->device_addr_bits() == 64)
		setupTypedArgs<ulcl>();
	else
		setupTypedArgs<uicl>();
}

}
} // namespace
