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
 * @brief Synchronize arrays between processes, sending information by
 * the network.
 * (See Aqua::CalcServer::MPISync for details)
 * @note Hardcoded versions of the files CalcServer/MPISync.cl.in and
 * CalcServer/MPISync.hcl.in are internally included as a text array.
 */

#include "aquagpusph/sphPrerequisites.hpp"

#ifdef HAVE_MPI

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "MPISync.hpp"
#include "Kernel.hpp"
#include "LinkList.hpp"
#include "CalcServer.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/MPISync.hcl"
#include "aquagpusph/CalcServer/MPISync.cl"
#endif
std::string MPISYNC_INC = xxd2string(MPISync_hcl_in, MPISync_hcl_in_len);
std::string MPISYNC_SRC = xxd2string(MPISync_cl_in, MPISync_cl_in_len);

#define cmd_queue_mpi CalcServer::cmd_queue::cmd_queue_mpi

MPISync::MPISync(const std::string name,
                 const std::string mask,
                 const std::vector<std::string> fields,
                 const std::vector<unsigned int> procs,
                 bool once)
  : Tool(name, once)
  , _mask_name(mask)
  , _mask(NULL)
  , _field_names(fields)
  , _procs(procs)
  , _unsorted_id(NULL)
  , _sorted_id(NULL)
  , _sort(NULL)
  , _n(0)
  , _n_offset_recv(NULL)
{
	int mpi_rank, mpi_size;
	mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
	mpi_size = Aqua::MPI::size(MPI_COMM_WORLD);
	assert(mpi_rank >= 0);
	assert(mpi_size > 0);

	// Check that we have a good list of processes
	if (_procs.size() == 0) {
		for (unsigned int proc = 0; proc < mpi_size; proc++) {
			_procs.push_back(proc);
		}
	}
	for (unsigned int i = 0; i < _procs.size(); i++) {
		if ((_procs.at(i) == mpi_rank) || (_procs.at(i) >= mpi_size)) {
			_procs.erase(_procs.begin() + i);
		}
	}
}

MPISync::~MPISync()
{
	if (_sort)
		delete _sort;
	_sort = NULL;
	if (_n_offset_recv_reinit)
		delete _n_offset_recv_reinit;
	_n_offset_recv_reinit = NULL;
	if (_mask_reinit)
		delete _mask_reinit;
	_mask_reinit = NULL;

	// We must unset the inner memory objects of the sorted fields to avoid
	// releasing them twice while deleting the sorters
	cl_mem inner_mem = NULL;
	for (auto field : _fields_sorted) {
		field->reallocatable(true);
		field->set((void*)(&inner_mem));
		field->reallocatable(false);
	}
	// Now we can delete the sorters
	for (auto sorter : _field_sorters) {
		delete sorter;
	}
	// And the senders
	for (auto field : _fields_send) {
		free(field);
	}
	for (auto sender : _senders) {
		delete sender;
	}
	// And the receivers
	for (auto field : _fields_recv) {
		free(field);
	}
	for (auto receiver : _receivers) {
		delete receiver;
	}
}

void
MPISync::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	std::vector<Profile*> profilers;

	// Get the involved variables
	variables();
	// Setup the mask sorting subtool
	try {
		setupSort();
	} catch (std::runtime_error& e) {
		std::stringstream msg;
		msg << "Error setting up the sorter for tool \"" << name() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
	// Setup the field sorting subtools
	for (auto field : _fields) {
		setupFieldSort(field);
	}

	std::vector<InputOutput::Variable*> deps;
	for (auto field : _fields) {
		deps.push_back((InputOutput::Variable*)field);
	}
	deps.push_back((InputOutput::Variable*)_mask);
	for (auto field : _fields_sorted) {
		deps.push_back((InputOutput::Variable*)field);
	}
	deps.push_back((InputOutput::Variable*)_unsorted_id);
	deps.push_back((InputOutput::Variable*)_sorted_id);
	setOutputDependencies(deps);

	setupSenders();
	for (auto s : _senders) {
		std::stringstream pname;
		pname << "sender " << s->proc();
		profilers.push_back(new EventProfile(pname.str(), this));
	}

	setupReceivers();
	for (auto s : _senders) {
		std::stringstream pname;
		pname << "receiver " << s->proc();
		profilers.push_back(new EventProfile(pname.str(), this));
	}

	Profiler::substages(profilers);
}

cl_event
MPISync::_execute(const std::vector<cl_event> events)
{
	unsigned int i;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	if (!_procs.size())
		return NULL;

	// Sort the mask
	try {
		_sort->execute();
	} catch (std::runtime_error& e) {
		std::stringstream msg;
		msg << "Error while sorting the mask for tool \"" << name() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
	// Sort the fields
	for (auto sorter : _field_sorters) {
		try {
			sorter->execute();
		} catch (std::runtime_error& e) {
			std::stringstream msg;
			msg << "Error while sorting \"" << sorter->input()->name()
			    << "\" for tool \"" << name() << "\"" << std::endl;
			LOG(L_ERROR, msg.str());
			throw;
		}
	}

	// Send the data to other processes
	unsigned int i_profiler = 0;
	for (auto sender : _senders) {
		auto profiler =
		    dynamic_cast<EventProfile*>(Profiler::substages().at(i_profiler));
		sender->execute(profiler);
		i_profiler++;
	}

	// Receive the data from the other processes
	_n_offset_recv_reinit->execute();
	_mask_reinit->execute();
	for (auto receiver : _receivers) {
		auto profiler =
		    dynamic_cast<EventProfile*>(Profiler::substages().at(i_profiler));
		receiver->execute(profiler);
		i_profiler++;
	}

	return NULL;
}

void
MPISync::variables()
{
	size_t n;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	if (!vars->get(_mask_name)) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is asking the undeclared variable \"" << _mask_name << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_mask_name)->type().compare("unsigned int*")) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is asking the variable \""
		    << _mask_name << "\", which has an invalid type" << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"unsigned int*\" was expected, but \""
		    << vars->get(_mask_name)->type() << "\" was found." << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_mask = (InputOutput::ArrayVariable*)vars->get(_mask_name);
	_n = _mask->size() / InputOutput::Variables::typeToBytes(_mask->type());

	for (auto var_name : _field_names) {
		if (!vars->get(var_name)) {
			std::stringstream msg;
			msg << "The tool \"" << name()
			    << "\" is asking the undeclared variable \"" << var_name
			    << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable");
		}
		if (vars->get(var_name)->type().find('*') == std::string::npos) {
			std::stringstream msg;
			msg << "The tool \"" << name()
			    << "\" may not use a scalar variable (\"" << var_name << "\")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
		InputOutput::ArrayVariable* var =
		    (InputOutput::ArrayVariable*)vars->get(var_name);
		n = var->size() / InputOutput::Variables::typeToBytes(var->type());
		if (n != _n) {
			std::stringstream msg;
			msg << "Wrong variable length in the tool \"" << name() << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			msg.str("");
			msg << "\t\"" << _mask_name << "\" has length " << _n << std::endl;
			LOG0(L_DEBUG, msg.str());
			msg.str("");
			msg << "\t\"" << var_name << "\" has length " << n << std::endl;
			LOG0(L_DEBUG, msg.str());
			throw std::runtime_error("Invalid variable length");
		}
		_fields.push_back((InputOutput::ArrayVariable*)vars->get(var_name));
	}
}

void
MPISync::setupSort()
{
	std::string var_name;
	std::ostringstream valstr;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Register the variables to store the permutations
	valstr << _n;
	var_name = varPrefix() + _mask_name + "_unsorted";
	vars->registerVariable(var_name, "unsigned int*", valstr.str(), "");
	_unsorted_id = (InputOutput::ArrayVariable*)vars->get(var_name);
	var_name = varPrefix() + _mask_name + "_sorted";
	vars->registerVariable(var_name, "unsigned int*", valstr.str(), "");
	_sorted_id = (InputOutput::ArrayVariable*)vars->get(var_name);

	// Create the sorter
	var_name = varPrefix() + _mask_name + "->Radix-Sort";
	_sort = new RadixSort(
	    var_name, _mask_name, _unsorted_id->name(), _sorted_id->name());
	_sort->setup();
}

void
MPISync::setupFieldSort(InputOutput::ArrayVariable* field)
{
	std::string var_name;
	std::ostringstream valstr;
	cl_int err_code;
	cl_mem inner_mem;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Register the variable to store the sorted copy
	valstr.str("");
	valstr << _n;
	var_name = varPrefix() + field->name() + "_sorted";
	vars->registerVariable(var_name, field->type(), valstr.str(), "");
	_fields_sorted.push_back((InputOutput::ArrayVariable*)vars->get(var_name));
	// Remove the inner memory object, since we are using the one computed by
	// the sorting tool
	_fields_sorted.back()->reallocatable(true);
	inner_mem = *(cl_mem*)_fields_sorted.back()->get();
	err_code = clReleaseMemObject(inner_mem);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure Releasing the inner memory object of \"") +
	        _fields_sorted.back()->name() + "\" in tool \"" + name() + "\".");

	// Create the sorter
	var_name = "__" + field->name() + "->Radix-Sort";
	_field_sorters.push_back(
	    new UnSort(var_name, field->name(), _sorted_id->name()));
	_field_sorters.back()->setup();

	// Replace the inner CL memory object of the sorted field by the one
	// computed by the sorter
	inner_mem = _field_sorters.back()->output();
	_fields_sorted.back()->set((void*)(&inner_mem));
	_fields_sorted.back()->reallocatable(false);
}

void
MPISync::setupSenders()
{
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Allocate the host memory to download the fields and subsequently send
	// them
	for (auto field : _fields) {
		void* data = malloc(field->size());
		if (!data) {
			std::stringstream msg;
			msg << "Failure allocating \"" << field->size()
			    << "\" bytes for the array \"" << field->name()
			    << "\" send in the tool \"" << name() << "\"" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		_fields_send.push_back(data);
	}

	// Create the senders
	for (auto proc : _procs) {
		Sender* sender = new Sender(name(),
		                            varPrefix(),
		                            _mask,
		                            _fields_sorted,
		                            _fields_send,
		                            proc);
		if (!sender) {
			std::stringstream msg;
			msg << "Failure Allocating memory for the process " << proc
			    << " sender in tool \"" << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		_senders.push_back(sender);
	}
}

void
MPISync::setupReceivers()
{
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Create a variable for the cumulative offset computation
	std::string var_name = varPrefix() + "mpi_offset";
	vars->registerVariable(var_name, "unsigned int", "", "0");
	_n_offset_recv = (InputOutput::UIntVariable*)vars->get(var_name);

	// Create a tool to reinit it to zero value
	_n_offset_recv_reinit = new SetScalar(
	    "__" + _n_offset_recv->name() + "->reset", _n_offset_recv->name(), "0");
	_n_offset_recv_reinit->setup();

	// Create a tool to reinit the mask
	std::ostringstream valstr;
	int mpi_rank, mpi_size;
	mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
	assert(mpi_rank >= 0);

	valstr << mpi_rank;
	_mask_reinit =
	    new Set("__" + _mask->name() + "->reset", _mask->name(), valstr.str());
	_mask_reinit->setup();

	// Allocate the host memory to receive the fields and subsequently upload
	// them
	for (auto field : _fields) {
		void* data = malloc(field->size());
		if (!data) {
			std::stringstream msg;
			msg << "Failure allocating \"" << field->size()
			    << "\" bytes for the array \"" << field->name()
			    << "\" in the tool \"" << name() << "\"" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		_fields_recv.push_back(data);
	}

	// Create the receivers
	for (auto proc : _procs) {
		Receiver* receiver = new Receiver(name(),
		                                  varPrefix(),
		                                  _mask,
		                                  _fields,
		                                  _fields_recv,
		                                  proc,
		                                  _n_offset_recv);
		if (!receiver) {
			std::stringstream msg;
			msg << "Failure Allocating memory for the process " << proc
			    << " receiver in tool \"" << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		_receivers.push_back(receiver);
	}
}

MPISync::Exchanger::Exchanger(
	const std::string tool_name,
	const std::string vars_prefix,
	InputOutput::ArrayVariable* mask,
	const std::vector<InputOutput::ArrayVariable*> fields,
	const std::vector<void*> field_hosts,
	const unsigned int proc)
  : _name(tool_name)
  , _var_prefix(vars_prefix)
  , _mask(mask)
  , _fields(fields)
  , _fields_host(field_hosts)
  , _proc(proc)
  , _n(0)
{
	_n = _mask->size() / InputOutput::Variables::typeToBytes(_mask->type());
}

const MPISync::Exchanger::MPIType
MPISync::Exchanger::typeToMPI(std::string t)
{
	MPISync::Exchanger::MPIType mpi_t;

	mpi_t.n = 1;
	mpi_t.t = MPI_DATATYPE_NULL;

	if (t.back() == '*') {
		t.pop_back();
	}
	if (hasSuffix(t, "vec")) {
#ifdef HAVE_3D
		mpi_t.n = 4;
#else
		mpi_t.n = 2;
#endif
	} else if (hasSuffix(t, "2")) {
		mpi_t.n = 2;
		t.pop_back();
	} else if (hasSuffix(t, "3")) {
		mpi_t.n = 3;
		t.pop_back();
	} else if (hasSuffix(t, "4")) {
		mpi_t.n = 4;
		t.pop_back();
	}

	if ((!t.compare("int")) || (!t.compare("ivec"))) {
		mpi_t.t = MPI_INT;
	} else if ((!t.compare("unsigned int")) || (!t.compare("uivec"))) {
		mpi_t.t = MPI_UNSIGNED;
	} else if ((!t.compare("float")) || (!t.compare("vec"))) {
		mpi_t.t = MPI_FLOAT;
	}

	return mpi_t;
}

MPISync::Sender::Sender(const std::string name,
                        const std::string vars_prefix,
                        InputOutput::ArrayVariable* mask,
                        const std::vector<InputOutput::ArrayVariable*> fields,
                        const std::vector<void*> field_hosts,
                        const unsigned int proc)
  : MPISync::Exchanger(name, vars_prefix, mask, fields, field_hosts, proc)
  , _n_offset(NULL)
  , _n_offset_mask(NULL)
  , _n_offset_kernel(NULL)
  , _n_offset_reduction(NULL)
  , _n_send(NULL)
  , _n_send_mask(NULL)
  , _n_send_kernel(NULL)
  , _n_send_reduction(NULL)
  , _global_work_size(0)
  , _local_work_size(0)
{
	setupSubMaskMems();
	setupOpenCL("n_offset_mask");
	setupReduction("n_offset");
	setupOpenCL("n_send_mask");
	setupReduction("n_send");
}

MPISync::Sender::~Sender()
{
	if (_n_offset_kernel)
		clReleaseKernel(_n_offset_kernel);
	_n_offset_kernel = NULL;
	if (_n_offset_reduction)
		delete _n_offset_reduction;
	_n_offset_reduction = NULL;
	if (_n_send_kernel)
		clReleaseKernel(_n_send_kernel);
	_n_send_kernel = NULL;
	if (_n_send_reduction)
		delete _n_send_reduction;
	_n_send_reduction = NULL;
}

typedef struct
{
	CalcServer* C;
	InputOutput::ArrayVariable* field;
	void* ptr;
	unsigned int proc;
	unsigned int* offset;
	unsigned int* n;
	int tag;
	cl_event user_event;
} MPISyncSendUserData;

void CL_CALLBACK
cbMPISend(cl_event n_event, cl_int cmd_exec_status, void* user_data)
{
	MPI_Request req;
	MPISyncSendUserData* data = (MPISyncSendUserData*)user_data;
	unsigned int offset = *(data->offset);
	unsigned int n = *(data->n);

	if (data->tag == 1) {
		Aqua::MPI::isend(&n, 1, MPI_UNSIGNED, data->proc, 0, MPI_COMM_WORLD);
	}

	if (!n) {
		free(data);
		return;
	}

	InputOutput::ArrayVariable* field = data->field;
	const size_t tsize = InputOutput::Variables::typeToBytes(field->type());
	void* ptr = (void*)((char*)(data->ptr) + offset * tsize);

	// Download the data
	cl_int err_code;
	err_code = clEnqueueReadBuffer(data->C->command_queue(cmd_queue_mpi),
	                               *(cl_mem*)field->get(),
	                               CL_TRUE,
	                               offset * tsize,
	                               n * tsize,
	                               ptr,
	                               0,
	                               NULL,
	                               NULL);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure downloading the variable \"") +
	                       field->name() + "\".");
	err_code = clSetUserEventStatus(data->user_event, CL_COMPLETE);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure setting as completed the reading event of \"") +
	        field->name() + "\".");
	err_code = clReleaseEvent(data->user_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the reading event of \"") +
	        field->name() + "\".");

	// Select the appropriate datatype (MPI needs its own type decriptor), and
	// addapt the array length (vectorial types)
	const MPISync::Exchanger::MPIType mpi_t =
	    MPISync::Exchanger::typeToMPI(field->type());
	if (mpi_t.t == MPI_DATATYPE_NULL) {
		std::ostringstream msg;
		msg << "Unrecognized type \"" << field->type() << "\" for variable \""
		    << field->name() << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		free(data);
		throw std::runtime_error("Invalid type");
	}

	// Launch the missiles. Again we can proceed in synchronous mode
	Aqua::MPI::isend(ptr, n * mpi_t.n, mpi_t.t, data->proc, data->tag,
	                 MPI_COMM_WORLD);
	free(data);
}

void
MPISync::Sender::execute(EventProfile* profiler)
{
	unsigned int i, j;
	cl_int err_code;
	cl_event event;
	std::vector<cl_event> event_wait_list;
	CalcServer* C = CalcServer::singleton();

	// Compute the offset of the first particle to be sent
	event_wait_list = { _mask->getWritingEvent(),
		                _n_offset_mask->getWritingEvent() };
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _n_offset_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  &_local_work_size,
	                                  event_wait_list.size(),
	                                  event_wait_list.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing the tool \"") + name() +
	                       "\".");
	profiler->start(event);

	// Mark the variables with an event to be subsequently waited for
	_mask->addReadingEvent(event);
	_n_offset_mask->setWritingEvent(event);
	_n_offset_reduction->execute();

	// Compute number of particles to be sent
	event_wait_list = { _mask->getWritingEvent(),
		                _n_send_mask->getWritingEvent() };
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _n_send_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  &_local_work_size,
	                                  event_wait_list.size(),
	                                  event_wait_list.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing the tool \"") + name() +
	                       "\".");
	// Mark the variables with an event to be subsequently waited for
	_mask->addReadingEvent(event);
	_n_send_mask->setWritingEvent(event);
	_n_send_reduction->execute();
	profiler->end(_n_send->getWritingEvent());

	for (unsigned int i = 0; i < _fields.size(); i++) {
		// Setup a events syncing point to can register a callback where the
		// fields can be sent
		event_wait_list = { _n_offset->getWritingEvent(),
			                _n_send->getWritingEvent(),
			                _fields.at(i)->getWritingEvent() };

		err_code = clEnqueueMarkerWithWaitList(C->command_queue(),
		                                       event_wait_list.size(),
		                                       event_wait_list.data(),
		                                       &event);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string(
		        "Failure creating send events syncing point in tool \"") +
		        name() + "\".");

		// Setup the data to pass to the callback
		MPISyncSendUserData* user_data =
		    (MPISyncSendUserData*)malloc(sizeof(MPISyncSendUserData));
		if (!user_data) {
			std::stringstream msg;
			msg << "Failure allocating " << sizeof(MPISyncSendUserData)
			    << " bytes for the variable \"" << _fields.at(i)->name()
			    << "\" user data" << std::endl;
			throw std::bad_alloc();
		}
		user_data->C = C;
		user_data->field = _fields.at(i);
		user_data->ptr = _fields_host.at(i);
		user_data->proc = _proc;
		user_data->offset = (unsigned int*)(_n_offset->get_async());
		user_data->n = (unsigned int*)(_n_send->get_async());
		user_data->tag = i + 1;
		user_data->user_event = clCreateUserEvent(C->context(), &err_code);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure creating the download & send event for \"") +
		        _fields.at(i)->name() + "\" in tool \"" + name() + "\".");
		user_data->field->addReadingEvent(user_data->user_event);

		// So we can asynchronously ask to dispatch the data send
		err_code = clSetEventCallback(
		    event, CL_COMPLETE, &cbMPISend, (void*)(user_data));
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure setting the download & send callback for \"") +
		        _fields.at(i)->name() + "\" in tool \"" + name() + "\".");
		err_code = clFlush(C->command_queue());
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure flushing the command queue for \"") +
		        _fields.at(i)->name() + "\" in tool \"" + name() + "\".");
	}
}

void
MPISync::Sender::setupSubMaskMems()
{
	unsigned int i;
	std::ostringstream name, valstr;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();
	valstr << _n;

	i = 0;
	name << _var_prefix << _mask->name() << "_n_offset_mask_" << i;
	while (vars->get(name.str()) != NULL) {
		name << _var_prefix << _mask->name() << "_n_offset_mask_" << ++i;
	}
	vars->registerVariable(name.str(), "unsigned int*", valstr.str(), "");
	_n_offset_mask = (InputOutput::ArrayVariable*)vars->get(name.str());

	i = 0;
	name.str("");
	name << _var_prefix << _mask->name() << "_n_send_mask_" << i;
	while (vars->get(name.str()) != NULL) {
		name << _var_prefix << _mask->name() << "_n_send_mask_" << ++i;
	}
	if (!vars->get(name.str()))
		vars->registerVariable(name.str(), "unsigned int*", valstr.str(), "");
	_n_send_mask = (InputOutput::ArrayVariable*)vars->get(name.str());
}

void
MPISync::Sender::setupOpenCL(const std::string kernel_name)
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();
	InputOutput::ArrayVariable* submask;

	if (!kernel_name.compare("n_offset_mask")) {
		submask = _n_offset_mask;
	} else if (!kernel_name.compare("n_send_mask")) {
		submask = _n_send_mask;
	} else {
		LOG(L_ERROR, std::string("Bad kernel name \"") + kernel_name + "\"\n");
		throw std::runtime_error("Invalid value");
	}

	std::ostringstream source;
	source << MPISYNC_INC << MPISYNC_SRC;
	cl_kernel kernel = compile_kernel(source.str(), kernel_name);

	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying the work group size in tool \"") +
	        name() + "\".");
	if (_local_work_size < __CL_MIN_LOCALSIZE__) {
		std::stringstream msg;
		LOG(L_ERROR, "UnSort cannot be performed.\n");
		msg << "\t" << _local_work_size
		    << " elements can be executed, but __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}

	_global_work_size = roundUp<size_t>(_n, _local_work_size);
	err_code = clSetKernelArg(kernel, 0, _mask->typesize(), _mask->get());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the mask argument in tool \"") + name() +
	        "\".");
	err_code = clSetKernelArg(kernel, 1, submask->typesize(), submask->get());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the submask argument in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(kernel, 2, sizeof(unsigned int), (void*)&_proc);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the proc argument in tool \"") + name() +
	        "\".");
	err_code = clSetKernelArg(kernel, 3, sizeof(unsigned int), (void*)&_n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the array size argument in tool \"") +
	        name() + "\".");

	if (!kernel_name.compare("n_offset_mask")) {
		_n_offset_kernel = kernel;
	} else if (!kernel_name.compare("n_send_mask")) {
		_n_send_kernel = kernel;
	}
}

void
MPISync::Sender::setupReduction(const std::string var_name)
{
	unsigned int i;
	std::ostringstream name;
	InputOutput::ArrayVariable* submask;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	if (!var_name.compare("n_offset")) {
		submask = _n_offset_mask;
	} else if (!var_name.compare("n_send")) {
		submask = _n_send_mask;
	} else {
		LOG(L_ERROR, std::string("Bad variable name \"") + var_name + "\"\n");
		throw std::runtime_error("Invalid value");
	}

	// Register a variable where we can reduce the result
	i = 0;
	name << _var_prefix << var_name << "_" << i;
	while (vars->get(name.str()) != NULL) {
		name << _var_prefix << var_name << "_" << ++i;
	}
	if (!vars->get(name.str()))
		vars->registerVariable(name.str(), "unsigned int", "", "0");
	InputOutput::Variable* var = vars->get(name.str());
	if (!var_name.compare("n_offset")) {
		_n_offset = (InputOutput::UIntVariable*)var;
	} else if (!var_name.compare("n_send")) {
		_n_send = (InputOutput::UIntVariable*)var;
	}

	name << "->Sum";
	std::string op = "c = a + b;\n";
	Reduction* reduction =
	    new Reduction(name.str(), submask->name(), var->name(), op, "0");
	reduction->setup();

	if (!var_name.compare("n_offset")) {
		_n_offset_reduction = reduction;
	} else if (!var_name.compare("n_send")) {
		_n_send_reduction = reduction;
	}
}

MPISync::Receiver::Receiver(
    const std::string name,
    const std::string vars_prefix,
    InputOutput::ArrayVariable* mask,
    const std::vector<InputOutput::ArrayVariable*> fields,
    const std::vector<void*> field_hosts,
    const unsigned int proc,
    InputOutput::UIntVariable* n_offset)
  : MPISync::Exchanger(name, vars_prefix, mask, fields, field_hosts, proc)
  , _kernel(NULL)
  , _n_offset(n_offset)
  , _local_work_size(0)
{
	setupOpenCL();
}

MPISync::Receiver::~Receiver()
{
	if (_kernel)
		clReleaseKernel(_kernel);
	_kernel = NULL;
}

typedef struct
{
	CalcServer* C;
	size_t n_fields;
	InputOutput::ArrayVariable** fields;
	void** ptrs;
	unsigned int proc;
	InputOutput::UIntVariable* offset;
	// To set the mask
	InputOutput::ArrayVariable* mask;
	cl_kernel kernel;
	size_t local_work_size;
	// We shall store the events to avoid someone change them while the
	// callback has not already processed the data
	cl_event* field_events;
	cl_event offset_event;
	cl_event mask_event;
	// Store the profiler so we can feed it up on the callback
	EventProfile* profiler;
} MPISyncRecvUserData;

void CL_CALLBACK
cbMPIRecv(cl_event n_event, cl_int cmd_exec_status, void* user_data)
{
	cl_int err_code;
	MPI_Request req;
	MPI_Status status;
	cl_event mask_event = NULL, field_event = NULL;
	MPISyncRecvUserData* data = (MPISyncRecvUserData*)user_data;
	unsigned int offset = *(unsigned int*)data->offset->get_async();
	unsigned int n;

	// We need to receive the number of transmitted elements
	Aqua::MPI::recv(&n, 1, MPI_UNSIGNED, data->proc, 0, MPI_COMM_WORLD);

	// So we can set the offset for the next receivers
	unsigned int next_offset = offset + n;
	data->offset->set_async((void*)(&offset));
	err_code = clSetUserEventStatus(data->offset_event, CL_COMPLETE);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure setting user event status for \"") +
	                       data->offset->name() + "\" variable.");
	err_code = clReleaseEvent(data->offset_event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure releasing user event for \"") +
	                       data->offset->name() + "\" variable.");

	if (!n) {
		err_code = clSetUserEventStatus(data->mask_event, CL_COMPLETE);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure setting user event status for \"") +
		        data->mask->name() + "\" variable.");
		err_code = clReleaseEvent(data->mask_event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure releasing user event for \"") +
		                       data->mask->name() + "\" variable.");
		for (unsigned int i = 0; i < data->n_fields; i++) {
			InputOutput::ArrayVariable* field = data->fields[i];
			err_code = clSetUserEventStatus(data->field_events[i], CL_COMPLETE);
			CHECK_OCL_OR_THROW(
			    err_code,
			    std::string("Failure setting user event status for \"") +
			        field->name() + "\" variable.");
			err_code = clReleaseEvent(data->field_events[i]);
			CHECK_OCL_OR_THROW(
			    err_code,
			    std::string("Failure releasing user event status for \"") +
			        field->name() + "\" variable.");
		}
		free(data);

		return;
	}

	// We can now set the mask values
	err_code =
	    clSetKernelArg(data->kernel, 2, sizeof(unsigned int), (void*)&offset);
	CHECK_OCL_OR_THROW(err_code, "Failure sending the offset argument");
	err_code = clSetKernelArg(data->kernel, 3, sizeof(unsigned int), (void*)&n);
	CHECK_OCL_OR_THROW(err_code, "Failure sending the n argument");
	size_t global_work_size = roundUp<size_t>(n, data->local_work_size);
	err_code = clEnqueueNDRangeKernel(data->C->command_queue(cmd_queue_mpi),
	                                  data->kernel,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &(data->local_work_size),
	                                  0,
	                                  NULL,
	                                  &mask_event);
	CHECK_OCL_OR_THROW(err_code, "Failure setting the mask");
	data->profiler->start(mask_event);

	// Synchronize the mask writing event with the locking one
	sync_user_event(data->mask_event, mask_event);

	for (unsigned int i = 0; i < data->n_fields; i++) {
		InputOutput::ArrayVariable* field = data->fields[i];
		const size_t tsize = InputOutput::Variables::typeToBytes(field->type());
		void* ptr = (void*)((char*)(data->ptrs[i]) + offset * tsize);

		// Select the appropriate datatype (MPI needs its own type decriptor),
		// and addapt the array length (vectorial types)
		const MPISync::Exchanger::MPIType mpi_t =
		    MPISync::Exchanger::typeToMPI(field->type());
		if (mpi_t.t == MPI_DATATYPE_NULL) {
			std::ostringstream msg;
			msg << "Unrecognized type \"" << field->type()
			    << "\" for variable \"" << field->name() << "\"" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("MPI receiving error");
		}

		// Get the data in synchronous mode
		Aqua::MPI::recv(ptr, n * mpi_t.n, mpi_t.t, data->proc, i + 1,
		                MPI_COMM_WORLD);

		// But upload it in asynchronous mode
		err_code = clEnqueueWriteBuffer(data->C->command_queue(cmd_queue_mpi),
		                                *(cl_mem*)field->get(),
		                                CL_FALSE,
		                                offset * tsize,
		                                n * tsize,
		                                ptr,
		                                0,
		                                NULL,
		                                &field_event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure uploading the variable \"") +
		                       field->name() + "\".");

		sync_user_event(data->field_events[i], field_event);
	}
	data->profiler->end(field_event);

	free(data->field_events);
	free(data);
}

void
MPISync::Receiver::execute(EventProfile* profiler)
{
	unsigned int i;
	cl_int err_code;
	cl_event event;
	auto C = CalcServer::singleton();

	// Setup a syncing event to trigger the callback
	std::vector<cl_event> event_wait_list = _n_offset->getReadingEvents();
	event_wait_list.push_back(_n_offset->getWritingEvent());
	for (auto field : _fields) {
		for (auto e : field->getReadingEvents())
			event_wait_list.push_back(e);
		event_wait_list.push_back(field->getWritingEvent());
	}

	err_code = clEnqueueMarkerWithWaitList(C->command_queue(),
	                                       event_wait_list.size(),
	                                       event_wait_list.data(),
	                                       &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure creating recv events syncing point for tool \"") +
	        name() + "\".");

	// Setup the data to pass to the callback
	MPISyncRecvUserData* user_data =
	    (MPISyncRecvUserData*)malloc(sizeof(MPISyncRecvUserData));
	if (!user_data) {
		std::stringstream msg;
		msg << "Failure allocating " << sizeof(MPISyncSendUserData)
		    << " bytes for user data structure" << std::endl;
		throw std::bad_alloc();
	}

	user_data->C = C;
	user_data->n_fields = _fields.size();
	user_data->fields = _fields.data();
	user_data->ptrs = _fields_host.data();
	user_data->proc = _proc;
	user_data->offset = _n_offset;
	user_data->mask = _mask;
	user_data->kernel = _kernel;
	user_data->local_work_size = _local_work_size;
	user_data->field_events =
	    (cl_event*)malloc(_fields.size() * sizeof(cl_event));
	if (!user_data->field_events) {
		std::stringstream msg;
		msg << "Failure allocating " << _fields.size() * sizeof(cl_event)
		    << " bytes for user data structure" << std::endl;
		throw std::bad_alloc();
	}
	user_data->profiler = profiler;

	// We specifically want to lock some members until they (or their host
	// mirrors) are ready
	for (i = 0; i < _fields.size() + 2; i++) {
		cl_event user_event = clCreateUserEvent(C->context(), &err_code);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure creating recv user event in tool \"") +
		        name() + "\".");
		if (i < _fields.size()) {
			// Lock the field to be uploaded until the data has been correctly
			// received
			_fields.at(i)->setWritingEvent(user_event);
			user_data->field_events[i] = user_event;
		}
		if (i == _fields.size()) {
			// Lock the following receivers until the offset is known
			_n_offset->setWritingEvent(user_event);
			user_data->offset_event = user_event;
		} else {
			// Lock the mask write by further receivers until this instance has
			// wrote its part
			_mask->setWritingEvent(user_event);
			user_data->mask_event = user_event;
		}
	}

	// So we can asynchronously ask to dispatch the data reception
	err_code =
	    clSetEventCallback(event, CL_COMPLETE, &cbMPIRecv, (void*)(user_data));
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure setting the download & receive callback for process") +
	        std::to_string(_proc) + " in tool \"" + name() + "\"");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure flushing the command queue for process") +
			std::to_string(_proc) + " in tool \"" + name() + "\".");
}

void
MPISync::Receiver::setupOpenCL()
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	std::ostringstream source;
	source << MPISYNC_INC << MPISYNC_SRC;
	_kernel = compile_kernel(source.str(), "set_mask");

	err_code = clGetKernelWorkGroupInfo(_kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_local_work_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying the work group size in tool \"") +
	        name() + "\"");
	if (_local_work_size < __CL_MIN_LOCALSIZE__) {
		std::stringstream msg;
		LOG(L_ERROR, "UnSort cannot be performed.\n");
		msg << "\t" << _local_work_size
		    << " elements can be executed, but __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}

	err_code = clSetKernelArg(_kernel, 0, _mask->typesize(), _mask->get());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the mask argument in tool \"") + name() +
	        "\"");
	err_code = clSetKernelArg(_kernel, 1, sizeof(unsigned int), (void*)&_proc);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the proc argument in tool \"") + name() +
	        "\"");
}

}
} // namespaces

#endif
