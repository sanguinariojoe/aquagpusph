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

#include <sphPrerequisites.h>

#ifdef HAVE_MPI

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer/MPISync.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/MPISync.hcl"
#include "CalcServer/MPISync.cl"
#endif
std::string MPISYNC_INC = xxd2string(MPISync_hcl_in, MPISync_hcl_in_len);
std::string MPISYNC_SRC = xxd2string(MPISync_cl_in, MPISync_cl_in_len);

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
    , _n_offset(NULL)
{
    int mpi_rank, mpi_size;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
        mpi_size = MPI::COMM_WORLD.Get_size();
    } catch(MPI::Exception e){
        std::ostringstream msg;
        msg << "Error getting MPI rank and size. " << std::endl
            << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    assert(mpi_rank >= 0);
    assert(mpi_size > 0);

    // Check that we have a good list of processes
    if(_procs.size() == 0) {
        for(unsigned int proc=0; proc < mpi_size; proc++) {
            _procs.push_back(proc);
        }
    }
    for(unsigned int i = 0; i < _procs.size(); i++) {
        if((_procs.at(i) == mpi_rank) || (_procs.at(i) >= mpi_size)) {
            _procs.erase(_procs.begin() + i);
        }
    }
}

MPISync::~MPISync()
{
    if(_sort) delete _sort; _sort=NULL;
    if(_n_offset_reinit) delete _n_offset_reinit; _n_offset_reinit=NULL;
    if(_mask_reinit) delete _mask_reinit; _mask_reinit=NULL;

    // We must unset the inner memory objects of the sorted fields to avoid
    // releasing them twice while deleting the sorters
    cl_mem inner_mem = NULL;
    for(auto field : _fields_sorted) {
        field->set((void*)(&inner_mem));
    }
    // Now we can delete the sorters
    for(auto sorter : _field_sorters) {
        delete sorter;
    }
    // And the senders
    for(auto field : _fields_send) {
        free(field);
    }
    for(auto sender : _senders) {
        delete sender;
    }    
    // And the receivers
    for(auto field : _fields_recv) {
        free(field);
    }
    for(auto receiver : _receivers) {
        delete receiver;
    }    
}

void MPISync::setup()
{
    // Get the involved variables
    variables();
    // Setup the mask sorting subtool
    try {
        setupSort();
    } catch (std::runtime_error &e) {
        std::stringstream msg;
        msg << "Error setting up the sorter for tool \"" << name() << "\""
            << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    // Setup the field sorting subtools
    for(auto field : _fields) {
        setupFieldSort(field);
    }

    setupSenders();
    setupReceivers();
}

cl_event MPISync::_execute(const std::vector<cl_event> events)
{
    int mpi_rank, mpi_size;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
        mpi_size = MPI::COMM_WORLD.Get_size();
    } catch(MPI::Exception e){
        std::ostringstream msg;
        msg << "Error getting MPI rank and size. " << std::endl
            << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    assert(mpi_rank >= 0);
    assert(mpi_size > 0);



    if(!_procs.size())
        return NULL;

    // Sort the mask
    std::cout << "*** [" << mpi_rank << "] "
              << "Sort mask..."
              << std::endl;
    try {
        _sort->execute();
    } catch (std::runtime_error &e) {
        std::stringstream msg;
        msg << "Error while sorting the mask for tool \"" << name() << "\""
            << std::endl;
        LOG(L_ERROR, msg.str());
        throw;        
    }
    // Sort the fields
    std::cout << "*** [" << mpi_rank << "] "
              << "Sort fields..."
              << std::endl;
    for(auto sorter : _field_sorters) {
        try {
            sorter->execute();
        } catch (std::runtime_error &e) {
            std::stringstream msg;
            msg << "Error while sorting \"" << sorter->input()->name()
                << "\" for tool \"" << name() << "\"" << std::endl;
            LOG(L_ERROR, msg.str());
            throw;        
        }
    }

    // Send the data to other processes
    std::cout << "*** [" << mpi_rank << "] "
              << "Reinit offset..."
              << std::endl;
    _n_offset_reinit->execute();
    CalcServer *C = CalcServer::singleton();
    clFinish(C->command_queue());
    std::cout << "*** [" << mpi_rank << "] "
              << "Sending..."
              << std::endl;
    for(auto sender : _senders) {
        sender->execute();
    }

    // Receive the data from the other processes
    std::cout << "*** [" << mpi_rank << "] "
              << "Reinit offset..."
              << std::endl;
    _n_offset_reinit->execute();
    std::cout << "*** [" << mpi_rank << "] "
              << "Reinit mask..."
              << std::endl;
    _mask_reinit->execute();
    std::cout << "*** [" << mpi_rank << "] "
              << "Receiving..."
              << std::endl;    
    for(auto receiver : _receivers) {
        receiver->execute();
    }
    std::cout << "*** [" << mpi_rank << "] "
              << "Waiting..."
              << std::endl;    

    return NULL;
}

void MPISync::variables()
{
    size_t n;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    if(!vars->get(_mask_name)) {
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared variable \"" <<
            _mask_name << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_mask_name)->type().compare("unsigned int*")) {
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the variable \"" << _mask_name
            << "\", which has an invalid type" << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\t\"unsigned int*\" was expected, but \""
            << vars->get(_mask_name)->type() << "\" was found." << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _mask = (InputOutput::ArrayVariable *)vars->get(_mask_name);
    _n = _mask->size() / InputOutput::Variables::typeToBytes(_mask->type());

    for(auto var_name : _field_names) {
        if(!vars->get(var_name)){
            std::stringstream msg;
            msg << "The tool \"" << name()
                << "\" is asking the undeclared variable \""
                << var_name << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable");
        }
        if(vars->get(var_name)->type().find('*') == std::string::npos){
            std::stringstream msg;
            msg << "The tool \"" << name()
                << "\" may not use a scalar variable (\""
                << var_name << "\")." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
        InputOutput::ArrayVariable *var =
            (InputOutput::ArrayVariable *)vars->get(var_name);
        n = var->size() / InputOutput::Variables::typeToBytes(var->type());
        if(n != _n) {
            std::stringstream msg;
            msg << "Wrong variable length in the tool \"" << name()
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            msg.str("");
            msg << "\t\"" << _mask_name << "\" has length " << _n << std::endl;
            LOG0(L_DEBUG, msg.str());
            msg.str("");
            msg << "\t\"" << var_name << "\" has length " << n << std::endl;
            LOG0(L_DEBUG, msg.str());
            throw std::runtime_error("Invalid variable length");
        }
        _fields.push_back((InputOutput::ArrayVariable *)vars->get(var_name));
    }

    std::vector<InputOutput::Variable*> deps;
    for(auto field : _fields) {
        deps.push_back((InputOutput::Variable*)field);
    }
    deps.push_back((InputOutput::Variable*)_mask);
    setDependencies(deps);
}

void MPISync::setupSort()
{
    std::string var_name;
    std::ostringstream valstr;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Register the variables to store the permutations
    valstr << _n;
    var_name = "__" + _mask_name + "_unsorted";
    vars->registerVariable(var_name, "unsigned int*", valstr.str(), "");
    _unsorted_id = (InputOutput::ArrayVariable*)vars->get(var_name);
    var_name = "__" + _mask_name + "_sorted";
    vars->registerVariable(var_name, "unsigned int*", valstr.str(), "");
    _sorted_id = (InputOutput::ArrayVariable*)vars->get(var_name);

    // Create the sorter
    var_name = "__" + _mask_name + "->Radix-Sort";
    _sort = new RadixSort(var_name,
                          _mask_name,
                          _unsorted_id->name(),
                          _sorted_id->name());
    _sort->setup();
}

void MPISync::setupFieldSort(InputOutput::ArrayVariable* field)
{
    std::string var_name;
    std::ostringstream valstr;
    cl_int err_code;
    cl_mem inner_mem;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Register the variable to store the sorted copy
    valstr.str(""); valstr << _n;
    var_name = "__" + field->name() + "_sorted";
    vars->registerVariable(var_name, field->type(), valstr.str(), "");
    _fields_sorted.push_back((InputOutput::ArrayVariable*)vars->get(var_name));
    // Remove the inner memory object, since we are using the one computed by
    // the sorting tool
    inner_mem = *(cl_mem*)_fields_sorted.back()->get();
    err_code = clReleaseMemObject(inner_mem);
    if(err_code != CL_SUCCESS) {
        std::stringstream msg;
        msg << "Failure Releasing the inner memory object of \""
            << _fields_sorted.back()->name() << "\" for tool \"" << name() 
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Create the sorter
    var_name = "__" + field->name() + "->Radix-Sort";
    _field_sorters.push_back(
        new UnSort(var_name, field->name(), _sorted_id->name()));
    _field_sorters.back()->setup();

    // Replace the inner CL memory object of the sorted field by the one
    // computed by the sorter
    inner_mem = _field_sorters.back()->output();
    _fields_sorted.back()->set((void*)(&inner_mem));
}

void MPISync::setupSenders()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a variable for the cumulative offset computation
    vars->registerVariable("__mpi_offset", "unsigned int", "", "0");
    _n_offset = (InputOutput::UIntVariable*)vars->get("__mpi_offset");

    // Create a tool to reinit it to zero value
    _n_offset_reinit = new SetScalar("__" + _n_offset->name() + "->reset",
                                     _n_offset->name(),
                                     "0");
    _n_offset_reinit->setup();

    // Allocate the host memory to download the fields and subsequently send
    // them
    for(auto field : _fields) {
        void *data = malloc(field->size());
        if(!data) {
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
    for(auto proc : _procs) {
        Sender* sender = new Sender(name(),
                                    _mask,
                                    _fields_sorted,
                                    _fields_send,
                                    proc,
                                    _n_offset);
        if(!sender) {
            std::stringstream msg;
            msg << "Failure Allocating memory for the process "
                << proc << " sender in tool \"" << name() 
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::bad_alloc();
        }
        _senders.push_back(sender);
    }
}

void MPISync::setupReceivers()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a tool to reinit the mask
    std::ostringstream valstr;
    int mpi_rank, mpi_size;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
    } catch(MPI::Exception e){
        std::ostringstream msg;
        msg << "Error getting MPI rank. " << std::endl
            << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    assert(mpi_rank >= 0);

    valstr << mpi_rank;
    _mask_reinit = new Set("__" + _mask->name() + "->reset",
                           _mask->name(),
                           valstr.str());
    _mask_reinit->setup();
    
    // Allocate the host memory to receive the fields and subsequently upload
    // them
    for(auto field : _fields) {
        void *data = malloc(field->size());
        if(!data) {
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
    for(auto proc : _procs) {
        Receiver* receiver = new Receiver(name(),
                                          _mask,
                                          _fields,
                                          _fields_recv,
                                          proc,
                                          _n_offset);
        if(!receiver) {
            std::stringstream msg;
            msg << "Failure Allocating memory for the process "
                << proc << " sender in tool \"" << name() 
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::bad_alloc();
        }
        _receivers.push_back(receiver);
    }
}

MPISync::Exchanger::Exchanger(const std::string tool_name,
                              InputOutput::ArrayVariable *mask,
                              const std::vector<InputOutput::ArrayVariable*> fields,
                              const std::vector<void*> field_hosts,
                              const unsigned int proc)
    : _name(tool_name)
    , _mask(mask)
    , _fields(fields)
    , _fields_host(field_hosts)
    , _proc(proc)
    , _n(0)
{
    _n = _mask->size() / InputOutput::Variables::typeToBytes(_mask->type());
}

MPISync::Exchanger::~Exchanger()
{
}

const MPISync::Exchanger::MPIType MPISync::Exchanger::typeToMPI(std::string t)
{
    MPISync::Exchanger::MPIType mpi_t;

    mpi_t.n = 1;
    mpi_t.t = MPI::DATATYPE_NULL;

    if((t.back() == '*')){
        t.pop_back();
    }
    if(hasSuffix(t, "vec")) {
#ifdef HAVE_3D
        mpi_t.n = 4;
#else
        mpi_t.n = 2;
#endif
    } else if(hasSuffix(t, "2")) {
        mpi_t.n = 2;
        t.pop_back();
    } else if(hasSuffix(t, "3")) {
        mpi_t.n = 3;
        t.pop_back();
    } else if(hasSuffix(t, "4")) {
        mpi_t.n = 4;
        t.pop_back();
    }

    if((!t.compare("int")) || (!t.compare("ivec"))) {
        mpi_t.t = MPI::INT;
    } else if((!t.compare("unsigned int")) || (!t.compare("uivec"))) {
        mpi_t.t = MPI::UNSIGNED;
    } else if((!t.compare("float")) || (!t.compare("vec"))) {
        mpi_t.t = MPI::FLOAT;
    }

    return mpi_t;
}

MPISync::Sender::Sender(const std::string name,
                        InputOutput::ArrayVariable *mask,
                        const std::vector<InputOutput::ArrayVariable*> fields,
                        const std::vector<void*> field_hosts,
                        const unsigned int proc,
                        InputOutput::UIntVariable *n_offset)
    : MPISync::Exchanger(name, mask, fields, field_hosts, proc)
    , _submask(NULL)
    , _kernel(NULL)
    , _n_send(NULL)
    , _n_offset(n_offset)
    , _global_work_size(0)
    , _local_work_size(0)
{
    setupSubMaskMem();
    setupOpenCL();
    setupNSend();
}

MPISync::Sender::~Sender()
{
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
    if(_n_send_reduction) delete _n_send_reduction; _n_send_reduction=NULL;
}

typedef struct {
    CalcServer *C;
    InputOutput::ArrayVariable *field;
    InputOutput::UIntVariable *n;
    bool send_n;
    unsigned int proc;
    InputOutput::UIntVariable *offset;
    void* ptr;
    unsigned int tag;
} MPISyncSendUserData;

void CL_CALLBACK cbDownloadAndSend(cl_event n_event,
                                   cl_int cmd_exec_status,
                                   void *user_data)
{
    int mpi_rank, mpi_size;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
        mpi_size = MPI::COMM_WORLD.Get_size();
    } catch(MPI::Exception e){
        std::ostringstream msg;
        msg << "Error getting MPI rank and size. " << std::endl
            << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    assert(mpi_rank >= 0);
    assert(mpi_size > 0);

    cl_int err_code;
    MPISyncSendUserData *data = (MPISyncSendUserData*)user_data;
    unsigned int n = *(unsigned int*)data->n->get_async();
    unsigned int offset = *(unsigned int*)data->offset->get_async();

    if(data->send_n) {
        // We need to report the proc the number of particles to be sent. We
        // don't want to hang here, so we aqsk for async sending
        MPI::COMM_WORLD.Isend(&n, 1, MPI::UNSIGNED, data->proc, 0);
        // We know this is the very only sending instance (to this proc)
        // executing this piece of code, so we can take advantage to increase
        // the offset variable
        offset += n;
        data->offset->set_async((void*)(&offset));
    }
    // Reached this point, the offset is always n elements ahead
    offset -= n;

    if(!n) {
        // There is nothing to send
        std::cout << "############# " << mpi_rank << " " << data->field->name() << " no send???" << std::endl;
        err_code = clSetUserEventStatus(data->field->getEvent(), CL_COMPLETE);
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure setting user event status for \"" <<
                data->field->name() << "\" variable." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        std::cout << "############# " << mpi_rank << " " << data->field->name() << " OK" << std::endl;
        return;
    }

    // We can now proceed to download the data. We are in a parallel thread, so
    // we can do this synchronously
    const size_t typesize = InputOutput::Variables::typeToBytes(
        data->field->type());
    void *ptr = (void*)((char*)data->ptr + offset * typesize);
    err_code = clEnqueueReadBuffer(data->C->command_queue(),
                                   *(cl_mem*)data->field->get(),
                                   CL_TRUE,
                                   offset * typesize,
                                   n * typesize,
                                   ptr,
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure downloading the variable \"" << data->field->name()
            << "\" to send that by MPI." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // If we reach this point, then all the information we need which can be
    // shared with other threads has been already consumed, so we can mark as
    // complete the events affecting both the field and the offset (eventually)
    // variables
    err_code = clSetUserEventStatus(data->field->getEvent(), CL_COMPLETE);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure setting user event status for \"" <<
               data->field->name() << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Select the appropriate datatype (MPI needs its own type decriptor), and
    // addapt the array length (vectorial types)
    const MPISync::Exchanger::MPIType mpi_t = 
        MPISync::Exchanger::typeToMPI(data->field->type());
    if(mpi_t.t == MPI::DATATYPE_NULL) {
        std::ostringstream msg;
        msg << "Unrecognized type \"" << data->field->type()
            << "\" for variable \"" << data->field->name() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        return;
    }

    // Launch the missiles. Again we can proceed in synchronous mode
    MPI::COMM_WORLD.Send(ptr, n * mpi_t.n, mpi_t.t, data->proc, data->tag);
}

void MPISync::Sender::execute()
{
    unsigned int i, j;
    cl_int err_code;
    cl_event event;
    CalcServer *C = CalcServer::singleton();

    int mpi_rank, mpi_size;
    try {
        mpi_rank = MPI::COMM_WORLD.Get_rank();
        mpi_size = MPI::COMM_WORLD.Get_size();
    } catch(MPI::Exception e){
        std::ostringstream msg;
        msg << "Error getting MPI rank and size. " << std::endl
            << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    assert(mpi_rank >= 0);
    assert(mpi_size > 0);


    // Extract the submask
    std::cout << "****** [" << mpi_rank << "] "
              << "Submask..."
              << std::endl;
    const cl_event event_wait = _mask->getEvent();
    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _kernel,
                                      1,
                                      NULL,
                                      &_global_work_size,
                                      &_local_work_size,
                                      1,
                                      &event_wait,
                                      &event);
    if(err_code != CL_SUCCESS) {
        std::stringstream msg;
        msg << "Failure executing the tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Mark the variables with an event to be subsequently waited for
    _mask->setEvent(event);
    _submask->setEvent(event);

    // Compute the number of elements to download
    clFinish(C->command_queue());
    std::cout << "****** [" << mpi_rank << "] "
              << "send reduction..."
              << std::endl;
    _n_send_reduction->execute();

    clFinish(C->command_queue());
    std::cout << "****** [" << mpi_rank << "] "
              << "send callbacks..."
              << std::endl;
    for(i = 0; i < _fields.size(); i++) {
        std::cout << "******** [" << mpi_rank << "] " << "(" << i << "/" << _fields.size() << ") "
                  << "start..."
                  << std::endl;
        // At this point we cannot do anything else but waiting for the
        // reduction to be completed, since we need to know the number of
        // elements to be downloaded and sent.
        // Along this line, we are collecting all the events involved in the
        // process to wait for them, and then we are setting up a new syncing
        // event to trigger the donwloading and send callback execution.
        // However, we should ensure that the field remain untouched until it is
        // downloaded, operation carried out inside the callback itself.
        // Therefore we are feeding the field with a new user event, which we
        // can control.
        // Along the same line, offset shall be computed in strict order, that
        // is we shall grant that the callbacks are executed in the same order
        // they were registered. For this reason we are setting the same custom
        // event to the offset variable
        const cl_event event_wait_list[3] = {
            _n_send->getEvent(),
            _n_offset->getEvent(),
            _fields.at(i)->getEvent()
        };
        // We better create the syncing event ASAP, so we don't need to retain
        // events before setting new ones to _fields.at(i) and _n_offset
        std::cout << "******** [" << mpi_rank << "] " << "(" << _fields.at(i)->name() << ") "
                  << "syncing event..."
                  << std::endl;
        err_code = clEnqueueMarkerWithWaitList(C->command_queue(),
                                               3,
                                               event_wait_list,
                                               &event);
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure creating send events syncing point for field \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }

        cl_event field_event = clCreateUserEvent(C->context(), &err_code);
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure creating send user event for \"" <<
                _fields.at(i)->name() << "\" variable in tool \"" << name()
                << "\"" << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        _fields.at(i)->setEvent(field_event);
        if(i == 0)
            _n_offset->setEvent(field_event);

        MPISyncSendUserData user_data;
        user_data.C = C;
        user_data.field = _fields.at(i);
        user_data.n = _n_send;
        user_data.send_n = i == 0;
        user_data.proc = _proc;
        user_data.offset = _n_offset;
        user_data.ptr = _fields_host.at(i);
        user_data.tag = i + 1;
        // So we can asynchronously ask to dispatch the data download and send
        // when reduction is done.
        std::cout << "******** [" << mpi_rank << "] " << "(" << _fields.at(i)->name() << ") "
                  << "callback..."
                  << std::endl;
        err_code = clSetEventCallback(event,
                                      CL_COMPLETE,
                                      &cbDownloadAndSend,
                                      (void*)(&user_data));
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure setting the download & send callback for \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        clFinish(C->command_queue());
        err_code = clReleaseEvent(field_event);
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure releasing user event for \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
    }
    std::cout << "****** [" << mpi_rank << "] "
              << "OK..."
              << std::endl;
}

void MPISync::Sender::setupSubMaskMem()
{
    unsigned int i;
    std::ostringstream name, valstr;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Register a variable where we can reduce the result
    i = 0;
    name << "__" << _mask->name() << "_submask_" << i;
    while(vars->get(name.str()) != NULL) {
        name << "__" << _mask->name() << "_submask_" << ++i;
    }
    valstr << _n;
    vars->registerVariable(name.str(), "unsigned int*", valstr.str(), "");
    _submask = (InputOutput::ArrayVariable*)vars->get(name.str());
}

void MPISync::Sender::setupOpenCL()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream source;
    source << MPISYNC_INC << MPISYNC_SRC;
    _kernel = compile_kernel(source.str(), "extract_submask");

    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(_local_work_size < __CL_MIN_LOCALSIZE__){
        std::stringstream msg;
        LOG(L_ERROR, "UnSort cannot be performed.\n");
        msg << "\t" << _local_work_size
            << " elements can be executed, but __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
    }

    _global_work_size = roundUp(_n, _local_work_size);
    err_code = clSetKernelArg(_kernel,
                              0,
                              _mask->typesize(),
                              _mask->get());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the mask argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_kernel,
                              1,
                              _submask->typesize(),
                              _submask->get());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the submask argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_kernel,
                              2,
                              sizeof(unsigned int),
                              (void*)&_proc);
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the proc argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_kernel,
                              3,
                              sizeof(unsigned int),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the array size argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
}

void MPISync::Sender::setupNSend()
{
    unsigned int i;
    std::ostringstream name;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Register a variable where we can reduce the result
    i = 0;
    name << "__n_send_" << i;
    while(vars->get(name.str()) != NULL) {
        name << "__n_send_" << ++i;
    }
    vars->registerVariable(name.str(), "unsigned int", "", "0");
    _n_send = (InputOutput::UIntVariable*)vars->get(name.str());

    name << "->Sum";
    std::string op = "c = a + b;\n";
    _n_send_reduction = new Reduction(name.str(),
                                      _submask->name(),
                                      _n_send->name(),
                                      op,
                                      "0");
    _n_send_reduction->setup();
}




MPISync::Receiver::Receiver(const std::string name,
                        InputOutput::ArrayVariable *mask,
                        const std::vector<InputOutput::ArrayVariable*> fields,
                        const std::vector<void*> field_hosts,
                        const unsigned int proc,
                        InputOutput::UIntVariable *n_offset)
    : MPISync::Exchanger(name, mask, fields, field_hosts, proc)
    , _kernel(NULL)
    , _n_offset(n_offset)
    , _n_recv(0)
    , _local_work_size(0)
{
    setupOpenCL();
}

MPISync::Receiver::~Receiver()
{
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
}

typedef struct {
    CalcServer *C;
    InputOutput::ArrayVariable *field;
    bool recv_n;
    unsigned int proc;
    InputOutput::UIntVariable *offset;
    unsigned int *n;
    void* ptr;
    unsigned int tag;
    // To set the mask
    InputOutput::ArrayVariable *mask;
    cl_kernel kernel;
    size_t local_work_size;

} MPISyncRecvUserData;

void CL_CALLBACK cbReceiveAndUpload(cl_event n_event,
                                   cl_int cmd_exec_status,
                                   void *user_data)
{
    cl_int err_code;
    cl_event mask_event=NULL;
    MPISyncRecvUserData *data = (MPISyncRecvUserData*)user_data;
    unsigned int offset = *(unsigned int*)data->offset->get();
    if(data->recv_n) {
        // We need to receive the number of transmitted elements
        MPI::COMM_WORLD.Recv(data->n, 1, MPI::UNSIGNED, data->proc, 0);
        // We know this is the very only sending instance (for this proc)
        // executing this piece of code, so we can take advantage to increase
        // the offset variable
        offset += *(data->n);
        data->offset->set_async((void*)(&offset));
        // We also take advantage to set the mask values
        err_code = clSetKernelArg(data->kernel,
                                  2,
                                  sizeof(unsigned int),
                                  (void*)&offset);
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure sending the offset argument\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        err_code = clSetKernelArg(data->kernel,
                                  3,
                                  sizeof(unsigned int),
                                  (void*)data->n);
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure sending the n argument\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        size_t global_work_size = roundUp(*(data->n), data->local_work_size);
        err_code = clEnqueueNDRangeKernel(data->C->command_queue(),
                                          data->kernel,
                                          1,
                                          NULL,
                                          &global_work_size,
                                          &(data->local_work_size),
                                          0,
                                          NULL,
                                          &mask_event);
        if(err_code != CL_SUCCESS) {
            LOG(L_ERROR, "Failure setting the mask");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
    }
    unsigned int n = *(data->n);
    // Reached this point, the offset is always n elements ahead
    offset -= n;

    if(!n) {
        // There is nothing to receive
        err_code = clSetUserEventStatus(data->field->getEvent(), CL_COMPLETE);
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure setting user event status for \"" <<
                data->field->name() << "\" variable." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        return;
    }

    // Select the appropriate datatype (MPI needs its own type decriptor), and
    // addapt the array length (vectorial types)
    const MPISync::Exchanger::MPIType mpi_t = 
        MPISync::Exchanger::typeToMPI(data->field->type());
    if(mpi_t.t == MPI::DATATYPE_NULL) {
        std::ostringstream msg;
        msg << "Unrecognized type \"" << data->field->type()
            << "\" for variable \"" << data->field->name() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        return;
    }

    // We can now proceed to receive the data. We are in a parallel thread, so
    // we can do this synchronously
    const size_t typesize = InputOutput::Variables::typeToBytes(
        data->field->type());
    void *ptr = (void*)((char*)data->ptr + offset * typesize);
    MPI::COMM_WORLD.Recv(ptr, n * mpi_t.n, mpi_t.t, data->proc, data->tag);

    // We are set, time to upload the data to the 
    err_code = clEnqueueWriteBuffer(data->C->command_queue(),
                                    *(cl_mem*)data->field->get(),
                                    CL_TRUE,
                                    offset * typesize,
                                    n * typesize,
                                    ptr,
                                    0,
                                    NULL,
                                    NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure uploading the variable \"" << data->field->name()
            << "\" received from MPI." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // This is not the ideal situation, but the safest way to avoid racing
    // conditions is to synchronouly wait until the mask has been completely set
    if(mask_event != NULL) {
        err_code = clWaitForEvents(1, &mask_event);
        if(err_code != CL_SUCCESS) {
            LOG(L_ERROR, "Failure waiting mask to be set");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
    }
    
    // Unlock the resources
    err_code = clSetUserEventStatus(data->field->getEvent(), CL_COMPLETE);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure setting user event status for \"" <<
               data->field->name() << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void MPISync::Receiver::execute()
{
    unsigned int i;
    cl_int err_code;
    cl_event event;
    CalcServer *C = CalcServer::singleton();

    for(i = 0; i <= _fields.size(); i++) {
        // Nothing can be done until we have free way
        const cl_event event_wait_list[3] = {
            _n_offset->getEvent(),
            _fields.at(i)->getEvent(),
            _mask->getEvent(),            
        };
        err_code = clEnqueueMarkerWithWaitList(C->command_queue(),
                                               3,
                                               event_wait_list,
                                               &event);
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure creating recv events syncing point for field \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }

        // As it happened with the sender, we are setting up a user event that we
        // can use to lock the simultaneous access to the fields and offset
        // variables
        cl_event field_event = clCreateUserEvent(C->context(), &err_code);
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure creating recv user event for \"" <<
                _fields.at(i)->name() << "\" variable in tool \"" << name()
                << "\"" << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        _fields.at(i)->setEvent(field_event);
        if(i == 0) {
            _n_offset->setEvent(field_event);
            _mask->setEvent(field_event);
        }

        // Asynchronously receive and upload the data
        MPISyncRecvUserData user_data;
        user_data.C = C;
        user_data.field = _fields.at(i);
        user_data.recv_n = i == 0;
        user_data.proc = _proc;
        user_data.offset = _n_offset;
        user_data.n = &_n_recv;
        user_data.ptr = _fields_host.at(i);
        user_data.tag = i + 1;
        user_data.mask = _mask;
        user_data.kernel = _kernel;
        user_data.local_work_size = _local_work_size;
        // So we can asynchronously ask to dispatch the data download and send
        // when reduction is done.
        err_code = clSetEventCallback(event,
                                      CL_COMPLETE,
                                      &cbReceiveAndUpload,
                                      (void*)(&user_data));
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure setting the receive & upload callback for \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        err_code = clReleaseEvent(field_event);
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure releasing user event for \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
    }
}

void MPISync::Receiver::setupOpenCL()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream source;
    source << MPISYNC_INC << MPISYNC_SRC;
    _kernel = compile_kernel(source.str(), "set_mask");

    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(_local_work_size < __CL_MIN_LOCALSIZE__){
        std::stringstream msg;
        LOG(L_ERROR, "UnSort cannot be performed.\n");
        msg << "\t" << _local_work_size
            << " elements can be executed, but __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
    }

    err_code = clSetKernelArg(_kernel,
                              0,
                              _mask->typesize(),
                              _mask->get());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the mask argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_kernel,
                              1,
                              sizeof(unsigned int),
                              (void*)&_proc);
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the proc argument\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
}



}}  // namespaces

#endif
