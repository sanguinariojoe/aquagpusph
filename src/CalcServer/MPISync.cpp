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
#include <mpi.h>

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
    for(auto sender : _senders) {
        delete sender;
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
}

cl_event MPISync::_execute(const std::vector<cl_event> events)
{
    if(!_procs.size())
        return NULL;

    // Sort the mask
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
    for(auto sender : _senders) {
        sender->execute();
    }

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

    // Create the senders
    for(auto proc : _procs) {
        Sender* sender = new Sender(name(),
                                    _mask,
                                    _fields_sorted,
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

MPISync::Exchanger::Exchanger(const std::string tool_name,
                              InputOutput::ArrayVariable *mask,
                              const std::vector<InputOutput::ArrayVariable*> fields,
                              const unsigned int proc)
    : _name(tool_name)
    , _mask(mask)
    , _fields(fields)
    , _proc(proc)
    , _n(0)
{
    _n = _mask->size() / InputOutput::Variables::typeToBytes(_mask->type());

    for(auto field : fields) {
        void *data = malloc(field->size());
        if(!data) {
            std::stringstream msg;
            msg << "Failure allocating \"" << field->size()
                << "\" bytes for the array \"" << field->name() 
                << "\" in the tool \"" << name() << "\"" << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::bad_alloc();
        }
        _fields_host.push_back(data);
    }
}

MPISync::Exchanger::~Exchanger()
{
    for(auto field : _fields_host) {
        free(field);
    }
}

MPISync::Sender::Sender(const std::string name,
                        InputOutput::ArrayVariable *mask,
                        const std::vector<InputOutput::ArrayVariable*> fields,
                        const unsigned int proc,
                        InputOutput::UIntVariable *n_offset)
    : MPISync::Exchanger(name, mask, fields, proc)
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
} MPISyncDownloadUserData;

void CL_CALLBACK cbDownloadAndSend(cl_event n_event,
                                   cl_int cmd_exec_status,
                                   void *user_data)
{
    cl_int err_code;
    MPISyncDownloadUserData *data = (MPISyncDownloadUserData*)user_data;
    unsigned int n = *(unsigned int*)data->n->get();
    const unsigned int offset = *(unsigned int*)data->offset->get();

    if(data->send_n) {
        // We need to report the proc the number of particles to be sent. We
        // don't want to hang here, so we aqsk for async sending
        MPI::COMM_WORLD.Isend(&n, 1, MPI::UNSIGNED, data->proc, 0);
        // We know this is the very only sending instance (to this proc)
        // executing this piece of code, so we can take advantage to increase
        // the offset variable
        const unsigned int next_offset = offset + n;
        data->offset->set((void*)(&next_offset));
    }

    if(!n) {
        // There is nothing to send
        return;
    }

    // We can now proceed to download the data. We are in a parallel thread, so
    // we can do this synchronously
    size_t typesize = InputOutput::Variables::typeToBytes(data->field->type());
    err_code = clEnqueueReadBuffer(data->C->command_queue(),
                                   *(cl_mem*)data->field->get(),
                                   CL_TRUE,
                                   offset * typesize,
                                   n * typesize,
                                   data->ptr,
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure downloading the variable \"" << data->field->name()
            << "\" to send that by MPI." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        return;
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
    std::string t = trimCopy(data->field->type());
    MPI_Datatype datatype;

    if((t.back() == '*')){
        t.pop_back();
    }
    if(hasSuffix(t, "vec")) {
#ifdef HAVE_3D
        n *= 4;
#else
        n *= 2;
#endif
    } else if(hasSuffix(t, "2")) {
        n *= 2;
        t.pop_back();
    } else if(hasSuffix(t, "3")) {
        n *= 3;
        t.pop_back();
    } else if(hasSuffix(t, "4")) {
        n *= 4;
        t.pop_back();
    }

    if((!t.compare("int")) || (!t.compare("ivec"))) {
        datatype = MPI::INT;
    } else if((!t.compare("unsigned int")) || (!t.compare("uivec"))) {
        datatype = MPI::UNSIGNED;
    } else if((!t.compare("float")) || (!t.compare("vec"))) {
        datatype = MPI::FLOAT;
    } else {
        std::ostringstream msg;
        msg << "Unrecognized type \"" << data->field->type() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        return;
    }

    // Launch the missiles. Again we can proceed in synchronous mode
    MPI::COMM_WORLD.Send(data->ptr, n, datatype, data->proc, data->tag);
}

void MPISync::Sender::execute()
{
    unsigned int i;
    cl_int err_code;
    cl_event event;
    CalcServer *C = CalcServer::singleton();

    // Extract the submask
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
    _n_send_reduction->execute();

    for(i = 0; i <= _fields.size(); i++) {
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
        
        cl_event field_event = clCreateUserEvent(C->context(), &err_code);
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure creating user event for \"" <<
                _fields.at(i)->name() << "\" variable in tool \"" << name()
                << "\"" << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
        _fields.at(i)->setEvent(field_event);
        if(i == 0)
            _n_offset->setEvent(field_event);

        MPISyncDownloadUserData user_data;
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
        err_code = clEnqueueMarkerWithWaitList(C->command_queue(),
                                               3,
                                               event_wait_list,
                                               &event);
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure creating events syncing point for field \""
                << _fields.at(i)->name() << "\" in tool \"" << name() << "\""
                << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
        }
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
    }
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
    _kernel = compile(source.str());

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

cl_kernel MPISync::Sender::compile(const std::string source)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream flags;
    #ifdef AQUA_DEBUG
        flags << " -DDEBUG ";
    #else
        flags << " -DNDEBUG ";
    #endif
    flags << " -cl-mad-enable -cl-fast-relaxed-math";
    #ifdef HAVE_3D
        flags << " -DHAVE_3D";
    #else
        flags << " -DHAVE_2D";
    #endif
    size_t source_length = source.size();
    const char* source_cstr = source.c_str();
    program = clCreateProgramWithSource(C->context(),
                                        1,
                                        &source_cstr,
                                        &source_length,
                                        &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the OpenCL program\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clBuildProgram(program, 0, NULL, flags.str().c_str(), NULL, NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Error compiling the OpenCL script\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        LOG0(L_ERROR, "--- Build log ---------------------------------\n");
        size_t log_size = 0;
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              0,
                              NULL,
                              &log_size);
        char *log = (char*)malloc(log_size + sizeof(char));
        if(!log){
            std::stringstream msg;
            msg << "Failure allocating " << log_size
                << " bytes for the building log" << std::endl;
            LOG0(L_ERROR, msg.str());
            LOG0(L_ERROR, "--------------------------------- Build log ---\n");
            throw std::bad_alloc();
        }
        strcpy(log, "");
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              log_size,
                              log,
                              NULL);
        strcat(log, "\n");
        LOG0(L_DEBUG, log);
        LOG0(L_ERROR, "--------------------------------- Build log ---\n");
        free(log); log=NULL;
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL compilation error");
    }
    kernel = clCreateKernel(program, "extract_submask", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the OpenCL kernel\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    return kernel;
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
}

}}  // namespaces

#endif
