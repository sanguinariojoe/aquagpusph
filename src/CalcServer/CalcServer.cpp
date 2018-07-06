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
 * @brief The calculation main entry point.
 * (See Aqua::CalcServer::CalcServer for details)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <sys/time.h>
#include <unistd.h>

#include <CalcServer.h>
#include <AuxiliarMethods.h>
#include <FileManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <ScreenManager.h>
#include <CalcServer/Tool.h>
#include <CalcServer/Assert.h>
#include <CalcServer/Copy.h>
#include <CalcServer/Kernel.h>
#include <CalcServer/LinkList.h>
#include <CalcServer/Python.h>
#include <CalcServer/RadixSort.h>
#include <CalcServer/Reduction.h>
#include <CalcServer/Set.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer/UnSort.h>
#include <CalcServer/Reports/Performance.h>
#include <CalcServer/Reports/Screen.h>
#include <CalcServer/Reports/TabFile.h>
#include <CalcServer/Reports/SetTabFile.h>

namespace Aqua{ namespace CalcServer{

CalcServer::CalcServer(const Aqua::InputOutput::ProblemSetup& sim_data)
    : _num_platforms(0)
    , _platforms(NULL)
    , _num_devices(0)
    , _devices(NULL)
    , _context(NULL)
    , _command_queues(NULL)
    , _platform(NULL)
    , _device(NULL)
    , _command_queue(NULL)
    , _current_tool_name(NULL)
    , _sim_data(sim_data)
{
    unsigned int i, j;

    if(setupOpenCL()) {
        exit(EXIT_FAILURE);
    }

    _base_path = _sim_data.settings.base_path;
    _current_tool_name = new char[256];
    strcpy(_current_tool_name, "");

    unsigned int num_sets = _sim_data.sets.size();
    unsigned int N = 0;
    for(i = 0; i < _sim_data.sets.size(); i++) {
        N += _sim_data.sets.at(i)->n();
    }

    unsigned int num_icell = nextPowerOf2(N);
    num_icell = roundUp(num_icell, _ITEMS*_GROUPS);

    // Register default scalars
    std::ostringstream valstr;
    #ifdef HAVE_3D
        unsigned int dims = 3;
    #else
        unsigned int dims = 2;
    #endif
    valstr.str(""); valstr << dims;
    _vars.registerVariable("dims", "unsigned int", "", valstr.str());
    _vars.registerVariable("t", "float", "", "0");
    _vars.registerVariable("dt", "float", "", "0");
    _vars.registerVariable("iter", "unsigned int", "", "0");
    _vars.registerVariable("frame", "unsigned int", "", "0");
    valstr.str(""); valstr << std::numeric_limits<float>::max();
    _vars.registerVariable("end_t", "float", "", valstr.str());
    valstr.str(""); valstr << std::numeric_limits<unsigned int>::max();
    _vars.registerVariable("end_iter", "unsigned int", "", valstr.str());
    _vars.registerVariable("end_frame", "unsigned int", "", valstr.str());
    valstr.str(""); valstr << N;
    _vars.registerVariable("N", "unsigned int", "", valstr.str());
    valstr.str(""); valstr << _sim_data.sets.size();
    _vars.registerVariable("n_sets", "unsigned int", "", valstr.str());
    valstr.str(""); valstr << num_icell;
    _vars.registerVariable("n_radix", "unsigned int", "", valstr.str());
    // Number of cells in x, y, z directions, and the total (n_x * n_y * n_z)
    _vars.registerVariable("n_cells", "uivec4", "", "0, 0, 0, 0");
    // Kernel support
    _vars.registerVariable("support", "float", "", "2");

    // Register default arrays
    valstr.str(""); valstr << N;
    _vars.registerVariable("id", "unsigned int*", valstr.str(), "");
    _vars.registerVariable("r", "vec*", valstr.str(), "");
    _vars.registerVariable("iset", "unsigned int*", valstr.str(), "");
    valstr.str(""); valstr << num_icell;
    _vars.registerVariable("id_sorted", "unsigned int*", valstr.str(), "");
    _vars.registerVariable("id_unsorted", "unsigned int*", valstr.str(), "");
    _vars.registerVariable("icell", "unsigned int*", valstr.str(), "");
    _vars.registerVariable("ihoc", "unsigned int*", "n_cells_w", "");

    // Register the user variables and arrays
    for(i = 0; i < _sim_data.variables.names.size(); i++){
        _vars.registerVariable(_sim_data.variables.names.at(i),
                                _sim_data.variables.types.at(i),
                                _sim_data.variables.lengths.at(i),
                                _sim_data.variables.values.at(i));
    }

    // Register the user definitions
    for(i = 0; i < _sim_data.definitions.names.size(); i++){
        valstr.str("");
        valstr << "-D" << _sim_data.definitions.names.at(i);
        if(!_sim_data.definitions.values.at(i).compare("")){
            if(!_sim_data.definitions.evaluations.at(i)) {
                valstr << "=" << _sim_data.definitions.values.at(i);
            }
            else {
                float defval = 0.f;
                _vars.solve("float",
                            _sim_data.definitions.values.at(i),
                            &defval);
                // We need to specify the format, to ensure that a decimal
                // number (i.e. a number with decimal point) is retrieved, so
                // C++ streamer can't be directly applied
                char defvalstr[100];
                snprintf(defvalstr, sizeof(defvalstr), "%#G", defval);
                valstr << "=" << defvalstr << "f";
            }
        }
        _definitions.push_back(valstr.str());
    }

    // Register the tools
    std::deque<Aqua::InputOutput::ProblemSetup::sphTool*>::iterator t_it;
    for(t_it = _sim_data.tools.begin(); t_it < _sim_data.tools.end(); t_it++){
        if(!(*t_it)->get("type").compare("kernel")){
            Kernel *tool = new Kernel((*t_it)->get("name"),
                                      (*t_it)->get("path"),
                                      (*t_it)->get("entry_point"),
                                      (*t_it)->get("n"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("copy")){
            Copy *tool = new Copy((*t_it)->get("name"),
                                  (*t_it)->get("in"),
                                  (*t_it)->get("out"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("python")){
            Python *tool = new Python((*t_it)->get("name"),
                                      (*t_it)->get("path"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("set")){
            Set *tool = new Set((*t_it)->get("name"),
                                (*t_it)->get("in"),
                                (*t_it)->get("value"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("set_scalar")){
            SetScalar *tool = new SetScalar((*t_it)->get("name"),
                                            (*t_it)->get("in"),
                                            (*t_it)->get("value"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("reduction")){
            Reduction *tool = new Reduction((*t_it)->get("name"),
                                            (*t_it)->get("in"),
                                            (*t_it)->get("out"),
                                            (*t_it)->get("operation"),
                                            (*t_it)->get("null"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("link-list")){
            LinkList *tool = new LinkList((*t_it)->get("name"),
                                          (*t_it)->get("in"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("radix-sort")){
            RadixSort *tool = new RadixSort((*t_it)->get("name"),
                                            (*t_it)->get("in"),
                                            (*t_it)->get("perm"),
                                            (*t_it)->get("inv_perm"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("assert")){
            Assert *tool = new Assert((*t_it)->get("name"),
                                      (*t_it)->get("condition"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("dummy")){
            Tool *tool = new Tool((*t_it)->get("name"));
            _tools.push_back(tool);
        }
        else{
            std::ostringstream msg;
            msg << "Unrecognized tool type \""
                << (*t_it)->get("type")
                << "\" when parsing the tool \""
                << (*t_it)->get("name") << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid tool type");
        }
    }

    // Register the reporters
    for(t_it = _sim_data.reports.begin(); t_it < _sim_data.reports.end(); t_it++){
        if(!(*t_it)->get("type").compare("screen")){
            bool bold = false;
            if(!(*t_it)->get("bold").compare("true") ||
               !(*t_it)->get("bold").compare("True")){
               bold = true;
            }
            Reports::Screen *tool = new Reports::Screen(
                (*t_it)->get("name"),
                (*t_it)->get("fields"),
                (*t_it)->get("color"),
                bold);
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("file")){
            Reports::TabFile *tool = new Reports::TabFile(
                (*t_it)->get("name"),
                (*t_it)->get("fields"),
                (*t_it)->get("path"));
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("particles")){
            // Get the first particle associated to this set
            unsigned int set_id = std::stoi((*t_it)->get("set"));
            if(set_id >= _sim_data.sets.size()){
                std::ostringstream msg;
                msg << "Report \"" << (*t_it)->get("name")
                    << "\" requested the particles set " << set_id
                    << " but just " << _sim_data.sets.size()
                    << " can be found." << std::endl;
                LOG(L_ERROR, msg.str());
                throw std::runtime_error("particles set out of bounds");
            }
            unsigned int first = 0;
            for(j = 0; j < set_id; j++){
                first += _sim_data.sets.at(j)->n();
            }

            // And the ipf and fps
            unsigned int ipf = std::stoi((*t_it)->get("ipf"));
            float fps = std::stof((*t_it)->get("fps"));

            Reports::SetTabFile *tool = new Reports::SetTabFile(
                (*t_it)->get("name"),
                (*t_it)->get("fields"),
                first,
                _sim_data.sets.at(set_id)->n(),
                (*t_it)->get("path"),
                ipf,
                fps);
            _tools.push_back(tool);
        }
        else if(!(*t_it)->get("type").compare("performance")){
            bool bold = false;
            if(!(*t_it)->get("bold").compare("true") ||
               !(*t_it)->get("bold").compare("True")){
               bold = true;
            }
            Reports::Performance *tool = new Reports::Performance(
                (*t_it)->get("name"),
                (*t_it)->get("color"),
                bold,
                (*t_it)->get("path"));
            _tools.push_back(tool);
        }
        else{
            std::ostringstream msg;
            msg << "Unrecognized report type \"" << (*t_it)->get("type")
                << "\" when parsing the report \"" << (*t_it)->get("name")
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid report type");
        }
    }

    std::ostringstream msg;
    msg << "Allocated memory = " << _vars.allocatedMemory()
        << " bytes" << std::endl;
    LOG(L_INFO, msg.str());

    setup();
}

CalcServer::~CalcServer()
{
    unsigned int i;
    delete[] _current_tool_name;

    if(_context) clReleaseContext(_context); _context = NULL;
    for(i = 0; i < _num_devices; i++){
        if(_command_queues[i]) clReleaseCommandQueue(_command_queues[i]);
        _command_queues[i] = NULL;
    }

    if(_platforms) delete[] _platforms; _platforms=NULL;
    if(_devices) delete[] _devices; _devices=NULL;
    if(_command_queues) delete[] _command_queues; _command_queues=NULL;

    for(i = 0; i < _tools.size(); i++){
        delete _tools.at(i);
        _tools.at(i) = NULL;
    }
    _tools.clear();

    std::map<std::string, UnSort*>::iterator it = unsorters.begin();
    while(it != unsorters.end())
    {
        delete it->second;
        it++;
    }
}

void CalcServer::update()
{
    InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
    unsigned int i;
    while(!T->mustPrintOutput() && !T->mustStop()){
        InputOutput::ScreenManager::singleton()->initFrame();

        // Execute the tools
        strcpy(_current_tool_name, "__pre execution__");
        for(i = 0; i < _tools.size(); i++){
            strncpy(_current_tool_name, _tools.at(i)->name().c_str(), 255);
            _current_tool_name[255] = '\0';
            try {
                _tools.at(i)->execute();
            } catch (std::runtime_error &e) {
                sleep(__ERROR_SHOW_TIME__);
                throw;
            }
        }
        strcpy(_current_tool_name, "__post execution__");

        // Key events
        while(isKeyPressed()){
            if(getchar() == 'c'){
                LOG(L_WARNING, "Interrumption request detected.\n");
                sleep(__ERROR_SHOW_TIME__);
                throw std::runtime_error("Simulation interrupted by the user");
            }
        }

        clFinish(command_queue());

        InputOutput::ScreenManager::singleton()->endFrame();
    }
}

cl_event CalcServer::getUnsortedMem(const std::string var_name,
                                    size_t offset,
                                    size_t cb,
                                    void *ptr)
{
    cl_int err_code;

    // Generate the unsorted if it does not exist yet
    UnSort *unsorter = NULL;
    if(unsorters.find(var_name) == unsorters.end()){
        unsorter = new UnSort(var_name.c_str(), var_name.c_str());
        try {
            unsorter->setup();
        } catch(std::runtime_error &e) {
            delete unsorter;
            return NULL;
        }
        unsorters.insert(std::make_pair(var_name, unsorter));
    }
    // Get the unsorter
    unsorter = unsorters[var_name];
    try {
        unsorter->execute();
    } catch (std::runtime_error &e) {
        return NULL;
    }
    cl_mem mem = unsorter->output();
    cl_event event = NULL;
    err_code = clEnqueueReadBuffer(command_queue(),
                                   mem,
                                   CL_FALSE,
                                   offset,
                                   cb,
                                   ptr,
                                   0,
                                   NULL,
                                   &event);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure receiving the variable \"" << var_name
            << "\" from server." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        return NULL;
    }
    return event;
}

bool CalcServer::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    LOG(L_INFO, "Initializating OpenCL...\n");
    if(queryOpenCL()){
        return true;
    }
    if(setupPlatform()){
        return true;
    }
    if(setupDevices()){
        return true;
    }
    LOG(L_INFO, "OpenCL is ready to work!\n");
    return false;
}

bool CalcServer::queryOpenCL()
{
    cl_int err_code;
    cl_uint i, j, num_devices=0;
    cl_device_id *devices;
    char msg[1024], aux[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    _platforms = NULL;
    strcpy(msg, "");
    // Gets the total number of platforms
    err_code = clGetPlatformIDs(0, NULL, &_num_platforms);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting the number of platforms.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Get the array of platforms
    _platforms = new cl_platform_id[_num_platforms];
    if(!_platforms) {
        LOG(L_ERROR, "Allocation memory error.\n");
        LOG0(L_DEBUG, "\tPlatforms array cannot be allocated\n");
        return true;
    }
    err_code = clGetPlatformIDs(_num_platforms, _platforms, NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting the platforms list.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    for(i = 0; i < _num_platforms; i++){
        // Get the number of devices
        err_code = clGetDeviceIDs(_platforms[i],
                                  CL_DEVICE_TYPE_ALL,
                                  0,
                                  NULL,
                                  &num_devices);
        if(err_code != CL_SUCCESS) {
            LOG(L_ERROR, "Failure getting the number of devices.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        // Gets the devices array
        devices = new cl_device_id[num_devices];
        if(!devices) {
            LOG(L_ERROR, "Allocation memory error.\n");
            LOG0(L_DEBUG, "\tDevices array cannot be allocated\n");
            return true;
        }
        err_code = clGetDeviceIDs(_platforms[i], CL_DEVICE_TYPE_ALL,
                                  num_devices, devices, NULL);
        if(err_code != CL_SUCCESS) {
            LOG(L_ERROR, "Failure getting the devices list.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        // Shows device arrays
        for(j = 0; j < num_devices; j++){
            // Identifier
            strcpy(msg, "");
            sprintf(msg, "\tDevice %u, Platform %u...\n", j, i);
            LOG0(L_INFO, msg);
            // Device name
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_NAME,
                                       1024 * sizeof(char),
                                       &aux,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                LOG(L_ERROR, "Failure getting the device name.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            strcpy(msg, "");
            sprintf(msg, "\t\tDEVICE: %s\n", aux);
            LOG0(L_DEBUG, msg);
            // Platform vendor
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_VENDOR,
                                       1024 * sizeof(char),
                                       &aux,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                LOG(L_ERROR, "Failure getting the device vendor.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            strcpy(msg, "");
            sprintf(msg, "\t\tVENDOR: %s\n", aux);
            LOG0(L_DEBUG, msg);
            // Device type
            cl_device_type dType;
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_TYPE,
                                       sizeof(cl_device_type),
                                       &dType,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                LOG(L_ERROR, "Failure getting the device type.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            if(dType == CL_DEVICE_TYPE_CPU)
                LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_CPU\n");
            else if(dType == CL_DEVICE_TYPE_GPU)
                LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_GPU\n");
            else if(dType == CL_DEVICE_TYPE_ACCELERATOR)
                LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_ACCELERATOR\n");
            else if(dType == CL_DEVICE_TYPE_DEFAULT)
                LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_DEFAULT\n");
            #ifdef CL_DEVICE_TYPE_CUSTOM
            else if(dType == CL_DEVICE_TYPE_CUSTOM)
                LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_CUSTOM\n");
            #endif
            else {
                sprintf(msg, "\t\tTYPE: %ul\n", dType);
                LOG0(L_DEBUG, msg);
            }
        }
        delete[] devices; devices = NULL;
    }
    return false;
}

bool CalcServer::setupPlatform()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    if(_sim_data.settings.platform_id >= _num_platforms){
        LOG(L_ERROR, "Impossible to use the requested platform.\n");
        strcpy(msg, "");
        sprintf(msg, "\t%u platform has been selected, but just %u platforms have been found\n",
                _sim_data.settings.platform_id, _num_platforms);
        LOG0(L_DEBUG, msg);
        return true;
    }
    _platform = _platforms[_sim_data.settings.platform_id];
    return false;
}

/** @brief Runtime error reporting tool
 *
 * Errors reported in this way directly depends on the implementation.
 * @param errinfo is a pointer to an error string.
 * @param private_info pointer to binary data that is returned by the OpenCL
 * implementation that can be used to log additional information helpful in
 * debugging the error.
 * @param cb Size of the binary data, #private_info.
 * @param user_data Current tool name.
 */
void CL_CALLBACK context_error_notify(const char *errinfo,
                                      const void *private_info,
                                      size_t cb,
                                      void *user_data)
{
    const char* current_tool_name = (const char*)user_data;
    char msg[512];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    strcpy(msg, "OpenCL implementation reported a runtime error");
    if(strlen(current_tool_name)){
        strcat(msg, " (");
        strcat(msg, current_tool_name);
        strcat(msg, ")");
    }
    strcat(msg, ":\n");

    LOG(L_ERROR, msg);
    LOG0(L_DEBUG, errinfo);
    LOG0(L_DEBUG, "\n");
} 

bool CalcServer::setupDevices()
{
    cl_int err_code;
    cl_uint i;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    _devices = NULL;
    // Gets the number of valid devices
    err_code = clGetDeviceIDs(_platform,
                              _sim_data.settings.device_type,
                              0,
                              NULL,
                              &_num_devices);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting the number of devices.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(_sim_data.settings.device_id >= _num_devices) {
        LOG(L_ERROR, "Impossible to use the selected device.\n");
        strcpy(msg, "");
        sprintf(msg, "\t%u device has been selected, but just %u devices are available\n",
                _sim_data.settings.device_id, _num_devices);
        LOG0(L_DEBUG, msg);
        if(_sim_data.settings.device_type == CL_DEVICE_TYPE_ALL)
            LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_ALL filter activated\n");
        else if(_sim_data.settings.device_type == CL_DEVICE_TYPE_CPU)
            LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_CPU filter activated\n");
        else if(_sim_data.settings.device_type == CL_DEVICE_TYPE_GPU)
            LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_GPU filter activated\n");
        else if(_sim_data.settings.device_type == CL_DEVICE_TYPE_ACCELERATOR)
            LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_ACCELERATOR filter activated\n");
        else if(_sim_data.settings.device_type == CL_DEVICE_TYPE_DEFAULT)
            LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_DEFAULT filter activated\n");

        return true;
    }
    // Gets the devices array
    _devices = new cl_device_id[_num_devices];
    if(!_devices){
        LOG(L_ERROR, "Allocation memory error.\n");
        LOG0(L_DEBUG, "\tDevices array can't be allocated\n");
        return true;
    }
    err_code = clGetDeviceIDs(_platform,
                              _sim_data.settings.device_type,
                              _num_devices,
                              _devices,
                              &_num_devices);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting the devices list.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Create a devices context
    _context = clCreateContext(0,
                               _num_devices,
                               _devices,
                               context_error_notify,
                               _current_tool_name,
                               &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating an OpenCL context.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Create command queues
    _command_queues = new cl_command_queue[_num_devices];
    if(_command_queues == NULL){
        LOG(L_ERROR, "Allocation memory error.\n");
        LOG0(L_DEBUG, "\tCommand queues array cannot be allocated\n");
        return true;
    }
    for(i = 0; i < _num_devices; i++) {
        _command_queues[i] = clCreateCommandQueue(_context,
                                                  _devices[i],
                                                  0,
                                                  &err_code);
        if(err_code != CL_SUCCESS) {
            strcpy(msg, "");
            sprintf(msg, "Can't create a command queue for the device %u.\n",i);
            LOG(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
    }
    // Store the selected ones
    _device = _devices[_sim_data.settings.device_id];
    _command_queue = _command_queues[_sim_data.settings.device_id];
    return false;
}

void CalcServer::setup()
{
    unsigned int i, j;
    cl_uint err_code=0;

    // Check for the required variables that must be defined by the user
    if(!_vars.get("h")){
        LOG(L_ERROR, "Missing kernel length \"h\".\n");
        throw std::runtime_error("Undeclared kernel length variable");
    }
    if(_vars.get("h")->type().compare("float")){
        std::ostringstream msg;
        msg << "Kernel length \"h\" must be of type \"float\", not \""
            << _vars.get("h")->type() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid kernel length variable type");
    }
    float h = *(float *)_vars.get("h")->get();
    if(h <= 0.f){
        std::ostringstream msg;
        msg << "Kernel length \"h\" must be a positive value (h="
            << h << ")" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid kernel length value");
    }

    std::ostringstream msg;
    msg << "Found kernel length, h = " << h << " [m]" << std::endl;
    LOG(L_INFO, msg.str());

    // Setup the scalar data per particles sets
    for(i = 0; i < _sim_data.sets.size(); i++){
        InputOutput::ProblemSetup::sphParticlesSet* set = _sim_data.sets.at(i);
        for(j = 0; j < set->scalarNames().size(); j++){
            std::string name = set->scalarNames().at(j);
            std::string val = set->scalarValues().at(j);
            if(!_vars.get(name)){
                std::ostringstream msg;
                msg << "Particles set " << i
                    << " asks for the undeclared variable \"" << name
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                throw std::runtime_error("Invalid variable");
            }
            if(_vars.get(name)->type().find('*') == std::string::npos){
                std::ostringstream msg;
                msg << "Particles set " << i
                    << " can't set the value of a scalar variable, \"" << name
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                throw std::runtime_error("Invalid variable type");
            }
            InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable *)_vars.get(name);
            size_t typesize = _vars.typeToBytes(_vars.get(name)->type());
            size_t len = _vars.get(name)->size() / typesize;
            if(len != _sim_data.sets.size()){
                std::ostringstream msg;
                msg << "Particles set " << i
                    << " can't set the value of the array variable, \"" << name
                    << "\", because its length, " << len
                    << "differs from the number of sets, "
                    << _sim_data.sets.size() << std::endl;
                LOG(L_ERROR, msg.str());
                throw std::runtime_error("Invalid variable length");
            }
            void *data = malloc(typesize);
            try {
                _vars.solve(_vars.get(name)->type(), val, data);
            } catch(...) {
                std::ostringstream msg;
                msg << "Particles set " << i
                    << " failed evaluating variable, \"" << name
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                throw;
            }
            cl_mem mem = *(cl_mem*)_vars.get(name.c_str())->get();
            err_code = clEnqueueWriteBuffer(_command_queue, mem, CL_TRUE,
                                            i * typesize, typesize, data,
                                            0, NULL, NULL);
            free(data); data = NULL;
            if(err_code != CL_SUCCESS) {
                std::ostringstream msg;
                msg << "Particles set " << i
                    << " failed sending variable, \"" << name
                    << "\" to the computational device." << std::endl;
                LOG(L_ERROR, msg.str());
                InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
                throw std::runtime_error("OpenCL error");
            }
        }
    }

    // Setup the tools
    for(i = 0; i < _tools.size(); i++){
        _tools.at(i)->setup();
    }
}

}}  // namespace
