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

#include <CalcServer.h>
#include <AuxiliarMethods.h>
#include <FileManager.h>
#include <ProblemSetup.h>
#include <Fluid.h>
#include <TimeManager.h>
#include <ScreenManager.h>
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

CalcServer::CalcServer()
    : _num_platforms(0)
    , _platforms(NULL)
    , _num_devices(0)
    , _devices(NULL)
    , _context(NULL)
    , _command_queues(NULL)
    , _platform(NULL)
    , _device(NULL)
    , _command_queue(NULL)
    , _vars(NULL)
{
    unsigned int i;
    char msg[1024];
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if(setupOpenCL()) {
        exit(EXIT_FAILURE);
    }

    unsigned int num_sets = P->sets.size();
    unsigned int num_sensors = P->SensorsParameters.pos.size();
    unsigned int n = 0;
    for(i = 0; i < P->sets.size(); i++) {
        n += P->sets.at(i)->n();
    }
    unsigned int N = n + num_sensors;

    unsigned int num_icell = nextPowerOf2(N);
    num_icell = roundUp(num_icell, _ITEMS*_GROUPS);

    _vars = new InputOutput::Variables();

    // Register default scalars
    char val[16];
    char len[16];
    strcpy(len, "");
    sprintf(val, "%g", 0.f);
    if(_vars->registerVariable("t", "float", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%g", 0.f);
    if(_vars->registerVariable("dt", "float", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", 0);
    if(_vars->registerVariable("step", "unsigned int", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", n);
    if(_vars->registerVariable("n", "unsigned int", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", num_sensors);
    if(_vars->registerVariable("n_sensors", "unsigned int", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", N);
    if(_vars->registerVariable("N", "unsigned int", len, val, true))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", P->sets.size());
    if(_vars->registerVariable("n_sets", "unsigned int", len, val, false))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", num_icell);
    if(_vars->registerVariable("n_radix", "unsigned int", len, val, false))
        exit(EXIT_FAILURE);
    // Number of cells in x, y, z directions, and the total (n_x * n_y * n_z)
    strcpy(val, "0, 0, 0, 0");
    if(_vars->registerVariable("n_cells", "uivec4", len, val, false))
        exit(EXIT_FAILURE);
    // Register default arrays
    strcpy(val, "");
    sprintf(len, "%u", N);
    if(_vars->registerVariable("pos", "vec*", len, val, true))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("iset", "unsigned int*", len, val, true))
        exit(EXIT_FAILURE);
    sprintf(len, "%u", num_icell);
    if(_vars->registerVariable("id_sorted", "unsigned int*", len, val, true))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("id_unsorted", "unsigned int*", len, val, true))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("icell", "unsigned int*", len, val, true))
        exit(EXIT_FAILURE);
    sprintf(len, "n_cells_w");
    if(_vars->registerVariable("ihoc", "unsigned int*", len, val, true))
        exit(EXIT_FAILURE);

    // Register the user variables and arrays
    for(i = 0; i < P->variables.names.size(); i++){
        bool flag = _vars->registerVariable(P->variables.names.at(i),
                                            P->variables.types.at(i),
                                            P->variables.lengths.at(i),
                                            P->variables.values.at(i),
                                            P->variables.saves.at(i));
        if(flag){
            exit(EXIT_FAILURE);
        }
    }

    // Register the tools
    for(i = 0; i < P->tools.size(); i++){
        if(!strcmp(P->tools.at(i)->get("type"), "kernel")){
            Kernel *tool = new Kernel(P->tools.at(i)->get("name"),
                                      P->tools.at(i)->get("path"));
            _tools.push_back(tool);
        }
        else{
            sprintf(msg,
                    "Unrecognized tool type \"%s\" when parsing the tool \"%s\".\n",
                    P->tools.at(i)->get("type"),
                    P->tools.at(i)->get("name"));
            S->addMessageF(3, msg);
            exit(EXIT_FAILURE);
        }
    }

    sprintf(msg, "Allocated memory = %lu bytes\n", _vars->allocatedMemory());
    S->addMessageF(1, msg);
}

CalcServer::~CalcServer()
{
    unsigned int i;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    S->addMessageF(1, "Sutting down the OpenCL context...\n");
    if(_context) clReleaseContext(_context); _context = NULL;
    for(i = 0; i < _num_devices; i++){
        if(_command_queues[i]) clReleaseCommandQueue(_command_queues[i]);
        _command_queues[i] = NULL;
    }

    S->addMessageF(1, "Deallocating host memory...\n");
    if(_platforms) delete[] _platforms; _platforms=NULL;
    if(_devices) delete[] _devices; _devices=NULL;
    if(_command_queues) delete[] _command_queues; _command_queues=NULL;

    S->addMessageF(1, "Destroying tools...\n");
    for(i = 0; i < _tools.size(); i++){
        delete _tools.at(i);
        _tools.at(i) = NULL;
    }
    _tools.clear();

    S->addMessageF(1, "Destroying variables manager...\n");
    if(_vars) delete _vars; _vars=NULL;
}

bool CalcServer::update()
{
    InputOutput::TimeManager *T   = InputOutput::TimeManager::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    unsigned int i;
    while(!T->mustPrintOutput() && !T->mustStop()){
        // Execute the tools
        for(i = 0; i < _tools.size(); i++){
            if(_tools.at(i)->execute()){
                return true;
            }
        }
        // Key events
        while(isKeyPressed()){
            if(getchar() == 'c'){
                S->addMessageF(1, "Interrumption request detected.\n");
                return true;
            }
        }

        /// @todo let the tool to continue computing
        return true;
    }
    return false;
}

bool CalcServer::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    S->addMessageF(1, "Initializating OpenCL...\n");
    if(queryOpenCL()){
        return true;
    }
    if(setupPlatform()){
        return true;
    }
    if(setupDevices()){
        return true;
    }
    S->addMessageF(1, "OpenCL is ready to work!\n");
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
        S->addMessageF(3, "Failure getting the number of platforms.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Get the array of platforms
    _platforms = new cl_platform_id[_num_platforms];
    if(!_platforms) {
        S->addMessageF(3, "Allocation memory error.\n");
        S->addMessage(0, "\tPlatforms array cannot be allocated\n");
        return true;
    }
    err_code = clGetPlatformIDs(_num_platforms, _platforms, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure getting the platforms list.\n");
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
            S->addMessageF(3, "Failure getting the number of devices.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        // Gets the devices array
        devices = new cl_device_id[num_devices];
        if(!devices) {
            S->addMessageF(3, "Allocation memory error.\n");
            S->addMessage(0, "\tDevices array cannot be allocated\n");
            return true;
        }
        err_code = clGetDeviceIDs(_platforms[i], CL_DEVICE_TYPE_ALL,
                                  num_devices, devices, NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessageF(3, "Failure getting the devices list.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        // Shows device arrays
        for(j = 0; j < num_devices; j++){
            // Identifier
            strcpy(msg, "");
            sprintf(msg, "\tDevice %u, Platform %u...\n", j, i);
            S->addMessage(1, msg);
            // Device name
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_NAME,
                                       1024 * sizeof(char),
                                       &aux,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                S->addMessageF(3, "Failure getting the device name.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            strcpy(msg, "");
            sprintf(msg, "\t\tDEVICE: %s\n", aux);
            S->addMessage(0, msg);
            // Platform vendor
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_VENDOR,
                                       1024 * sizeof(char),
                                       &aux,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                S->addMessageF(3, "Failure getting the device vendor.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            strcpy(msg, "");
            sprintf(msg, "\t\tVENDOR: %s\n", aux);
            S->addMessage(0, msg);
            // Device type
            cl_device_type dType;
            err_code = clGetDeviceInfo(devices[j],
                                       CL_DEVICE_TYPE,
                                       sizeof(cl_device_type),
                                       &dType,
                                       NULL);
            if(err_code != CL_SUCCESS) {
                S->addMessageF(3, "Failure getting the device type.\n");
                S->printOpenCLError(err_code);
                return true;
            }
            if(dType == CL_DEVICE_TYPE_CPU)
                S->addMessage(0, "\t\tTYPE: CL_DEVICE_TYPE_CPU\n");
            else if(dType == CL_DEVICE_TYPE_GPU)
                S->addMessage(0, "\t\tTYPE: CL_DEVICE_TYPE_GPU\n");
            else if(dType == CL_DEVICE_TYPE_ACCELERATOR)
                S->addMessage(0, "\t\tTYPE: CL_DEVICE_TYPE_ACCELERATOR\n");
            else if(dType == CL_DEVICE_TYPE_DEFAULT)
                S->addMessage(0, "\t\tTYPE: CL_DEVICE_TYPE_DEFAULT\n");
        }
        delete[] devices; devices = NULL;
    }
    return false;
}

bool CalcServer::setupPlatform()
{
    char msg[1024];
    InputOutput::ProblemSetup  *P = InputOutput::ProblemSetup::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    if(P->settings.platform_id >= _num_platforms){
        S->addMessageF(3, "Impossible to use the requested platform.\n");
        strcpy(msg, "");
        sprintf(msg, "\t%u platform has been selected, but just %u platforms have been found\n",
                P->settings.platform_id, _num_platforms);
        S->addMessage(0, msg);
        return true;
    }
    _platform = _platforms[P->settings.platform_id];
    return false;
}

bool CalcServer::setupDevices()
{
    cl_int err_code;
    cl_uint i;
    char msg[1024];
    InputOutput::ProblemSetup  *P = InputOutput::ProblemSetup::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    _devices = NULL;
    // Gets the number of valid devices
    err_code = clGetDeviceIDs(_platform,
                              P->settings.device_type,
                              0,
                              NULL,
                              &_num_devices);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure getting the number of devices.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(P->settings.device_id >= _num_devices) {
        S->addMessageF(3, "Impossible to use the selected device.\n");
        strcpy(msg, "");
        sprintf(msg, "\t%u device has been selected, but just %u devices are available\n",
                P->settings.device_id, _num_devices);
        S->addMessage(0, msg);
        if(P->settings.device_type == CL_DEVICE_TYPE_ALL)
            S->addMessage(0, "\t\tCL_DEVICE_TYPE_ALL filter activated\n");
        else if(P->settings.device_type == CL_DEVICE_TYPE_CPU)
            S->addMessage(0, "\t\tCL_DEVICE_TYPE_CPU filter activated\n");
        else if(P->settings.device_type == CL_DEVICE_TYPE_GPU)
            S->addMessage(0, "\t\tCL_DEVICE_TYPE_GPU filter activated\n");
        else if(P->settings.device_type == CL_DEVICE_TYPE_ACCELERATOR)
            S->addMessage(0, "\t\tCL_DEVICE_TYPE_ACCELERATOR filter activated\n");
        else if(P->settings.device_type == CL_DEVICE_TYPE_DEFAULT)
            S->addMessage(0, "\t\tCL_DEVICE_TYPE_DEFAULT filter activated\n");

        return true;
    }
    // Gets the devices array
    _devices = new cl_device_id[_num_devices];
    if(!_devices){
        S->addMessageF(3, "Allocation memory error.\n");
        S->addMessage(0, "\tDevices array can't be allocated\n");
        return true;
    }
    err_code = clGetDeviceIDs(_platform,
                              P->settings.device_type,
                              _num_devices,
                              _devices,
                              &_num_devices);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure getting the devices list.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Create a devices context
    _context = clCreateContext(0,
                               _num_devices,
                               _devices,
                               NULL,
                               NULL,
                               &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating an OpenCL context.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // Create command queues
    _command_queues = new cl_command_queue[_num_devices];
    if(_command_queues == NULL){
        S->addMessageF(3, "Allocation memory error.\n");
        S->addMessage(0, "\tCommand queues array cannot be allocated\n");
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
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }
    }
    // Store the selected ones
    _device = _devices[P->settings.device_id];
    _command_queue = _command_queues[P->settings.device_id];
    return false;
}

bool CalcServer::setup()
{
    unsigned int i, j;
    cl_uint err_code=0;
    char msg[512];
    InputOutput::Fluid *F = InputOutput::Fluid::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    strcpy(msg, "");

    // Check for the required variables that must be defined by the user
    if(!_vars->get("h")){
        S->addMessageF(3, "Missed kernel length \"h\".\n");
        return true;
    }
    if(strcmp(_vars->get("h")->type(), "float")){
        sprintf(msg,
                "Kernel length \"h\" must be of type \"float\", but \"%s\" has been specified\n",
                _vars->get("h")->type());
        S->addMessageF(3, msg);
        return true;
    }
    float h = *(float *)_vars->get("h")->get();
    if(h <= 0.f){
        sprintf(msg,
                "Kernel length \"h\" must be greater than 0, but \"%g\" has been set\n",
                h);
        S->addMessageF(3, msg);
        return true;
    }

    sprintf(msg, "Found kernel height: h = %g [m]\n", h);
    S->addMessageF(1, msg);

    // Setup the scalar data per particles sets
    for(i = 0; i < P->sets.size(); i++){
        InputOutput::ProblemSetup::sphParticlesSet* set = P->sets.at(i);
        for(j = 0; j < set->scalarNames().size(); j++){
            const char *name = set->scalarNames().at(j);
            const char *val = set->scalarValues().at(j);
            if(!_vars->get(name)){
                sprintf(msg,
                        "Variable \"%s\" has not been registered\n",
                        name);
                S->addMessageF(3, msg);
                sprintf(msg, "Particles set: %u\n", i);
                S->addMessageF(0, msg);
                return true;
            }
            if(!strchr(_vars->get(name)->type(), '*')){
                sprintf(msg,
                        "Variable \"%s\" has been registered as a scalar, however it is established per particles set\n",
                        name);
                S->addMessageF(3, msg);
                sprintf(msg, "Particles set: %u\n", i);
                S->addMessageF(0, msg);
                return true;
            }
            InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable *)_vars->get(name);
            size_t typesize = _vars->typeToBytes(_vars->get(name)->type());
            size_t len = _vars->get(name)->size() / typesize;
            if(len != P->sets.size()){
                sprintf(msg,
                        "Variable \"%s\" is an array of %u components, which does not match with the number of particles sets (n_sets = %u)\n",
                        name,
                        len,
                        P->sets.size());
                S->addMessageF(3, msg);
                return true;
            }
            void *data = malloc(typesize);
            if(_vars->solve(_vars->get(name)->type(), val, data)){
                sprintf(msg, "Particles set: %u\n", i);
                S->addMessageF(0, msg);
                return true;
            }
            cl_mem mem = *(cl_mem*)_vars->get(name)->get();
            cl_int status;
            status = clEnqueueWriteBuffer(_command_queue, mem, CL_TRUE,
                                          i * typesize, typesize, data,
                                          0, NULL, NULL);
            free(data); data = NULL;
            if(status != CL_SUCCESS) {
                sprintf(msg,
                        "Failure sending variable \"%s\" to particles set %u\n",
                        name,
                        i);
                S->addMessageF(3, msg);
                S->printOpenCLError(status);
                return true;
            }
        }
    }

    // Setup the tools
    for(i = 0; i < _tools.size(); i++){
        if(_tools.at(i)->setup()){
            return true;
        }
    }

    S->addMessageF(1, "Calculation server is ready! ;-) \n");
    return false;
}

}}  // namespace
