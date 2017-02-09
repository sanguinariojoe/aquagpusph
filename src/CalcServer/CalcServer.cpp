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

CalcServer::CalcServer()
    : _base_path("")
    , _num_platforms(0)
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
    unsigned int i, j;
    char msg[1024];
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    _current_tool_name = new char[256];
    strcpy(_current_tool_name, "");
    if(setupOpenCL()) {
        exit(EXIT_FAILURE);
    }

    _base_path = P->settings.base_path;

    unsigned int num_sets = P->sets.size();
    unsigned int N = 0;
    for(i = 0; i < P->sets.size(); i++) {
        N += P->sets.at(i)->n();
    }

    unsigned int num_icell = nextPowerOf2(N);
    num_icell = roundUp(num_icell, _ITEMS*_GROUPS);

    _vars = new InputOutput::Variables();

    // Register default scalars
    char val[64];
    char len[16];
    strcpy(len, "");
    unsigned int dims = 2;
    #ifdef HAVE_3D
        dims = 3;
    #endif
    sprintf(val, "%u", dims);
    if(_vars->registerVariable("dims", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%g", 0.f);
    if(_vars->registerVariable("t", "float", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%g", 0.f);
    if(_vars->registerVariable("dt", "float", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", 0);
    if(_vars->registerVariable("iter", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", 0);
    if(_vars->registerVariable("frame", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%g", std::numeric_limits<float>::max());
    if(_vars->registerVariable("end_t", "float", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", std::numeric_limits<unsigned int>::max());
    if(_vars->registerVariable("end_iter", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", std::numeric_limits<unsigned int>::max());
    if(_vars->registerVariable("end_frame", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", N);
    if(_vars->registerVariable("N", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", P->sets.size());
    if(_vars->registerVariable("n_sets", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    sprintf(val, "%u", num_icell);
    if(_vars->registerVariable("n_radix", "unsigned int", len, val))
        exit(EXIT_FAILURE);
    // Number of cells in x, y, z directions, and the total (n_x * n_y * n_z)
    strcpy(val, "0, 0, 0, 0");
    if(_vars->registerVariable("n_cells", "uivec4", len, val))
        exit(EXIT_FAILURE);
    // Kernel support
    sprintf(val, "%g", 2.f);
    if(_vars->registerVariable("support", "float", len, val))
        exit(EXIT_FAILURE);

    // Register default arrays
    strcpy(val, "");
    sprintf(len, "%u", N);
    if(_vars->registerVariable("id", "unsigned int*", len, val))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("r", "vec*", len, val))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("iset", "unsigned int*", len, val))
        exit(EXIT_FAILURE);
    sprintf(len, "%u", num_icell);
    if(_vars->registerVariable("id_sorted", "unsigned int*", len, val))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("id_unsorted", "unsigned int*", len, val))
        exit(EXIT_FAILURE);
    if(_vars->registerVariable("icell", "unsigned int*", len, val))
        exit(EXIT_FAILURE);
    sprintf(len, "n_cells_w");
    if(_vars->registerVariable("ihoc", "unsigned int*", len, val))
        exit(EXIT_FAILURE);

    // Register the user variables and arrays
    for(i = 0; i < P->variables.names.size(); i++){
        bool flag = _vars->registerVariable(P->variables.names.at(i),
                                            P->variables.types.at(i),
                                            P->variables.lengths.at(i),
                                            P->variables.values.at(i));
        if(flag){
            exit(EXIT_FAILURE);
        }
    }

    // Register the user definitions
    for(i = 0; i < P->definitions.names.size(); i++){
        size_t deflen=0;
        char *defstr = NULL;
        if(!strcmp(P->definitions.values.at(i), "")){
            // First case, named definitions
            deflen = strlen(P->definitions.names.at(i));
            defstr = new char[deflen + 3];
            if(!defstr){
                S->addMessageF(
                    3, "Failure allocating memory for the definition\n");
                sprintf(msg, "\t\"%s\"\n", P->definitions.names.at(i));
                S->addMessage(0, msg);
                exit(EXIT_FAILURE);
            }
            strcpy(defstr, "-D");
            strcat(defstr, P->definitions.names.at(i));
        }
        else if(!P->definitions.evaluations.at(i)){
            // Second case, valued definitions
            deflen = strlen(P->definitions.names.at(i)) +
                     strlen(P->definitions.values.at(i));
            defstr = new char[deflen + 4];
            if(!defstr){
                S->addMessageF(
                    3, "Failure allocating memory for the definition\n");
                sprintf(msg, "\t\"%s\"\n", P->definitions.names.at(i));
                S->addMessage(0, msg);
                exit(EXIT_FAILURE);
            }
            strcpy(defstr, "-D");
            strcat(defstr, P->definitions.names.at(i));
            strcat(defstr, "=");
            strcat(defstr, P->definitions.values.at(i));
        }
        else{
            // Third case, evaluated definitions
            deflen = strlen(P->definitions.names.at(i));
            defstr = new char[deflen + 1 + 4 + 128];
            if(!defstr){
                S->addMessageF(
                    3, "Failure allocating memory for the definition\n");
                sprintf(msg, "\t\"%s\"\n", P->definitions.names.at(i));
                S->addMessage(0, msg);
                exit(EXIT_FAILURE);
            }
            float defval = 0.f;
            if(_vars->solve("float",
                            P->definitions.values.at(i),
                            &defval))
            {
                exit(EXIT_FAILURE);
            }
            strcpy(defstr, "-D");
            strcat(defstr, P->definitions.names.at(i));
            sprintf(defstr + strlen(defstr), "=%#Gf", defval);
        }
        _definitions.push_back(defstr);
    }

    // Register the tools
    for(i = 0; i < P->tools.size(); i++){
        if(!strcmp(P->tools.at(i)->get("type"), "kernel")){
            Kernel *tool = new Kernel(P->tools.at(i)->get("name"),
                                      P->tools.at(i)->get("path"),
                                      P->tools.at(i)->get("entry_point"),
                                      P->tools.at(i)->get("n"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "copy")){
            Copy *tool = new Copy(P->tools.at(i)->get("name"),
                                  P->tools.at(i)->get("in"),
                                  P->tools.at(i)->get("out"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "python")){
            Python *tool = new Python(P->tools.at(i)->get("name"),
                                      P->tools.at(i)->get("path"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "set")){
            Set *tool = new Set(P->tools.at(i)->get("name"),
                                P->tools.at(i)->get("in"),
                                P->tools.at(i)->get("value"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "set_scalar")){
            SetScalar *tool = new SetScalar(P->tools.at(i)->get("name"),
                                            P->tools.at(i)->get("in"),
                                            P->tools.at(i)->get("value"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "reduction")){
            Reduction *tool = new Reduction(P->tools.at(i)->get("name"),
                                            P->tools.at(i)->get("in"),
                                            P->tools.at(i)->get("out"),
                                            P->tools.at(i)->get("operation"),
                                            P->tools.at(i)->get("null"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "link-list")){
            LinkList *tool = new LinkList(P->tools.at(i)->get("name"),
                                          P->tools.at(i)->get("in"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "radix-sort")){
            RadixSort *tool = new RadixSort(P->tools.at(i)->get("name"),
                                            P->tools.at(i)->get("in"),
                                            P->tools.at(i)->get("perm"),
                                            P->tools.at(i)->get("inv_perm"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "assert")){
            Assert *tool = new Assert(P->tools.at(i)->get("name"),
                                      P->tools.at(i)->get("condition"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->tools.at(i)->get("type"), "dummy")){
            Tool *tool = new Tool(P->tools.at(i)->get("name"));
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
    // Register the reporters
    for(i = 0; i < P->reports.size(); i++){
        if(!strcmp(P->reports.at(i)->get("type"), "screen")){
            bool bold = false;
            if(!strcmp(P->reports.at(i)->get("bold"), "true") ||
               !strcmp(P->reports.at(i)->get("bold"), "True")){
               bold = true;
            }
            Reports::Screen *tool = new Reports::Screen(
                P->reports.at(i)->get("name"),
                P->reports.at(i)->get("fields"),
                P->reports.at(i)->get("color"),
                bold);
            _tools.push_back(tool);
        }
        else if(!strcmp(P->reports.at(i)->get("type"), "file")){
            Reports::TabFile *tool = new Reports::TabFile(
                P->reports.at(i)->get("name"),
                P->reports.at(i)->get("fields"),
                P->reports.at(i)->get("path"));
            _tools.push_back(tool);
        }
        else if(!strcmp(P->reports.at(i)->get("type"), "particles")){
            // Get the first particle associated to this set
            unsigned int set_id = atoi(P->reports.at(i)->get("set"));
            if(set_id >= P->sets.size()){
                sprintf(msg,
                        "Report \"%s\" requested the particles set %u but just %lu can be found.\n",
                        P->reports.at(i)->get("name"),
                        set_id,
                        P->sets.size());
                S->addMessageF(3, msg);
                exit(EXIT_FAILURE);
            }
            unsigned int first = 0;
            for(j = 0; j < set_id; j++){
                first += P->sets.at(j)->n();
            }

            // And the ipf and fps
            unsigned int ipf = atoi(P->reports.at(i)->get("ipf"));
            float fps = atof(P->reports.at(i)->get("fps"));

            Reports::SetTabFile *tool = new Reports::SetTabFile(
                P->reports.at(i)->get("name"),
                P->reports.at(i)->get("fields"),
                first,
                P->sets.at(set_id)->n(),
                P->reports.at(i)->get("path"),
                ipf,
                fps);
            _tools.push_back(tool);
        }
        else if(!strcmp(P->reports.at(i)->get("type"), "performance")){
            bool bold = false;
            if(!strcmp(P->reports.at(i)->get("bold"), "true") ||
               !strcmp(P->reports.at(i)->get("bold"), "True")){
               bold = true;
            }
            Reports::Performance *tool = new Reports::Performance(
                P->reports.at(i)->get("name"),
                P->reports.at(i)->get("color"),
                bold,
                P->reports.at(i)->get("path"));
            _tools.push_back(tool);
        }
        else{
            sprintf(msg,
                    "Unrecognized report type \"%s\" when parsing the \"%s\".\n",
                    P->reports.at(i)->get("type"),
                    P->reports.at(i)->get("name"));
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

    S->addMessageF(1, "Destroying definitions...\n");
    for(i = 0; i < _definitions.size(); i++){
        delete[] _definitions.at(i);
        _definitions.at(i) = NULL;
    }
    _definitions.clear();

    S->addMessageF(1, "Destroying variables manager...\n");
    if(_vars) delete _vars; _vars=NULL;

    if(_current_tool_name) delete[] _current_tool_name; _current_tool_name=NULL;

    std::map<std::string, UnSort*>::iterator it = unsorters.begin();
    while(it != unsorters.end())
    {
        delete it->second;
        it++;
    }
}

bool CalcServer::update()
{
    InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    unsigned int i;
    while(!T->mustPrintOutput() && !T->mustStop()){
        S->initFrame();

        // Execute the tools
        strcpy(_current_tool_name, "__pre execution__");
        for(i = 0; i < _tools.size(); i++){
            strcpy(_current_tool_name, _tools.at(i)->name());
            if(_tools.at(i)->execute()){
                sleep(__ERROR_SHOW_TIME__);
                return true;
            }
        }
        strcpy(_current_tool_name, "__post execution__");

        // Key events
        while(isKeyPressed()){
            if(getchar() == 'c'){
                S->addMessageF(2, "Interrumption request detected.\n");
                sleep(__ERROR_SHOW_TIME__);
                return true;
            }
        }

        clFinish(command_queue());

        S->endFrame();
	}
	return false;
}

cl_event CalcServer::getUnsortedMem(const char* var_name,
                                    size_t offset,
                                    size_t cb,
                                    void *ptr)
{
    cl_int err_code;
    char msg[256];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    // Generate the unsorted if it does not exist yet
    UnSort *unsorter = NULL;
    if(unsorters.find(var_name) == unsorters.end()){
        unsorter = new UnSort(var_name, var_name);
        if(unsorter->setup()){
            delete unsorter;
            return NULL;
        }
        unsorters.insert(std::make_pair(var_name, unsorter));
    }
    // Get the unsorter
    unsorter = unsorters[var_name];
    if(unsorter->execute()){
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
        sprintf(msg,
                "Failure receiving the variable \"%s\" from server.\n",
                var_name);
        S->addMessageF(3, msg);
        S->printOpenCLError(err_code);
        return NULL;
    }
    return event;
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
            #ifdef CL_DEVICE_TYPE_CUSTOM
            else if(dType == CL_DEVICE_TYPE_CUSTOM)
                S->addMessage(0, "\t\tTYPE: CL_DEVICE_TYPE_CUSTOM\n");
            #endif
            else {
                sprintf(msg, "\t\tTYPE: %ul\n", dType);
                S->addMessage(0, msg);
            }
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

    S->addMessageF(3, msg);
    S->addMessage(0, errinfo);
    S->addMessage(0, "\n");
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
                               context_error_notify,
                               _current_tool_name,
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
                S->addMessage(0, msg);
                return true;
            }
            if(!strchr(_vars->get(name)->type(), '*')){
                sprintf(msg,
                        "Variable \"%s\" has been registered as a scalar, however it is established per particles set\n",
                        name);
                S->addMessageF(3, msg);
                sprintf(msg, "Particles set: %u\n", i);
                S->addMessage(0, msg);
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
                S->addMessage(0, msg);
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

    return false;
}

}}  // namespace
