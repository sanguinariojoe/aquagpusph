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

#include <limits>
#include <sys/time.h>

#include <CalcServer.h>
#include <AuxiliarMethods.h>
#include <FileManager.h>
#include <ProblemSetup.h>
#include <Fluid.h>
#include <TimeManager.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{

CalcServer::CalcServer()
	: platforms(NULL)
	, devices(NULL)
	, command_queues(NULL)
	, energy_computed(false)
	, bounds_computed(false)
	, fluid_mass(0.f)
{
	unsigned int i;
	char msg[1024];
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	verbose_level = P->settings.verbose_level;

	if(setupOpenCL()) {
	    exit(255);
	}

	num_fluids  = P->n_fluids;
	num_sensors = P->SensorsParameters.pos.size();
	n = 0;
	for(i=0;i<P->n_fluids;i++) {
	    n += P->fluids[i].n;
	}
	N = n + num_sensors;

	allocated_mem = 0;
	imove = allocMemory(N * sizeof( cl_int ));
	if(!imove) exit(255);
	imovein = allocMemory(N * sizeof( cl_int ));
	if(!imovein) exit(255);
	ifluid = allocMemory(N * sizeof( cl_int ));
	if(!ifluid) exit(255);
	ifluidin = allocMemory(N * sizeof( cl_int ));
	if(!ifluidin) exit(255);
	pos = allocMemory(N * sizeof( vec ));
	if(!pos) exit(255);
	normal = allocMemory(N * sizeof( vec ));
	if(!normal) exit(255);
	v = allocMemory(N * sizeof( vec ));
	if(!v) exit(255);
	f = allocMemory(N * sizeof( vec ));
	if(!f) exit(255);
	dens = allocMemory(N * sizeof( cl_float ));
	if(!dens) exit(255);
	drdt = allocMemory(N * sizeof( cl_float ));
	if(!drdt) exit(255);
	drdt_F = allocMemory(N * sizeof( cl_float ));
	if(!drdt_F) exit(255);
	mass = allocMemory(N * sizeof( cl_float ));
	if(!mass) exit(255);
	hp = allocMemory(N * sizeof( cl_float ));
	if(!hp) exit(255);
	press = allocMemory(N * sizeof( cl_float ));
	if(!press) exit(255);
	posin = allocMemory(N * sizeof( vec ));
	if(!posin) exit(255);
	normalin = allocMemory(N * sizeof( vec ));
	if(!normalin) exit(255);
	vin = allocMemory(N * sizeof( vec ));
	if(!vin) exit(255);
	fin = allocMemory(N * sizeof( vec ));
	if(!fin) exit(255);
	densin = allocMemory(N * sizeof( cl_float ));
	if(!densin) exit(255);
	drdtin = allocMemory(N * sizeof( cl_float ));
	if(!drdtin) exit(255);
	massin = allocMemory(N * sizeof( cl_float ));
	if(!massin) exit(255);
	hpin = allocMemory(N * sizeof( cl_float ));
	if(!hpin) exit(255);
	pressin = allocMemory(N * sizeof( cl_float ));
	if(!pressin) exit(255);
	sigma = allocMemory(N * sizeof( cl_float ));
	if(!sigma) exit(255);
	dtconv = allocMemory(N * sizeof( cl_float ));
	if(!dtconv) exit(255);
	shepard = allocMemory(N * sizeof( cl_float ));
	if(!shepard) exit(255);
	shepard_gradient = allocMemory(N * sizeof( vec ));
	if(!shepard_gradient) exit(255);

	num_icell = nextPowerOf2(N);
	num_icell = roundUp(num_icell, _ITEMS*_GROUPS);
	permutation = allocMemory(num_icell * sizeof( cl_uint ));
	if(!permutation) exit(255);
	permutation_inverse = allocMemory(num_icell * sizeof( cl_uint ));
	if(!permutation_inverse) exit(255);
	icell = allocMemory(num_icell * sizeof( cl_uint ));
	if(!icell) exit(255);
	ihoc = NULL;               // ihoc must be allocated later
	cell_has_particles = NULL; // cell_has_particles must be allocated later

	gamma = allocMemory(num_fluids * sizeof( cl_float ));
	if(!gamma) exit(255);
	refd = allocMemory(num_fluids * sizeof( cl_float ));
	if(!refd) exit(255);
	visc_dyn = allocMemory(num_fluids * sizeof( cl_float ));
	if(!visc_dyn) exit(255);
	visc_kin = allocMemory(num_fluids * sizeof( cl_float ));
	if(!visc_kin) exit(255);
	visc_dyn_corrected = allocMemory(num_fluids * sizeof( cl_float ));
	if(!visc_dyn_corrected) exit(255);
	delta = allocMemory(num_fluids * sizeof( cl_float ));
	if(!delta) exit(255);

	DT = NULL;     // Reduction auxiliar time step must be allocated later
	sensor_mode = allocMemory(max(num_sensors,(uint)1) * sizeof( cl_ushort ));
	if(!sensor_mode) exit(255);

	sprintf(msg, "Allocated memory = %u bytes\n", (unsigned int)allocated_mem);
	S->addMessageF(1, msg);
	// Create the computation tools
	predictor       = new Predictor();
	grid            = new Grid();
	link_list       = new LinkList();
	rates           = new Rates();
	elastic_bounce  = new Boundary::ElasticBounce();
	de_Leffe        = new Boundary::DeLeffe();
	ghost_particles = new Boundary::GhostParticles();
	shepard_tool    = new Shepard();
	corrector       = new Corrector();
	domain          = new Domain();
	time_step       = new TimeStep();
	dens_int        = new DensityInterpolation();
	sensors         = new Sensors();
	energy_tool     = new Energy();
	bounds_tool     = new Bounds();
	motions.clear();
	for(i=0;i<P->motions.size();i++){
	    if(P->motions.at(i)->type == 0){
	        Movement::Quaternion *Move = new Movement::Quaternion();
	        motions.push_back(Move);
	    }
	    if(P->motions.at(i)->type == 1){
	        Movement::LIQuaternion *Move = new Movement::LIQuaternion();
	        motions.push_back(Move);
	    }
	    if(P->motions.at(i)->type == 2){
	        Movement::C1Quaternion *Move = new Movement::C1Quaternion();
	        motions.push_back(Move);
	    }
	    if(P->motions.at(i)->type == 3){
	        Movement::ScriptQuaternion *Move = new Movement::ScriptQuaternion();
	        motions.push_back(Move);
	    }
	}
	portals.clear();
	for(i=0;i<P->portals.size();i++){
	    Portal::Portal *portal = new Portal::Portal(P->portals.at(i));
	    portals.push_back(portal);
	}
}

CalcServer::~CalcServer()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    S->addMessageF(1, "Destroying the predictor manager...\n");
	delete predictor; predictor=NULL;
    S->addMessageF(1, "Destroying the grid manager...\n");
	delete grid; grid=NULL;
    S->addMessageF(1, "Destroying the link list manager...\n");
	delete link_list; link_list=NULL;
    S->addMessageF(1, "Destroying the rates manager...\n");
	delete rates; rates=NULL;
    S->addMessageF(1, "Destroying the elasticBounce boundary condition manager...\n");
	delete elastic_bounce; elastic_bounce=NULL;
    S->addMessageF(1, "Destroying the DeLeffe boundary condition manager...\n");
	delete de_Leffe; de_Leffe=NULL;
    S->addMessageF(1, "Destroying the ghost particles manager...\n");
	delete ghost_particles; ghost_particles=NULL;
    S->addMessageF(1, "Destroying the 0th order correction manager...\n");
	delete shepard_tool; shepard_tool=NULL;
    S->addMessageF(1, "Destroying the corrector manager...\n");
	delete corrector; corrector=NULL;
    S->addMessageF(1, "Destroying the domain bounds manager...\n");
	delete domain; domain=NULL;
    S->addMessageF(1, "Destroying the time step manager...\n");
	delete time_step; time_step=NULL;
    S->addMessageF(1, "Destroying the density interpolator...\n");
	delete dens_int; dens_int=NULL;
    S->addMessageF(1, "Destroying the sensors manager...\n");
	delete sensors; sensors=NULL;
    S->addMessageF(1, "Destroying the energy computation tool...\n");
	delete energy_tool; energy_tool=NULL;
    S->addMessageF(1, "Destroying the bounds computation tool...\n");
	delete bounds_tool; bounds_tool=NULL;
    S->addMessageF(1, "Destroying the movement managers...\n");
	for(i=0;i<motions.size();i++) {
	    delete motions.at(i);
	}
	motions.clear();
    S->addMessageF(1, "Sutting down the OpenCL context...\n");
	if(context)clReleaseContext(context);
	for(i=0;i<num_devices;i++) {
	    if(command_queues[i])clReleaseCommandQueue(command_queues[i]);
	}
    S->addMessageF(1, "Deallocating memory from the devices...\n");
	if(imove)clReleaseMemObject(imove); imove=NULL;
	if(imovein)clReleaseMemObject(imovein); imovein=NULL;
	if(ifluid)clReleaseMemObject(ifluid); ifluid=NULL;
	if(ifluidin)clReleaseMemObject(ifluidin); ifluidin=NULL;
	if(pos)clReleaseMemObject(pos); pos=NULL;
	if(normal)clReleaseMemObject(normal); normal=NULL;
	if(v)clReleaseMemObject(v); v=NULL;
	if(f)clReleaseMemObject(f); f=NULL;
	if(dens)clReleaseMemObject(dens); dens=NULL;
	if(drdt)clReleaseMemObject(drdt); drdt=NULL;
	if(drdt_F)clReleaseMemObject(drdt_F); drdt_F=NULL;
	if(mass)clReleaseMemObject(mass); mass=NULL;
	if(hp)clReleaseMemObject(hp); hp=NULL;
	if(press)clReleaseMemObject(press); press=NULL;
	if(posin)clReleaseMemObject(posin); posin=NULL;
	if(normalin)clReleaseMemObject(normalin); normalin=NULL;
	if(vin)clReleaseMemObject(vin); vin=NULL;
	if(fin)clReleaseMemObject(fin); fin=NULL;
	if(densin)clReleaseMemObject(densin); densin=NULL;
	if(drdtin)clReleaseMemObject(drdtin); drdtin=NULL;
	if(massin)clReleaseMemObject(massin); massin=NULL;
	if(hpin)clReleaseMemObject(hpin); hpin=NULL;
	if(pressin)clReleaseMemObject(pressin); pressin=NULL;
	if(sigma)clReleaseMemObject(sigma); sigma=NULL;
	if(dtconv)clReleaseMemObject(dtconv); dtconv=NULL;
	if(shepard)clReleaseMemObject(shepard); shepard=NULL;
	if(shepard_gradient)clReleaseMemObject(shepard_gradient); shepard_gradient=NULL;
	if(permutation)clReleaseMemObject(permutation); permutation=NULL;
	if(permutation_inverse)clReleaseMemObject(permutation_inverse); permutation_inverse=NULL;
	if(icell)clReleaseMemObject(icell); icell=NULL;
	if(ihoc)clReleaseMemObject(ihoc); ihoc=NULL;
	if(cell_has_particles)clReleaseMemObject(cell_has_particles); cell_has_particles=NULL;
	if(gamma)clReleaseMemObject(gamma); gamma=NULL;
	if(refd)clReleaseMemObject(refd); refd=NULL;
	if(visc_dyn)clReleaseMemObject(visc_dyn); visc_dyn=NULL;
	if(visc_kin)clReleaseMemObject(visc_kin); visc_kin=NULL;
	if(visc_dyn_corrected)clReleaseMemObject(visc_dyn_corrected); visc_dyn_corrected=NULL;
	if(delta)clReleaseMemObject(delta); delta=NULL;
	if(DT)clReleaseMemObject(DT); DT=NULL;
	if(sensor_mode)clReleaseMemObject(sensor_mode); sensor_mode=NULL;
    S->addMessageF(1, "Deallocating host memory...\n");
	if(platforms) delete[] platforms; platforms=NULL;
	if(devices) delete[] devices; devices=NULL;
	if(command_queues) delete[] command_queues; command_queues=NULL;
}

bool CalcServer::update()
{
	InputOutput::TimeManager *T   = InputOutput::TimeManager::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	unsigned int i;
	while(!T->mustPrintOutput() && !T->mustStop()){
	    if(predictor->execute())
	        return true;
	    if(link_list_step >= link_list_steps){
	        link_list_step = 0;
	        if(grid->execute())
	            return true;
	        if(link_list->execute())
	            return true;
	    }
	    link_list_step++;
	    if(rates->execute())
	        return true;
        // Since the density interpolation will not take into account the
        // boundaries, we must perform it before the shepard term is corrected.
	    if(T->time() >= 0.f){
            if(dens_int_step >= dens_int_steps){
                dens_int_step = 0;
                if(dens_int->execute())
                    return true;
            }
            dens_int_step++;
	    }
	    if(de_Leffe->execute())
	        return true;
	    if(ghost_particles->execute())
	        return true;
	    if(shepard_tool->execute())
	        return true;
	    if(elastic_bounce->execute())
	        return true;
	    if(corrector->execute())
	        return true;
	    for(i=0;i<portals.size();i++){
	        if(portals.at(i)->execute())
	            return true;
	    }
	    if(domain->execute())
	        return true;
	    if(time_step->execute())
	        return true;
	    if(T->time() >= 0.f){
	        for(i=0;i<motions.size();i++){
	            if(motions.at(i)->execute())
	                return true;
	        }
	        if(sensors->execute())
	            return true;
	        if(T->mustPrintLog()) {
	            printLog();
	        }
	        if(T->mustPrintEnergy()) {
	            printEnergy();
	        }
	        if(T->mustPrintBounds()) {
	            printBounds();
	        }
	    }
		T->update(dt);
		S->update();
		energy_computed=false;
		bounds_computed=false;
		// Key events
		while(isKeyPressed()){
	        if(getchar() == 'c'){
	            S->addMessageF(1, "Interrumption request detected.\n");
	            return true;
	        }
		}
	}
	return false;
}

bool CalcServer::getData(void *dest, cl_mem orig, size_t size, size_t offset)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	cl_int err_code = clEnqueueReadBuffer(command_queue, orig, CL_TRUE, offset,
                                          size, dest, 0, NULL, NULL);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Failure retrieving memory from the server.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	return false;
}

bool CalcServer::sendData(cl_mem dest, void* orig, size_t size)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	cl_int err_code;
	err_code  = clEnqueueWriteBuffer(command_queue, dest, CL_TRUE, 0, size, orig, 0, NULL, NULL);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Failure sending memory to the server.\n");
	    S->printOpenCLError(err_code);
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
	if(getPlatform()){
	    return true;
	}
	if(getDevices()){
	    if(platforms) delete[] platforms; platforms=0;
	    return true;
	}
	S->addMessageF(1, "OpenCL is ready to work!\n");
	return false;
}

bool CalcServer::queryOpenCL()
{
	cl_int err_code;
	cl_uint i,j,num_devices=0;
	cl_device_id *devices;
	char msg[1024], aux[1024];
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	platforms = NULL;
	strcpy(msg, "");
	// Gets the total number of platforms
	err_code = clGetPlatformIDs(0, NULL, &num_platforms);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't take the number of platforms.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	// Get array of patforms
	platforms = new cl_platform_id[num_platforms];
	if(!platforms) {
	    S->addMessageF(3, "Allocation memory error.\n");
	    S->addMessage(0, "\tPlatforms array can't be allocated\n");
	    return true;
	}
	err_code = clGetPlatformIDs(num_platforms, platforms, NULL);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't take the platforms list.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	for(i=0;i<num_platforms;i++){
	    // Get the number of devices
	    err_code = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL,
                                  &num_devices);
	    if(err_code != CL_SUCCESS) {
	        S->addMessageF(3, "Can't take the number of devices.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    // Gets the devices array
	    devices = new cl_device_id[num_devices];
	    if(!devices) {
	        S->addMessageF(3, "Allocation memory error.\n");
	        S->addMessage(0, "\tDevices array can't be allocated\n");
	        return true;
	    }
	    err_code = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                                  num_devices, devices, &num_devices);
	    if(err_code != CL_SUCCESS) {
	        S->addMessageF(3, "Can't take the devices list.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    // Shows device arrays
	    for(j=0;j<num_devices;j++){
	        // Identifier
	        strcpy(msg, "");
	        sprintf(msg, "\tDevice %u, Platform %u...\n", j, i);
	        S->addMessage(1, msg);
	        // Device name
	        err_code = clGetDeviceInfo(devices[j], CL_DEVICE_NAME,
                                       1024*sizeof(char), &aux, NULL);
	        if(err_code != CL_SUCCESS) {
	            S->addMessageF(3, "Can't get the device name.\n");
                S->printOpenCLError(err_code);
	            return true;
	        }
	        strcpy(msg, "");
	        sprintf(msg, "\t\tDEVICE: %s\n", aux);
	        S->addMessage(0, msg);
	        // Platform vendor
	        err_code = clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR,
                                       1024*sizeof(char), &aux, NULL);
	        if(err_code != CL_SUCCESS) {
	            S->addMessageF(3, "Can't get the device vendor.\n");
                S->printOpenCLError(err_code);
	            return true;
	        }
	        strcpy(msg, "");
	        sprintf(msg, "\t\tVENDOR: %s\n", aux);
	        S->addMessage(0, msg);
	        // Device type
	        cl_device_type dType;
	        err_code = clGetDeviceInfo(devices[j], CL_DEVICE_TYPE,
                                       sizeof(cl_device_type), &dType, NULL);
	        if(err_code != CL_SUCCESS) {
	            S->addMessageF(3, "Can't get the device type.\n");
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

bool CalcServer::getPlatform()
{
	char msg[1024];
	InputOutput::ProblemSetup  *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	if(P->settings.platform_id >= num_platforms){
	    S->addMessageF(3, "Can't find the requested platform.\n");
	    strcpy(msg, "");
	    sprintf(msg, "\t%u platform requested, but just %u platforms can be found\n",
                P->settings.platform_id, num_platforms);
	    S->addMessage(0, msg);
	    return true;
	}
	platform = platforms[P->settings.platform_id];
	return false;
}

bool CalcServer::getDevices()
{
	cl_int err_code;
	cl_uint i;
	char msg[1024];
	InputOutput::ProblemSetup  *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	devices = NULL;
	// Gets the number of valid devices
	err_code = clGetDeviceIDs(platform, P->settings.device_type, 0, NULL,
                              &num_devices);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't take the number of devices.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	if(P->settings.device_id >= num_devices) {
	    S->addMessageF(3, "Can't find the requested device.\n");
	    strcpy(msg, "");
	    sprintf(msg, "\t%u device requested, but just %u devices have been found\n",
                P->settings.device_id, num_devices);
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
	devices = new cl_device_id[num_devices];
	if(!devices) {
	    S->addMessageF(3, "Allocation memory error.\n");
	    S->addMessage(0, "\tDevices array can't be allocated\n");
	    return true;
	}
	err_code = clGetDeviceIDs(platform, P->settings.device_type, num_devices,
                              devices, &num_devices);
	if(err_code != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't take the devices list.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	// Create a devices context
	context = clCreateContext(0, num_devices, devices, NULL, NULL, &err_code);
	if(err_code != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't create an OpenCL context.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	// Create command queues
	command_queues = new cl_command_queue[num_devices];
	if(command_queues == NULL) {
	    S->addMessage(3, "(CalcServer::getDevices): Allocation memory error.\n");
	    S->addMessage(0, "\tCommand queues array can't be allocated\n");
	    return true;
	}
	for(i=0;i<num_devices;i++) {
	    #ifdef HAVE_GPUPROFILE
	        command_queues[i] = clCreateCommandQueue(context, devices[i],
                                                     CL_QUEUE_PROFILING_ENABLE,
                                                     &err_code);
	    #else
	        command_queues[i] = clCreateCommandQueue(context, devices[i], 0,
                                                     &err_code);
	    #endif
	    if(err_code != CL_SUCCESS) {
	        strcpy(msg, "");
	        sprintf(msg, "Can't create a command queue for the device %u.\n",i);
	        S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
	        return true;
	    }
	}
	// Store the selected ones
	device   = devices[P->settings.device_id];
	command_queue = command_queues[P->settings.device_id];
	return false;
}

cl_mem CalcServer::allocMemory(size_t size)
{
	int err_code;
	cl_mem mem_obj;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &err_code);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Allocation failure.\n");
	    S->printOpenCLError(err_code);
	    return NULL;
	}

	allocated_mem += size;
	return mem_obj;
}

bool CalcServer::fillFluid()
{
	unsigned int i;
	int err_code = CL_SUCCESS;
	char msg[512], file_name[256];
	cl_program program;
	cl_kernel kernel;
	size_t global_work_size;
	size_t local_work_size;
	unsigned int num_parts_parsed=0;
	InputOutput::Fluid *F = InputOutput::Fluid::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	strcpy(msg, "");

	for(i=0;i<num_fluids;i++) {
	    // Set the script file
	    if(!strlen(P->fluids[i].Script)) {
	        continue;
	    }
	    else {
	        strcpy(P->fluids[i].path, "");
	    }
	    strcpy(file_name, P->fluids[i].Path);
	    strcat(file_name, "/");
	    strcat(file_name, P->fluids[i].Script);
	    strcat(file_name, ".cl");
	    if(!loadKernelFromFile(&kernel, &program, context, device, file_name, "ParticlesDistribution", ""))
	        return true;
	    err_code |= sendArgument(kernel,  0, sizeof(cl_mem ), (void*)&imove);
	    err_code |= sendArgument(kernel,  1, sizeof(cl_mem ), (void*)&ifluid);
	    err_code |= sendArgument(kernel,  2, sizeof(cl_mem ), (void*)&pos);
	    err_code |= sendArgument(kernel,  3, sizeof(cl_mem ), (void*)&normal);
	    err_code |= sendArgument(kernel,  4, sizeof(cl_mem ), (void*)&v);
	    err_code |= sendArgument(kernel,  5, sizeof(cl_mem ), (void*)&f);
	    err_code |= sendArgument(kernel,  6, sizeof(cl_mem ), (void*)&dens);
	    err_code |= sendArgument(kernel,  7, sizeof(cl_mem ), (void*)&drdt);
	    err_code |= sendArgument(kernel,  8, sizeof(cl_mem ), (void*)&mass);
	    err_code |= sendArgument(kernel,  9, sizeof(cl_mem ), (void*)&hp);
	    err_code |= sendArgument(kernel, 10, sizeof(cl_uint), (void*)&num_parts_parsed);
	    err_code |= sendArgument(kernel, 11, sizeof(cl_uint), (void*)&(P->fluids[i].n));
	    err_code |= sendArgument(kernel, 12, sizeof(cl_uint), (void*)&i);
	    err_code |= sendArgument(kernel, 13, sizeof(vec    ), (void*)&(P->SPH_opts.deltar));
	    err_code |= sendArgument(kernel, 14, sizeof(cl_float), (void*)&(P->fluids[i].refd));
	    err_code |= sendArgument(kernel, 15, sizeof(cl_float), (void*)&(P->SPH_opts.h));
	    err_code |= sendArgument(kernel, 16, sizeof(cl_float), (void*)&(P->SPH_opts.cs));
	    num_parts_parsed += P->fluids[i].n;
	    if(err_code)
	        return true;
	    local_work_size = 256;
	    global_work_size = getGlobalWorkSize(P->fluids[i].n, local_work_size);
	    sprintf(msg, "Launching the initialization OpenCL script.\n");
	    S->addMessageF(1, msg);
	    sprintf(msg, "\tLocal work size = %lu.\n", local_work_size);
	    S->addMessage(0, msg);
	    sprintf(msg, "\tGlobal work size = %lu.\n", global_work_size);
	    S->addMessage(0, msg);
	    err_code = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	    if(err_code != CL_SUCCESS) {
	        S->addMessageF(3, "Can't execute the kernel.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    if(kernel)clReleaseKernel(kernel);
	    if(program)clReleaseProgram(program);
	}
	S->addMessageF(1, "Finished the work, reading back the info...\n");

	err_code |= getData((void*)F->imove, imove, sizeof(cl_int)*n);
	err_code |= getData((void*)F->ifluid, ifluid, sizeof(cl_int)*n);
	err_code |= getData((void*)F->pos, pos, sizeof(vec)*n);
	err_code |= getData((void*)F->normal, normal, sizeof(vec)*n);
	err_code |= getData((void*)F->v, v, sizeof(vec)*n);
	err_code |= getData((void*)F->f, f, sizeof(vec)*n);
	err_code |= getData((void*)F->dens, dens, sizeof(cl_float)*n);
	err_code |= getData((void*)F->drdt, drdt, sizeof(cl_float)*n);
	err_code |= getData((void*)F->mass, mass, sizeof(cl_float)*n);
	err_code |= getData((void*)F->hp, hp, sizeof(cl_float)*n);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Failure retrieving memory from the server.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	sprintf(msg, "%u particles have been transfered from the server to the host.\n", n);
	S->addMessageF(1,msg);
	S->addMessageF(1, "All the fluids are running! ;-)\n");

	return false;
}

bool CalcServer::setup()
{
	unsigned int i;
	cl_uint err_code=0;
	char msg[512];
	InputOutput::Fluid *F = InputOutput::Fluid::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	strcpy(msg, "");
	// Fluids data
    S->addMessageF(1, "Sending fluids data to the server...\n");
	cl_float *Gamma   = new cl_float[num_fluids];
	cl_float *Refd    = new cl_float[num_fluids];
	cl_float *Viscdyn = new cl_float[num_fluids];
	cl_float *Visckin = new cl_float[num_fluids];
	cl_float *ViscdynCorr = new cl_float[num_fluids];
	cl_float *Delta = new cl_float[num_fluids];
	for(i=0;i<num_fluids;i++) {
		Gamma[i] = P->fluids[i].gamma;
		Refd[i] = P->fluids[i].refd;
		Viscdyn[i] = P->fluids[i].visc_dyn;
		Visckin[i] = P->fluids[i].visc_kin;
		ViscdynCorr[i] = P->fluids[i].visc_dyn_corrected;
		Delta[i] = P->fluids[i].delta;
	}
	err_code |= sendData(gamma, Gamma, sizeof(cl_float)*num_fluids);
	err_code |= sendData(refd, Refd, sizeof(cl_float)*num_fluids);
	err_code |= sendData(visc_dyn, Viscdyn, sizeof(cl_float)*num_fluids);
	err_code |= sendData(visc_kin, Visckin, sizeof(cl_float)*num_fluids);
	err_code |= sendData(visc_dyn_corrected, ViscdynCorr, sizeof(cl_float)*num_fluids);
	err_code |= sendData(delta, Delta, sizeof(cl_float)*num_fluids);
	delete[] Gamma; Gamma=NULL;
	delete[] Refd; Refd=NULL;
	delete[] Viscdyn; Viscdyn=NULL;
	delete[] Visckin; Visckin=NULL;
	delete[] ViscdynCorr; ViscdynCorr=NULL;
	delete[] Delta; Delta=NULL;
	if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't send the fluids data to the server.\n");
	    return true;
	}
	// Variables
	g          = P->SPH_opts.g;
	hfac       = P->SPH_opts.hfac;
	dt_divisor = P->SPH_opts.dt_divisor;
	h          = P->SPH_opts.h;
	cs         = P->SPH_opts.cs;
	link_list_steps = P->SPH_opts.link_list_steps;
	dens_int_steps  = P->SPH_opts.dens_int_steps;
	// Calculate cell_length_factor
	float sep;
	#ifdef __CUBIC_KERNEL_TYPE__
		sep = 2.0f;
	#elif defined(__GAUSS_KERNEL_TYPE__)
		sep = 3.0f;
	#else   // Wendland
		sep = 2.0f;
	#endif
	float dist = sep*h;             // Minimum cell size.
	float ddt  = 0.1f*h/dt_divisor; // Maximum distance that a particle can move in a step.
	cell_length_factor = 1.f + (link_list_steps-1) * ddt / dist;
	sprintf(msg, "Cells size increased with %g factor.\n", cell_length_factor);
    S->addMessageF(1, msg);
	link_list_step = link_list_steps;
	dens_int_step  = 0;
	// Particles data
    S->addMessageF(1, "Sending the particles data to server...\n");
	err_code  = sendData(imove, F->imove, sizeof(cl_int)*N);
	err_code |= sendData(imovein, F->imove, sizeof(cl_int)*N);
	err_code |= sendData(ifluid, F->ifluid, sizeof(cl_int)*N);
	err_code |= sendData(ifluidin, F->ifluid, sizeof(cl_int)*N);
	err_code |= sendData(posin, F->pos, sizeof(vec)*N);
	err_code |= sendData(normalin, F->normal, sizeof(vec)*N);
	err_code |= sendData(vin, F->v, sizeof(vec)*N);
	err_code |= sendData(densin, F->dens, sizeof(cl_float)*N);
	err_code |= sendData(hpin, F->hp, sizeof(cl_float)*N);
	err_code |= sendData(drdtin, F->drdt, sizeof(cl_float)*N);
	err_code |= sendData(fin, F->f, sizeof(vec)*N);
	err_code |= sendData(massin, F->mass, sizeof(cl_float)*N);
	err_code |= sendData(pressin, F->press, sizeof(cl_float)*N);
	err_code |= sendData(pos, F->pos, sizeof(vec)*N);
	err_code |= sendData(normal, F->normal, sizeof(vec)*N);
	err_code |= sendData(v, F->v, sizeof(vec)*N);
	err_code |= sendData(dens, F->dens, sizeof(cl_float)*N);
	err_code |= sendData(drdt, F->drdt, sizeof(cl_float)*N);
	err_code |= sendData(hp, F->hp, sizeof(cl_float)*N);
	err_code |= sendData(f, F->f, sizeof(vec)*N);
	err_code |= sendData(mass, F->mass, sizeof(cl_float)*N);
	err_code |= sendData(press, F->press, sizeof(cl_float)*N);
	cl_float *Sigma = new cl_float[N];
	for(i=0;i<N;i++)
	    Sigma[i] = 1000.f;
	err_code |= sendData(sigma, Sigma, N*sizeof(cl_float));
	delete[] Sigma; Sigma=0;
	if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't send the particles data to the server.\n");
	    return true;
	}
	cl_ushort *mode = new cl_ushort[num_sensors];
	for(i=0;i<N-n;i++)
	    mode[i] = P->SensorsParameters.mod.at(i);
	err_code |= sendData(sensor_mode, mode, sizeof(cl_ushort)*(N-n));
	delete[] mode; mode = 0;

	num_cells = 0;
	num_cells_allocated = 0;
	dt = 0;

	for(i=0;i<motions.size();i++){
	    if(motions.at(i)->parse(P->motions.at(i)->path))
	        return true;
	}

	for(i=0;i<N;i++) {
	    if(F->imove[i] > 0)
	        fluid_mass += F->mass[i];
	}
    S->addMessageF(1, "Calculation server is ready! ;-) \n");
	return false;
}

void CalcServer::printLog()
{
	struct timeval now_time;
	char date[64], *aux;
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *files = InputOutput::FileManager::singleton();

	energy();

	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	strcpy(date, ctime(&seconds));
	aux = strstr(date, "\n");
	if(aux) strcpy(aux, "");

	fprintf(files->logFile(),"<i>%s: Printed file (%d)</i><br>\n",date,T->frame());
	fprintf(files->logFile(),"<ul><li><i>nstep=%d, n=%d, time=%f, dt=%g</i></li>\n",T->step(),n,T->time() - T->dt(),dt);
	fprintf(files->logFile(),"<li><i>Epot=%g, Ekin=%g, U=%g, E=%g</i></li></ul>\n",energy_tool->potentialEnergy(),ekin,eint,etot);
	fflush(files->logFile());
}

void CalcServer::printEnergy()
{
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *files = InputOutput::FileManager::singleton();
	energy();

	fprintf(files->enFile(),"%g\t",T->time() - T->dt());
	fprintf(files->enFile(),"%g\t%g\t%g\t%g\t%g\t%g\n", energy_tool->potentialEnergy(),
                                                        energy_tool->kineticEnergy(),
                                                        energy_tool->internalEnergy(),
                                                        energy_tool->enthalpy(),
                                                        energy_tool->entropy(),
                                                        energy_tool->energy());
	fflush(files->enFile());
}

void CalcServer::printBounds()
{
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *files = InputOutput::FileManager::singleton();

	bounds();

	fprintf(files->boundsFile(),"%g\t",T->time() - T->dt());
	fprintf(files->boundsFile(),"%g\t%g\t",min_fluid_bound.x,min_fluid_bound.y);
	#ifdef HAVE_3D
	    fprintf(files->boundsFile(),"%g\t",min_fluid_bound.z);
	#endif
	fprintf(files->boundsFile(),"%g\t%g\t",max_fluid_bound.x,max_fluid_bound.y);
	#ifdef HAVE_3D
	    fprintf(files->boundsFile(),"%g\t",max_fluid_bound.z);
	#endif
	fprintf(files->boundsFile(),"%g\t%g\n",min_v, max_v);
	fflush(files->boundsFile());
}

void CalcServer::energy()
{
	if(energy_computed)
	    return;
    if(energy_tool->execute())
        return;
    eint = energy_tool->internalEnergy();
    ekin = energy_tool->kineticEnergy();
    etot = energy_tool->energy();
	energy_computed = true;
}

void CalcServer::bounds()
{
	if(bounds_computed)
	    return;
    if(bounds_tool->execute())
        return;
    min_fluid_bound = bounds_tool->minCoords();
    max_fluid_bound = bounds_tool->maxCoords();
    min_v = length(bounds_tool->minVel());
    max_v = length(bounds_tool->maxVel());
	bounds_computed = true;
}

}}  // namespace
