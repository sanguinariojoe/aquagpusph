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
	: clPlatforms(NULL)
	, clDevices(NULL)
	, clComQueues(NULL)
	, energyPerformed(false)
	, boundsPerformed(false)
	, fluidMass(0.f)
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	unsigned int i;
	int Error=0;
	verbose_level = P->settings.verbose_level;
	// Create the OpenCL working environment
	clDeviceType = CL_DEVICE_TYPE_ALL;
	if(setupOpenCL()) {
	    exit(1);
	}
	// Allocate memory for the data in the computational device
	nfluid   = P->n_fluids;
	n        = 0;
	nSensors = P->SensorsParameters.pos.size();
	for(i=0;i<P->n_fluids;i++) {
	    n += P->fluids[i].n;
	}
	N = n + nSensors;
	AllocatedMem = 0;
	Error |= allocMemory(&imove,   N * sizeof( cl_int ));
	Error |= allocMemory(&imovein, N * sizeof( cl_int ));
	Error |= allocMemory(&ifluid,  N * sizeof( cl_int ));
	Error |= allocMemory(&ifluidin,N * sizeof( cl_int ));
	Error |= allocMemory(&pos,     N * sizeof( vec ));
	Error |= allocMemory(&normal,  N * sizeof( vec ));
	Error |= allocMemory(&v,       N * sizeof( vec ));
	Error |= allocMemory(&f,       N * sizeof( vec ));
	Error |= allocMemory(&dens,    N * sizeof( cl_float ));
	Error |= allocMemory(&drdt,    N * sizeof( cl_float ));
	Error |= allocMemory(&drdt_F,  N * sizeof( cl_float ));
	Error |= allocMemory(&mass,    N * sizeof( cl_float ));
	Error |= allocMemory(&hp,      N * sizeof( cl_float ));
	Error |= allocMemory(&press,   N * sizeof( cl_float ));
	Error |= allocMemory(&posin,   N * sizeof( vec ));
	Error |= allocMemory(&normalin,N * sizeof( vec ));
	Error |= allocMemory(&vin,     N * sizeof( vec ));
	Error |= allocMemory(&fin,     N * sizeof( vec ));
	Error |= allocMemory(&densin,  N * sizeof( cl_float ));
	Error |= allocMemory(&drdtin,  N * sizeof( cl_float ));
	Error |= allocMemory(&massin,  N * sizeof( cl_float ));
	Error |= allocMemory(&hpin,    N * sizeof( cl_float ));
	Error |= allocMemory(&pressin, N * sizeof( cl_float ));
	Error |= allocMemory(&sigma,   N * sizeof( cl_float ));
	Error |= allocMemory(&dtconv,  N * sizeof( cl_float ));
	Error |= allocMemory(&shepard, N * sizeof( cl_float ));
	Error |= allocMemory(&gradShepard, N * sizeof( vec ));
	nLcell = nextPowerOf2(N);
	nLcell = roundUp(nLcell, _ITEMS*_GROUPS);
	Error |= allocMemory(&permutation,        nLcell * sizeof( cl_uint ));
	Error |= allocMemory(&reversePermutation, nLcell * sizeof( cl_uint ));
	Error |= allocMemory(&lcell,              nLcell * sizeof( cl_uint ));
	ihoc = 0;        // ihoc must be allocated later, when we have a tentative of number of cells
	isValidCell = 0; // isValidCell must be allocated later, when we have a tentative of number of cells
	Error |= allocMemory(&gamma,   nfluid * sizeof( cl_float ));
	Error |= allocMemory(&refd,    nfluid * sizeof( cl_float ));
	Error |= allocMemory(&visc_dyn, nfluid * sizeof( cl_float ));
	Error |= allocMemory(&visc_kin, nfluid * sizeof( cl_float ));
	Error |= allocMemory(&visc_dyn_corrected, nfluid * sizeof( cl_float ));
	Error |= allocMemory(&delta, nfluid * sizeof( cl_float ));
	DT = 0;     // Minimum time step must be allocated later, into the time step reduction phase.
	Error |= allocMemory(&sensorMode, max(nSensors,(uint)1) * sizeof( cl_ushort ));
	if(Error) {
	    exit(2);
	}
	printf("\tINFO (CalcServer::Init): Allocated memory = %u bytes\n", (unsigned int)AllocatedMem);
	// Create the computation tools
	mPredictor     = new Predictor();
	mGrid          = new Grid();
	mLinkList      = new LinkList();
	mRates         = new Rates();
	mElasticBounce = new Boundary::ElasticBounce();
	mDeLeffe       = new Boundary::DeLeffe();
	mGhost         = new Boundary::GhostParticles();
	mShepard       = new Shepard();
	mCorrector     = new Corrector();
	mDomain        = new Domain();
	mTimeStep      = new TimeStep();
	mDensInt       = new DensityInterpolation();
	mSensors       = new Sensors();
	mEnergy        = new Energy();
	mBounds        = new Bounds();
	mMoves.clear();
	for(i=0;i<P->MoveParameters.size();i++){
	    if(P->MoveParameters.at(i)->MoveType == 0){
	        Movement::Quaternion *Move = new Movement::Quaternion();
	        mMoves.push_back(Move);
	    }
	    if(P->MoveParameters.at(i)->MoveType == 1){
	        Movement::LIQuaternion *Move = new Movement::LIQuaternion();
	        mMoves.push_back(Move);
	    }
	    if(P->MoveParameters.at(i)->MoveType == 2){
	        Movement::C1Quaternion *Move = new Movement::C1Quaternion();
	        mMoves.push_back(Move);
	    }
	    if(P->MoveParameters.at(i)->MoveType == 3){
	        Movement::ScriptQuaternion *Move = new Movement::ScriptQuaternion();
	        mMoves.push_back(Move);
	    }
	}
	mPortals.clear();
	for(i=0;i<P->Portals.size();i++){
	    Portal::Portal *portal = new Portal::Portal(P->Portals.at(i));
	    mPortals.push_back(portal);
	}
}

CalcServer::~CalcServer()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the predictor manager...\n");
	delete mPredictor; mPredictor=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the grid manager...\n");
	delete mGrid; mGrid=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the link list manager...\n");
	delete mLinkList; mLinkList=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the rates manager...\n");
	delete mRates; mRates=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the elasticBounce boundary condition manager...\n");
	delete mElasticBounce; mElasticBounce=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the DeLeffe boundary condition manager...\n");
	delete mDeLeffe; mDeLeffe=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the ghost particles manager...\n");
	delete mGhost; mGhost=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the 0th order correction manager...\n");
	delete mShepard; mShepard=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the corrector manager...\n");
	delete mCorrector; mCorrector=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the domain bounds manager...\n");
	delete mDomain; mDomain=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the time step manager...\n");
	delete mTimeStep; mTimeStep=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the density interpolator...\n");
	delete mDensInt; mDensInt=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the sensors manager...\n");
	delete mSensors; mSensors=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the energy computation tool...\n");
	delete mEnergy; mEnergy=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the bounds computation tool...\n");
	delete mBounds; mBounds=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Destroying the movement managers...\n");
	for(i=0;i<mMoves.size();i++) {
	    delete mMoves.at(i);
	}
	mMoves.clear();
    S->addMessage(1, "(CalcServer::~CalcServer): Sutting down the OpenCL context...\n");
	if(clContext)clReleaseContext(clContext);
	for(i=0;i<clNDevices;i++) {
	    if(clComQueues[i])clReleaseCommandQueue(clComQueues[i]);
	}
    S->addMessage(1, "(CalcServer::~CalcServer): Deallocating memory from the devices...\n");
	if(imove)clReleaseMemObject(imove); imove=0;
	if(imovein)clReleaseMemObject(imovein); imovein=0;
	if(ifluid)clReleaseMemObject(ifluid); ifluid=0;
	if(ifluidin)clReleaseMemObject(ifluidin); ifluidin=0;
	if(pos)clReleaseMemObject(pos); pos=0;
	if(normal)clReleaseMemObject(normal); normal=0;
	if(v)clReleaseMemObject(v); v=0;
	if(f)clReleaseMemObject(f); f=0;
	if(dens)clReleaseMemObject(dens); dens=0;
	if(drdt)clReleaseMemObject(drdt); drdt=0;
	if(drdt_F)clReleaseMemObject(drdt_F); drdt_F=0;
	if(mass)clReleaseMemObject(mass); mass=0;
	if(hp)clReleaseMemObject(hp); hp=0;
	if(press)clReleaseMemObject(press); press=0;
	if(posin)clReleaseMemObject(posin); posin=0;
	if(normalin)clReleaseMemObject(normalin); normalin=0;
	if(vin)clReleaseMemObject(vin); vin=0;
	if(fin)clReleaseMemObject(fin); fin=0;
	if(densin)clReleaseMemObject(densin); densin=0;
	if(drdtin)clReleaseMemObject(drdtin); drdtin=0;
	if(massin)clReleaseMemObject(massin); massin=0;
	if(hpin)clReleaseMemObject(hpin); hpin=0;
	if(pressin)clReleaseMemObject(pressin); pressin=0;
	if(sigma)clReleaseMemObject(sigma); sigma=0;
	if(dtconv)clReleaseMemObject(dtconv); dtconv=0;
	if(shepard)clReleaseMemObject(shepard); shepard=0;
	if(gradShepard)clReleaseMemObject(gradShepard); gradShepard=0;
	if(permutation)clReleaseMemObject(permutation); permutation=0;
	if(reversePermutation)clReleaseMemObject(reversePermutation); reversePermutation=0;
	if(lcell)clReleaseMemObject(lcell); lcell=0;
	if(ihoc)clReleaseMemObject(ihoc); ihoc=0;
	if(isValidCell)clReleaseMemObject(isValidCell); isValidCell=0;
	if(gamma)clReleaseMemObject(gamma); gamma=0;
	if(refd)clReleaseMemObject(refd); refd=0;
	if(visc_dyn)clReleaseMemObject(visc_dyn); visc_dyn=0;
	if(visc_kin)clReleaseMemObject(visc_kin); visc_kin=0;
	if(visc_dyn_corrected)clReleaseMemObject(visc_dyn_corrected); visc_dyn_corrected=0;
	if(delta)clReleaseMemObject(delta); delta=0;
	if(DT)clReleaseMemObject(DT); DT=0;
	if(sensorMode)clReleaseMemObject(sensorMode); sensorMode=0;
    S->addMessage(1, "(CalcServer::~CalcServer): Deallocating host memory...\n");
	if(clPlatforms) delete[] clPlatforms; clPlatforms=0;
	if(clDevices) delete[] clDevices; clDevices=0;
	if(clComQueues) delete[] clComQueues; clComQueues=0;
}

bool CalcServer::update()
{
	InputOutput::TimeManager *T   = InputOutput::TimeManager::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	unsigned int i;
	while(!T->mustPrintOutput() && !T->mustStop())
	{
	    if(mPredictor->execute())
	        return true;
	    if(mLLStep >= link_list_steps){
	        mLLStep = 0;
	        if(mGrid->execute())
	            return true;
	        if(mLinkList->execute())
	            return true;
	    }
	    mLLStep++;
	    if(mRates->execute())
	        return true;
        // Since the density interpolation will not take into account the
        // boundaries, we must perform it before the shepard term is corrected.
	    if(T->time() >= 0.f){
            if(mDensStep >= dens_int_steps){
                mDensStep = 0;
                if(mDensInt->execute())
                    return true;
            }
            mDensStep++;
	    }
	    if(mDeLeffe->execute())
	        return true;
	    if(mGhost->execute())
	        return true;
	    if(mShepard->execute())
	        return true;
	    if(mElasticBounce->execute())
	        return true;
	    if(mCorrector->execute())
	        return true;
	    for(i=0;i<mPortals.size();i++){
	        if(mPortals.at(i)->execute())
	            return true;
	    }
	    if(mDomain->execute())
	        return true;
	    if(mTimeStep->execute())
	        return true;
	    if(T->time() >= 0.f){
	        for(i=0;i<mMoves.size();i++){
	            if(mMoves.at(i)->execute())
	                return true;
	        }
	        if(mSensors->execute())
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
		energyPerformed=false;
		boundsPerformed=false;
		// Key events
		while(isKeyPressed()){
	        if(getchar() == 'c'){
	            S->addMessage(1, (char *)"(CalcServer::Update): Interrumption request detected.\n");
	            return true;
	        }
		}
	}
	return false;
}

bool CalcServer::getData(void *Dest, cl_mem Orig, size_t Size, size_t Offset)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	cl_int clFlag;
	clFlag  = clEnqueueReadBuffer(clComQueue, Orig, CL_TRUE, Offset, Size, Dest, 0, NULL, NULL);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getData): Failure retrieving memory from the server.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	return false;
}

bool CalcServer::sendData(cl_mem Dest, void* Orig, size_t Size)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	cl_int clFlag;
	clFlag  = clEnqueueWriteBuffer(clComQueue, Dest, CL_TRUE, 0, Size, Orig, 0, NULL, NULL);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::sendData): Failure sending memory to the server.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	return false;
}

bool CalcServer::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "(CalcServer::setupOpenCL): Initializating OpenCL...\n");
	if(queryOpenCL()){
	    return true;
	}
	if(getPlatform()){
	    return true;
	}
	if(getDevices()){
	    if(clPlatforms) delete[] clPlatforms; clPlatforms=0;
	    return true;
	}
	S->addMessage(1, "(CalcServer::setupOpenCL): OpenCL is ready to work!\n");
	return false;
}

bool CalcServer::queryOpenCL()
{
	cl_device_id *devices;
	cl_uint i,j,nDevices=0;
	cl_int clFlag;
	char msg[1024], aux[1024];
	clPlatforms = NULL;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	strcpy(msg, "");
	// Gets the total number of platforms
	clFlag = clGetPlatformIDs(0, NULL, &clNPlatforms);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::queryOpenCL): Can't take the number of platforms.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	// Get array of patforms
	clPlatforms = new cl_platform_id[clNPlatforms];
	if(!clPlatforms) {
	    S->addMessage(3, "(CalcServer::queryOpenCL): Allocation memory error.\n");
	    S->addMessage(0, "\tPlatforms array can't be allocated\n");
	    return true;
	}
	clFlag = clGetPlatformIDs(clNPlatforms, clPlatforms, NULL);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::queryOpenCL): Can't take the platforms list.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	for(i=0;i<clNPlatforms;i++){
	    // Get the number of devices
	    clFlag = clGetDeviceIDs (clPlatforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &nDevices);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(CalcServer::queryOpenCL): Can't take the number of devices.\n");
            S->printOpenCLError(clFlag);
	        return true;
	    }
	    // Gets the devices array
	    devices = new cl_device_id[nDevices];
	    if(!devices) {
	        S->addMessage(3, "(CalcServer::queryOpenCL): Allocation memory error.\n");
	        S->addMessage(0, "\tDevices array can't be allocated\n");
	        return true;
	    }
	    clFlag = clGetDeviceIDs(clPlatforms[i], CL_DEVICE_TYPE_ALL, nDevices, devices, &nDevices);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(CalcServer::queryOpenCL): Can't take the devices list.\n");
            S->printOpenCLError(clFlag);
	        return true;
	    }
	    // Shows device arrays
	    for(j=0;j<nDevices;j++){
	        // Identifier
	        strcpy(msg, "");
	        sprintf(msg, "\tDevice %u, Platform %u...\n", j, i);
	        S->addMessage(1, msg);
	        // Device name
	        clFlag = clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 1024*sizeof(char), &aux, NULL);
	        if(clFlag != CL_SUCCESS) {
	            S->addMessage(3, "(CalcServer::queryOpenCL): Can't get the device name.\n");
                S->printOpenCLError(clFlag);
	            return true;
	        }
	        strcpy(msg, "");
	        sprintf(msg, "\t\tDEVICE: %s\n", aux);
	        S->addMessage(0, msg);
	        // Platform vendor
	        clFlag = clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, 1024*sizeof(char), &aux, NULL);
	        if(clFlag != CL_SUCCESS) {
	            S->addMessage(3, "(CalcServer::queryOpenCL): Can't get the device vendor.\n");
                S->printOpenCLError(clFlag);
	            return true;
	        }
	        strcpy(msg, "");
	        sprintf(msg, "\t\tVENDOR: %s\n", aux);
	        S->addMessage(0, msg);
	        // Device type
	        cl_device_type dType;
	        clFlag = clGetDeviceInfo(devices[j], CL_DEVICE_TYPE, sizeof(cl_device_type), &dType, NULL);
	        if(clFlag != CL_SUCCESS) {
	            S->addMessage(3, "(CalcServer::queryOpenCL): Can't get the device type.\n");
                S->printOpenCLError(clFlag);
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
	if(P->settings.platform_id >= clNPlatforms){
	    S->addMessage(3, "(CalcServer::getPlatform): Can't find the requested platform.\n");
	    strcpy(msg, "");
	    sprintf(msg, "\t%u platform requested, but just %u platforms can be found\n", P->settings.platform_id, clNPlatforms);
	    S->addMessage(0, msg);
	    return true;
	}
	clPlatform = clPlatforms[P->settings.platform_id];
	return false;
}

bool CalcServer::getDevices()
{
	cl_int clFlag;
	cl_uint i;
	char msg[1024];
	InputOutput::ProblemSetup  *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	clDevices = NULL;
	// Gets the number of valid devices
	clFlag = clGetDeviceIDs (clPlatform, P->settings.device_type, 0, NULL, &clNDevices);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't take the number of devices.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	if(P->settings.device_id >= clNDevices) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't find the requested device.\n");
	    strcpy(msg, "");
	    sprintf(msg, "\t%u device requested, but just %u devices have been found\n", P->settings.device_id, clNDevices);
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
	clDevices = new cl_device_id[clNDevices];
	if(!clDevices) {
	    S->addMessage(3, "(CalcServer::getDevices): Allocation memory error.\n");
	    S->addMessage(0, "\tDevices array can't be allocated\n");
	    return true;
	}
	clFlag = clGetDeviceIDs(clPlatform, P->settings.device_type, clNDevices, clDevices, &clNDevices);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't take the devices list.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	// Create a devices context
	clContext = clCreateContext(0, clNDevices, clDevices, NULL, NULL, &clFlag);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::getDevices): Can't create an OpenCL context.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	// Create command queues
	clComQueues = new cl_command_queue[clNDevices];
	if(clComQueues == NULL) {
	    S->addMessage(3, "(CalcServer::getDevices): Allocation memory error.\n");
	    S->addMessage(0, "\tCommand queues array can't be allocated\n");
	    return true;
	}
	for(i=0;i<clNDevices;i++) {
	    #ifdef HAVE_GPUPROFILE
	        clComQueues[i] = clCreateCommandQueue(clContext, clDevices[i], CL_QUEUE_PROFILING_ENABLE, &clFlag);
	    #else
	        clComQueues[i] = clCreateCommandQueue(clContext, clDevices[i], 0, &clFlag);
	    #endif
	    if(clFlag != CL_SUCCESS) {
	        strcpy(msg, "");
	        sprintf(msg, "(CalcServer::getDevices): Can't create a command queue for the device %u.\n",i);
	        S->addMessage(3, msg);
            S->printOpenCLError(clFlag);
	        return true;
	    }
	}
	// Store the selected ones
	clDevice   = clDevices[P->settings.device_id];
	clComQueue = clComQueues[P->settings.device_id];
	return false;
}

bool CalcServer::allocMemory(cl_mem *clID, size_t size)
{
	int clFlag;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	*clID = clCreateBuffer(clContext, CL_MEM_READ_WRITE, size, NULL, &clFlag);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::allocMemory): Allocation failure.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}

	AllocatedMem += size;
	return false;
}

bool CalcServer::fillFluid()
{
	unsigned int i;
	InputOutput::Fluid *F = InputOutput::Fluid::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	char *file_name = new char[256];
	int clFlag=0;
	cl_program clProgram;
	cl_kernel clKernel;
	size_t clGlobalWorkSize;
	size_t clLocalWorkSize;
	unsigned int AuxN=0;

	for(i=0;i<nfluid;i++) {
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
	    if(!loadKernelFromFile(&clKernel, &clProgram, clContext, clDevice, file_name, "ParticlesDistribution", ""))
	        return true;
	    clFlag |= sendArgument(clKernel,  0, sizeof(cl_mem ), (void*)&imove);
	    clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem ), (void*)&ifluid);
	    clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem ), (void*)&pos);
	    clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem ), (void*)&normal);
	    clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem ), (void*)&v);
	    clFlag |= sendArgument(clKernel,  5, sizeof(cl_mem ), (void*)&f);
	    clFlag |= sendArgument(clKernel,  6, sizeof(cl_mem ), (void*)&dens);
	    clFlag |= sendArgument(clKernel,  7, sizeof(cl_mem ), (void*)&drdt);
	    clFlag |= sendArgument(clKernel,  8, sizeof(cl_mem ), (void*)&mass);
	    clFlag |= sendArgument(clKernel,  9, sizeof(cl_mem ), (void*)&hp);
	    clFlag |= sendArgument(clKernel, 10, sizeof(cl_uint), (void*)&AuxN);
	    clFlag |= sendArgument(clKernel, 11, sizeof(cl_uint), (void*)&(P->fluids[i].n));
	    clFlag |= sendArgument(clKernel, 12, sizeof(cl_uint), (void*)&i);
	    clFlag |= sendArgument(clKernel, 13, sizeof(vec    ), (void*)&(P->SPH_opts.deltar));
	    clFlag |= sendArgument(clKernel, 14, sizeof(cl_float), (void*)&(P->fluids[i].refd));
	    clFlag |= sendArgument(clKernel, 15, sizeof(cl_float), (void*)&(P->SPH_opts.h));
	    clFlag |= sendArgument(clKernel, 16, sizeof(cl_float), (void*)&(P->SPH_opts.cs));
	    AuxN += P->fluids[i].n;
	    if(clFlag)
	        return true;
	    clLocalWorkSize = 256;
	    clGlobalWorkSize = getGlobalWorkSize(P->fluids[i].n, clLocalWorkSize);
	    sprintf(msg, "(CalcServer::fillFluid): Launching the initialization OpenCL script.\n");
	    S->addMessage(1, msg);
	    sprintf(msg, "\tLocal work size = %lu.\n", clLocalWorkSize);
	    S->addMessage(0, msg);
	    sprintf(msg, "\tGlobal work size = %lu.\n", clGlobalWorkSize);
	    S->addMessage(0, msg);
	    clFlag = clEnqueueNDRangeKernel(clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(CalcServer::fillFluid): Can't execute the kernel.\n");
            S->printOpenCLError(clFlag);
	        return true;
	    }
	    if(clKernel)clReleaseKernel(clKernel);
	    if(clProgram)clReleaseProgram(clProgram);
	}
	S->addMessage(1, "(CalcServer::fillFluid): Finished the work, reading back the info...\n");

	clFlag |= getData((void*)F->imove, imove, sizeof(cl_int)*n);
	clFlag |= getData((void*)F->ifluid, ifluid, sizeof(cl_int)*n);
	clFlag |= getData((void*)F->pos, pos, sizeof(vec)*n);
	clFlag |= getData((void*)F->normal, normal, sizeof(vec)*n);
	clFlag |= getData((void*)F->v, v, sizeof(vec)*n);
	clFlag |= getData((void*)F->f, f, sizeof(vec)*n);
	clFlag |= getData((void*)F->dens, dens, sizeof(cl_float)*n);
	clFlag |= getData((void*)F->drdt, drdt, sizeof(cl_float)*n);
	clFlag |= getData((void*)F->mass, mass, sizeof(cl_float)*n);
	clFlag |= getData((void*)F->hp, hp, sizeof(cl_float)*n);
	if(clFlag != CL_SUCCESS) {
	    S->addMessage(3, "(CalcServer::fillFluid): Failure retrieving memory from the server.\n");
	    S->printOpenCLError(clFlag);
	    return true;
	}
	sprintf(msg, "(CalcServer::fillFluid): %u particles have been transfered from the server to the host.\n", n);
	S->addMessage(1,msg);
	S->addMessage(1, "(CalcServer::fillFluid): All the fluids are running! ;-)\n");

	delete[] file_name;
	return false;
}

bool CalcServer::setup()
{
	unsigned int i;
	cl_uint clFlag=0;
	InputOutput::Fluid *F = InputOutput::Fluid::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	// Fluids data
    S->addMessage(1, "(CalcServer::setup): Sending fluids data to the server...\n");
	cl_float *Gamma   = new cl_float[nfluid];
	cl_float *Refd    = new cl_float[nfluid];
	cl_float *VIscdyn = new cl_float[nfluid];
	cl_float *VIsckin = new cl_float[nfluid];
	cl_float *VIscdynCorr = new cl_float[nfluid];
	cl_float *Delta = new cl_float[nfluid];
	for(i=0;i<nfluid;i++) {
		Gamma[i] = P->fluids[i].gamma;
		Refd[i] = P->fluids[i].refd;
		VIscdyn[i] = P->fluids[i].visc_dyn;
		VIsckin[i] = P->fluids[i].visc_kin;
		VIscdynCorr[i] = P->fluids[i].visc_dyn_corrected;
		Delta[i] = P->fluids[i].delta;
	}
	clFlag |= sendData(gamma, Gamma, sizeof(cl_float)*nfluid);
	clFlag |= sendData(refd, Refd, sizeof(cl_float)*nfluid);
	clFlag |= sendData(visc_dyn, VIscdyn, sizeof(cl_float)*nfluid);
	clFlag |= sendData(visc_kin, VIsckin, sizeof(cl_float)*nfluid);
	clFlag |= sendData(visc_dyn_corrected, VIscdynCorr, sizeof(cl_float)*nfluid);
	clFlag |= sendData(delta, Delta, sizeof(cl_float)*nfluid);
	delete[] Gamma;Gamma=0;
	delete[] Refd;Refd=0;
	delete[] VIscdyn;VIscdyn=0;
	delete[] VIsckin;VIsckin=0;
	delete[] VIscdynCorr;VIscdynCorr=0;
	delete[] Delta;Delta=0;
	if(clFlag != CL_SUCCESS) {
        S->addMessage(3, "(CalcServer::setup): Can't send the fluids data to the server.\n");
	    return true;
	}
	// Variables
	g     = P->SPH_opts.g;
	hfac  = P->SPH_opts.hfac;
	dt_divisor = P->SPH_opts.dt_divisor;
	h     = P->SPH_opts.h;
	cs    = P->SPH_opts.cs;
	link_list_steps = P->SPH_opts.link_list_steps;
	dens_int_steps = P->SPH_opts.dens_int_steps;
	// Calculate CellFac
	float sep;
	#ifdef __CUBIC_KERNEL_TYPE__
		sep = 2.0f;
	#elif defined(__GAUSS_KERNEL_TYPE__)
		sep = 3.0f;
	#else   // Wendland
		sep = 2.0f;
	#endif
	float dist  = sep*h;        // Minimum cell size.
	float ddt   = 0.1f*h/dt_divisor; // Maximum distance that a particle can move in a step.
	CellFac = 1.f + (link_list_steps-1) * ddt / dist;
	sprintf(msg, "(CalcServer::setup): Cells size increased with %g factor.\n", CellFac);
    S->addMessage(1, msg);
	mLLStep = link_list_steps;          // Initializating as LinkList needed
	mDensStep = 0;              // Initializating as Density interpolation not needed
	// Particles data
    S->addMessage(1, "(CalcServer::setup): Sending the particles data to server...\n");
	clFlag  = sendData(imove, F->imove, sizeof(cl_int)*N);
	clFlag |= sendData(imovein, F->imove, sizeof(cl_int)*N);
	clFlag |= sendData(ifluid, F->ifluid, sizeof(cl_int)*N);
	clFlag |= sendData(ifluidin, F->ifluid, sizeof(cl_int)*N);
	clFlag |= sendData(posin, F->pos, sizeof(vec)*N);
	clFlag |= sendData(normalin, F->normal, sizeof(vec)*N);
	clFlag |= sendData(vin, F->v, sizeof(vec)*N);
	clFlag |= sendData(densin, F->dens, sizeof(cl_float)*N);
	clFlag |= sendData(hpin, F->hp, sizeof(cl_float)*N);
	clFlag |= sendData(drdtin, F->drdt, sizeof(cl_float)*N);
	clFlag |= sendData(fin, F->f, sizeof(vec)*N);
	clFlag |= sendData(massin, F->mass, sizeof(cl_float)*N);
	clFlag |= sendData(pressin, F->press, sizeof(cl_float)*N);
	clFlag |= sendData(pos, F->pos, sizeof(vec)*N);
	clFlag |= sendData(normal, F->normal, sizeof(vec)*N);
	clFlag |= sendData(v, F->v, sizeof(vec)*N);
	clFlag |= sendData(dens, F->dens, sizeof(cl_float)*N);
	clFlag |= sendData(drdt, F->drdt, sizeof(cl_float)*N);
	clFlag |= sendData(hp, F->hp, sizeof(cl_float)*N);
	clFlag |= sendData(f, F->f, sizeof(vec)*N);
	clFlag |= sendData(mass, F->mass, sizeof(cl_float)*N);
	clFlag |= sendData(press, F->press, sizeof(cl_float)*N);
	cl_float *Sigma = new cl_float[N];
	for(i=0;i<N;i++)
	    Sigma[i] = 1000.f;
	clFlag |= sendData(sigma, Sigma, N*sizeof(cl_float));
	delete[] Sigma; Sigma=0;
	if(clFlag != CL_SUCCESS) {
        S->addMessage(3, "(CalcServer::setup): Can't send the particles data to the server.\n");
	    return true;
	}
	cl_ushort *mode = new cl_ushort[nSensors];
	for(i=0;i<N-n;i++)
	    mode[i] = P->SensorsParameters.mod.at(i);
	clFlag |= sendData(sensorMode, mode, sizeof(cl_ushort)*(N-n));
	delete[] mode; mode = 0;

	// Inital values
	lxy    = 0;
	lxydim = 0;
	dt     = 0;
	// Movements
	for(i=0;i<mMoves.size();i++){
	    if(mMoves.at(i)->parse(P->MoveParameters.at(i)->defFile))
	        return true;
	}
	// Compute the total fluid mass
	for(i=0;i<N;i++) {
	    if(F->imove[i] > 0)
	        fluidMass += F->mass[i];
	}
    S->addMessage(1, "(CalcServer::setup): Calculation server is ready! ;-) \n");
	return false;
}

void CalcServer::printLog()
{
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *Files = InputOutput::FileManager::singleton();
	energy();
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	char date[64], *aux;
	strcpy(date, ctime(&seconds));
	aux = strstr(date, "\n");
	if(aux) strcpy(aux, "");
	fprintf(Files->logFile(),"<i>%s: Printed file (%d)</i><br>\n",date,T->frame());
	fprintf(Files->logFile(),"<ul><li><i>nstep=%d, n=%d, time=%f, dt=%g</i></li>\n",T->step(),n,T->time() - T->dt(),dt);
	fprintf(Files->logFile(),"<li><i>Epot=%g, Ekin=%g, U=%g, E=%g</i></li></ul>\n",mEnergy->potentialEnergy(),ekin,eint,etot);
	fflush(Files->logFile());
}

void CalcServer::printEnergy()
{
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *Files = InputOutput::FileManager::singleton();
	energy();
	fprintf(Files->enFile(),"%g\t",T->time() - T->dt());
	fprintf(Files->enFile(),"%g\t%g\t%g\t%g\t%g\t%g\n", mEnergy->potentialEnergy(),
                                                        mEnergy->kineticEnergy(),
                                                        mEnergy->internalEnergy(),
                                                        mEnergy->enthalpy(),
                                                        mEnergy->entropy(),
                                                        mEnergy->energy());
	fflush(Files->enFile());
}

void CalcServer::printBounds()
{
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::FileManager *Files = InputOutput::FileManager::singleton();
	bounds();
	fprintf(Files->boundsFile(),"%g\t",T->time() - T->dt());
	fprintf(Files->boundsFile(),"%g\t%g\t",minCoords.x,minCoords.y);
	#ifdef HAVE_3D
	    fprintf(Files->boundsFile(),"%g\t",minCoords.z);
	#endif
	fprintf(Files->boundsFile(),"%g\t%g\t",maxCoords.x,maxCoords.y);
	#ifdef HAVE_3D
	    fprintf(Files->boundsFile(),"%g\t",maxCoords.z);
	#endif
	fprintf(Files->boundsFile(),"%g\t%g\n",minV, maxV);
	fflush(Files->boundsFile());
}

void CalcServer::energy()
{
	if(energyPerformed)
	    return;
    if(mEnergy->execute())
        return;
    eint = mEnergy->internalEnergy();
    ekin = mEnergy->kineticEnergy();
    etot = mEnergy->energy();
	energyPerformed = true;
}

void CalcServer::bounds()
{
	if(boundsPerformed)
	    return;
    if(mBounds->execute())
        return;
    minCoords = mBounds->minCoords();
    maxCoords = mBounds->maxCoords();
    minV = length(mBounds->minVel());
    maxV = length(mBounds->maxVel());
	boundsPerformed = true;
}

}}  // namespace
