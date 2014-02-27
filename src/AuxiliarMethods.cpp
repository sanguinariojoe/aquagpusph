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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <string>

#include <AuxiliarMethods.h>
#include <ProblemSetup.h>
#include <ScreenManager.h>

namespace Aqua{

int isKeyPressed()
{
	struct termios oldt, newt;
	int ch;
	int oldf;

	tcgetattr(STDIN_FILENO, &oldt);
	newt = oldt;
	newt.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(STDIN_FILENO, TCSANOW, &newt);
	oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
	fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

	ch = getchar();

	tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
	fcntl(STDIN_FILENO, F_SETFL, oldf);

	if(ch != EOF) {
	    ungetc(ch, stdin);
	    return 1;
	}

	return 0;
}

unsigned int nextPowerOf2( unsigned int n )
{
	if (n & !(n & (n - 1)))
        return n;

    unsigned int p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

unsigned int isPowerOf2( unsigned int n )
{
	return ((n&(n-1))==0);
}

unsigned int roundUp(unsigned int n, unsigned int divisor)
{
	unsigned int rest = n%divisor;
	if(rest) {
	    n -= rest;
	    n += divisor;
	}
	return n;
}

size_t loadKernelFromFile(cl_kernel* kernel, cl_program* program,
                          cl_context context, cl_device_id device,
                          const char* path, const char* entry_point,
                          const char* flags, const char* header)
{
	char* source = NULL;
	char* default_flags = NULL;
	size_t source_length = 0;
	int err_code = CL_SUCCESS;
	*program = NULL;
	*kernel = NULL;
	size_t work_group_size=0;
	InputOutput::ProblemSetup* problemSetup = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

	// Allocate the required memory for the source code
	source_length = readFile(NULL, path);
	if(!source_length)
	    return 0;
	source = new char[source_length+1];
	if(!source) {
	    S->addMessageF(3, "Imposible to allocate memory on host for the source code.\n");
	    return 0;
	}
    // Read the file
	source_length = readFile(source, path);
	if(!source_length){
	    delete[] source; source=NULL;
	    return 0;
	}
	// Append the header on top of the source code
	if(header){
        char *backup   = source;
        source_length += strlen(header)*sizeof(char);
        source        = new char[source_length+1];
        if(!source) {
            S->addMessageF(3, "Imposible to allocate memory on host to append header.\n");
            return 0;
        }
        strcpy(source, header);
        strcat(source, backup);
        delete[] backup;
	}

    // Create the program and compile it
	*program = clCreateProgramWithSource(context, 1, (const char **)&source,
                                         &source_length, &err_code);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "OpenCl source code send error.\n");
	    S->printOpenCLError(err_code);
	    delete[] source; source=NULL;
	    return 0;
	}
	default_flags = new char[1024];
	#ifdef AQUA_DEBUG
	    strcpy(default_flags, "-g -DDEBUG ");
	#else
	    strcpy(default_flags, "-DNDEBUG ");
	#endif
	strcat(default_flags, "-I");
	const char *folder = getFolderFromFilePath(path);
	strcat(default_flags, folder);

	strcat(default_flags, " -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-fast-relaxed-math ");
	#ifdef HAVE_3D
		strcat(default_flags, " -DHAVE_3D -Dvec=float4 -Divec=int4 -Duivec=uint4 -DVEC_ZERO=(float4)(0.f,0.f,0.f,0.f) ");
	#else
		strcat(default_flags, " -DHAVE_2D -Dvec=float2 -Divec=int2 -Duivec=uint2 -DVEC_ZERO=(float2)(0.f,0.f) ");
	#endif
	sprintf(default_flags, "%s -Dh=%ff ", default_flags,
            problemSetup->SPH_opts.h);
	sprintf(default_flags, "%s -D__BOUNDARY__=%u ", default_flags,
            problemSetup->SPH_opts.boundary_type);
	switch(problemSetup->SPH_opts.slip_condition){
	    case 1:
	        strcat(default_flags, "-D__FREE_SLIP__ ");
	        break;
	    default:
	        strcat(default_flags, "-D__NO_SLIP__ ");
	}

	strcat(default_flags, flags);
	char msg[1024];
	strcpy(msg, "\"");
	strcat(msg, path);
	strcat(msg, "\" :: \"");
	strcat(msg, entry_point);
	strcat(msg, "\"\n");
	S->addMessageF(1, msg);
	S->addMessage(0, "Flags = \"");
	S->addMessage(0, default_flags);
	S->addMessage(0, "\"\n");
	if(header){
        S->addMessageF(1, "Header append:\n");
        S->addMessage(0, header);
        S->addMessage(0, "\n");
	}
	err_code = clBuildProgram(*program, 0, NULL, default_flags, NULL, NULL);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't compile OpenCl source code.\n");
        S->printOpenCLError(err_code);
	    S->addMessage(3, "--- Build log ---------------------------------\n");
	    char log[10240];
	    clGetProgramBuildInfo(*program, device, CL_PROGRAM_BUILD_LOG,
                              10240*sizeof(char), log, NULL );
	    strcat(log, "\n");
	    S->addMessage(0, log);
	    S->addMessage(3, "--------------------------------- Build log ---\n");
	    delete[] source; source=NULL;
	    delete[] default_flags; default_flags=NULL;
	    return 0;
	}
	char log[10240];
	unsigned int log_start = 0;
	clGetProgramBuildInfo(*program, device, CL_PROGRAM_BUILD_LOG,
                          10240*sizeof(char), log, NULL );
	while(    (log[log_start] == ' ')
           || (log[log_start] == '\t')
           || (log[log_start] == '\n'))
    {
        log_start++;
    }
	if(strcmp(&(log[log_start]), "")){
	    S->addMessage(2, "--- Build log ---------------------------------\n");
	    strcat(&(log[log_start]), "\n");
	    S->addMessage(0, &(log[log_start]));
	    S->addMessage(2, "--------------------------------- Build log ---\n");
	}

	*kernel = clCreateKernel(*program, entry_point, &err_code);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't create kernel.\n");
	    S->printOpenCLError(err_code);
	    delete[] source; source=NULL;
	    delete[] default_flags; default_flags=NULL;
	    return 0;
	}
	delete[] source; source=NULL;
	delete[] default_flags; default_flags=NULL;
	err_code = clGetKernelWorkGroupInfo(*kernel, device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t), &work_group_size,
                                        NULL);
	if(err_code != CL_SUCCESS) {
	    S->addMessageF(3, "Can't get maximum work group size from kernel.\n");
        S->printOpenCLError(err_code);
	    return 0;
	}
	return work_group_size;
}

static char folder[1024];

const char* getFolderFromFilePath(const char* file_path)
{
	int i;
	strcpy(folder, "");
	if(file_path[0] != '/')
	    strcat(folder, "./");
	for(i=strlen(file_path);i>0;i--) {
	    if(file_path[i-1] == '/') {
	        break;
	    }
	}
	strncat(folder, file_path, i);
	return folder;
}

const char* getExtensionFromFilePath(const char* file_path)
{
	int i;
	for(i=strlen(file_path)-1;i>0;i--) {
	    if(file_path[i-1] == '.') {
	        break;
	    }
	}
	return strstr(file_path, &file_path[i]);
}

int isFile(const char* file_name)
{
	FILE *streamer=0;
	streamer = fopen(file_name, "rb");
	if(streamer){
	    fclose(streamer);
	    return 1;
	}
	return 0;
}

size_t readFile(char* source_code, const char* file_name)
{
	size_t length = 0;
	FILE *file = NULL;
    char msg[1024];
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

	file = fopen(file_name, "rb");
	if(file == NULL) {
	    sprintf(msg, "Imposible to open \"%s\".\n", file_name);
	    S->addMessageF(3, msg);
	    return 0;
	}
	fseek(file, 0, SEEK_END);
	length = ftell(file);
	fseek(file, 0, SEEK_SET);
	if(length == 0) {
	    strcpy(msg, "(ReadFile): \"");
	    strcat(msg,file_name);
	    strcat(msg,"\" is empty.\n");
	    S->addMessageF(3, msg);
	    fclose(file);
	    return 0;
	}
	if(!source_code) {
	    fclose(file);
	    return length;
	}
	size_t readed = fread(source_code, length, 1, file);
	source_code[length] = '\0';

	fclose(file);
	return length;
}

int sendArgument(cl_kernel kernel, int index, size_t size, void* ptr)
{
	int err_code;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	err_code = clSetKernelArg(kernel, index, size, ptr);
	if(err_code != CL_SUCCESS) {
	    char msg[1025]; strcpy(msg, "");
	    sprintf(msg, "The argument %d can not be set.\n", index);
	    S->addMessageF(3, msg);
	    S->printOpenCLError(err_code);
	    return 1;
	}
	return 0;
}

size_t getLocalWorkSize(cl_uint n, cl_command_queue queue)
{
	cl_int flag;
	cl_device_id d;
	flag = clGetCommandQueueInfo(queue,CL_QUEUE_DEVICE,
	                             sizeof(cl_device_id),&d, NULL);
	if(flag != CL_SUCCESS){
	    return 0;
	}
	// Start trying maximum local work size per dimension
	cl_uint dims;
	flag = clGetDeviceInfo(d,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
	                       sizeof(cl_uint),&dims, NULL);
	if(flag != CL_SUCCESS){
	    return 0;
	}
	size_t l[dims];
	flag = clGetDeviceInfo(d,CL_DEVICE_MAX_WORK_ITEM_SIZES,
	                       dims*sizeof(size_t),l, NULL);
	if(flag != CL_SUCCESS){
	    return 0;
	}
	// Correct it with maximum local size
	size_t max_l;
	flag = clGetDeviceInfo(d,CL_DEVICE_MAX_WORK_GROUP_SIZE,
	                       sizeof(size_t),&max_l, NULL);
	if(flag != CL_SUCCESS){
	    return 0;
	}
	if(max_l < l[0])
	    l[0] = max_l;
	return l[0];
}

size_t getGlobalWorkSize(cl_uint n, size_t local_work_size)
{
	return roundUp(n, local_work_size);
}

vec Vzero()
{
	vec r;
	r.x = 0.f; r.y = 0.f;
	#ifdef HAVE_3D
	    r.z = 0.f; r.w = 0.f;
	#endif
	return r;
}

vec Vx()
{
	vec r;
	r.x = 1.f; r.y = 0.f;
	#ifdef HAVE_3D
	    r.z = 0.f; r.w = 0.f;
	#endif
	return r;
}

vec Vy()
{
	vec r;
	r.x = 0.f; r.y = 1.f;
	#ifdef HAVE_3D
	    r.z = 0.f; r.w = 0.f;
	#endif
	return r;
}

#ifdef HAVE_3D
vec Vz()
{
	vec r;
	r.x = 0.f; r.y = 0.f;
	r.z = 1.f; r.w = 0.f;
	return r;
}
#endif

vec mult(float n, vec v)
{
	vec r;
	r.x = n*v.x;
	r.y = n*v.y;
	#ifdef HAVE_3D
	    r.z = n + v.z;
	    r.w = n + v.w;
	#endif
	return r;
}

vec add(vec a, vec b)
{
	vec r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	#ifdef HAVE_3D
	    r.z = a.z + b.z;
	    r.w = a.w + b.w;
	#endif
	return r;
}

vec sub(vec a, vec b)
{
	vec r;
	r.x = a.x - b.x;
	r.y = a.y - b.y;
	#ifdef HAVE_3D
	    r.z = a.z - b.z;
	    r.w = a.w - b.w;
	#endif
	return r;
}

float dot(vec a, vec b)
{
	float d = a.x*b.x + a.y*b.y;
	#ifdef HAVE_3D
	    d += a.z*b.z;
	    d += a.w*b.w;
	#endif
	return d;
}

float length(vec v)
{
	float m = v.x*v.x + v.y*v.y;
	#ifdef HAVE_3D
	    m += v.z*v.z;
	#endif
	return sqrt(m);
}

vec normalize(vec v)
{
    float m = length(v);
	vec n;
	n.x = v.x / m;
	n.y = v.y / m;
	#ifdef HAVE_3D
	    n.z = v.z / m;
	#endif
	return n;
}

#ifdef HAVE_3D
vec cross(vec a, vec b)
{
	vec c;
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;
	c.w = 0.f;
	return c;
}
#endif

unsigned int numberOfDigits(unsigned int number)
{
    return number > 0 ? (int)log10((double)number) + 1 : 1;
}

}   // namespace
