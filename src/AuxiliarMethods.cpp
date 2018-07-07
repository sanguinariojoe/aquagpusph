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
 * @brief Set of auxiliar functions.
 */

#include <algorithm> 
#include <cctype>
#include <locale>

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

bool hasSuffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void replaceAll(std::string &str,
                       const std::string &search,
                       const std::string &replace)
{
    size_t pos = 0;
    while((pos = str.find(search, pos)) != std::string::npos) {
        str.erase(pos, search.size());
        str.insert(pos, replace);        
    }
}

std::string replaceAllCopy(std::string str,
                                  std::string search,
                                  std::string replace)
{
    replaceAll(str, search, replace);
    return str;
}

void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

std::string ltrimCopy(std::string s) {
    ltrim(s);
    return s;
}

std::string rtrimCopy(std::string s) {
    rtrim(s);
    return s;
}

std::string trimCopy(std::string s) {
    trim(s);
    return s;
}

static std::string xxd_str;

std::string xxd2string(unsigned char* arr, unsigned int len)
{
    char txt[len + 1];
    strncpy(txt, (const char*)arr, len);
    txt[len] = '\0';
    xxd_str = txt;
    return xxd_str;
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

int round(float n)
{
    if(n < 0.f){
        return (int)(n - 0.5f);
    }
    return (int)(n + 0.5f);
}

static std::string folder;

const std::string getFolderFromFilePath(const std::string file_path)
{
    std::ostringstream str;
    if(file_path[0] != '/')
        str << "./";
    std::size_t last_sep = file_path.find_last_of("/\\");
    if(last_sep != std::string::npos)
        str << file_path.substr(0, last_sep);
    folder = str.str();
    return folder;
}

static std::string filename;

const std::string getFileNameFromFilePath(const std::string file_path)
{
    std::size_t last_sep = file_path.find_last_of("/\\");
    if(last_sep != std::string::npos)
        filename = file_path.substr(last_sep + 1);
    else
        filename = file_path;
    return filename;
}

static std::string extension;

const std::string getExtensionFromFilePath(const std::string file_path)
{
    std::size_t last_sep = file_path.find_last_of(".");
    if(last_sep != std::string::npos)
        extension = file_path.substr(last_sep + 1);
    else
        extension = "";
    return extension;
}

bool isFile(const std::string file_name)
{
    FILE *streamer=0;
    streamer = fopen(file_name.c_str(), "rb");
    if(streamer){
        fclose(streamer);
        return true;
    }
    return false;
}

size_t readFile(char* source_code, const std::string file_name)
{
    size_t length = 0;
    FILE *file = NULL;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    file = fopen(file_name.c_str(), "rb");
    if(file == NULL) {
        sprintf(msg, "Impossible to open \"%s\".\n", file_name.c_str());
        S->addMessageF(L_ERROR, msg);
        return 0;
    }
    fseek(file, 0, SEEK_END);
    length = ftell(file);
    fseek(file, 0, SEEK_SET);
    if(length == 0) {
        strcpy(msg, "(ReadFile): \"");
        strcat(msg,file_name.c_str());
        strcat(msg,"\" is empty.\n");
        S->addMessageF(L_ERROR, msg);
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
        S->addMessageF(L_ERROR, msg);
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
