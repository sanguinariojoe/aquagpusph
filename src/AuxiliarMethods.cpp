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
#include <sstream>

#include <AuxiliarMethods.h>
#include <ProblemSetup.h>
#include <InputOutput/Logger.h>

namespace Aqua{

const bool isKeyPressed()
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
        return true;
    }

    return false;
}

const bool hasSuffix(const std::string &str, const std::string &suffix)
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
                           const std::string &search,
                           const std::string &replace)
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

std::string xxd2string(const unsigned char* arr, const unsigned int &len)
{
    char txt[len + 1];
    strncpy(txt, (const char*)arr, len);
    txt[len] = '\0';
    xxd_str = txt;
    return xxd_str;
}

void toLower(std::string &str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

std::string toLowerCopy(std::string str)
{
    toLower(str);
    return str;
}

const unsigned int nextPowerOf2(const unsigned int &n)
{
    if(!n)
        return 1;
    if(isPowerOf2(n))
        return n;

    unsigned int p = 1;
    while(p < n) {
        p <<= 1;
    }
    return p;
}

const unsigned int roundUp(const unsigned int &n, const unsigned int &divisor)
{
    unsigned int result = n;
    const unsigned int rest = n % divisor;
    if(rest) {
        result -= rest;
        result += divisor;
    }
    return result;
}

const int round(const float &n)
{
    if(n < 0.f){
        return (int)(n - 0.5f);
    }
    return (int)(n + 0.5f);
}

static std::string folder;

const std::string getFolderFromFilePath(const std::string &file_path)
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

const std::string getFileNameFromFilePath(const std::string &file_path)
{
    std::size_t last_sep = file_path.find_last_of("/\\");
    if(last_sep != std::string::npos)
        filename = file_path.substr(last_sep + 1);
    else
        filename = file_path;
    return filename;
}

static std::string extension;

const std::string getExtensionFromFilePath(const std::string &file_path)
{
    std::size_t last_sep = file_path.find_last_of(".");
    if(last_sep != std::string::npos)
        extension = file_path.substr(last_sep + 1);
    else
        extension = "";
    return extension;
}

const bool isFile(const std::string &file_name)
{
    std::ifstream f(file_name);
    const bool good = f.good();
    f.close();
    return good;
}

const bool isRelativePath(const std::string &path)
{
    if (ltrimCopy(path).front() == '/')
        return false;
    return true;
}

const size_t getLocalWorkSize(const cl_command_queue &queue)
{
    cl_int err_code;
    cl_device_id d;
    err_code = clGetCommandQueueInfo(queue,
                                     CL_QUEUE_DEVICE,
                                     sizeof(cl_device_id),
                                     &d,
                                     NULL);
    if(err_code != CL_SUCCESS){
        return 0;
    }
    // Start trying maximum local work size per dimension
    cl_uint dims;
    err_code = clGetDeviceInfo(d,
                               CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                               sizeof(cl_uint),
                               &dims,
                               NULL);
    if(err_code != CL_SUCCESS){
        return 0;
    }
    size_t l[dims];
    err_code = clGetDeviceInfo(d,
                               CL_DEVICE_MAX_WORK_ITEM_SIZES,
                               dims * sizeof(size_t),
                               l,
                               NULL);
    if(err_code != CL_SUCCESS){
        return 0;
    }
    // Correct it with maximum local size
    size_t max_l;
    err_code = clGetDeviceInfo(d,
                               CL_DEVICE_MAX_WORK_GROUP_SIZE,
                               sizeof(size_t),
                               &max_l,
                               NULL);
    if(err_code != CL_SUCCESS){
        return 0;
    }
    if(max_l < l[0])
        l[0] = max_l;
    return l[0];
}

const vec mult(const float &n, const vec &v)
{
#ifdef HAVE_3D
    const vec r = {n * v.x, n * v.y, n * v.z, n * v.w};
#else
    const vec r = {n * v.x, n * v.y};
#endif
    return r;
}

const vec add(const vec &a, const vec &b)
{
#ifdef HAVE_3D
    const vec r = {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
#else
    const vec r = {a.x + b.x, a.y + b.y};
#endif
    return r;
}

const vec sub(const vec &a, const vec &b)
{
#ifdef HAVE_3D
    const vec r = {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
#else
    const vec r = {a.x - b.x, a.y - b.y};
#endif
    return r;
}

const float dot(const vec &a, vec &b)
{
    float d = a.x * b.x + a.y * b.y;
#ifdef HAVE_3D
    d += a.z * b.z + a.w * b.w;
#endif
    return d;
}

const float length(const vec &v)
{
    float m = v.x * v.x + v.y * v.y;
#ifdef HAVE_3D
    m += v.z * v.z;
#endif
    return sqrt(m);
}

const vec normalize(const vec &v)
{
    const float m = length(v);
#ifdef HAVE_3D
    const vec r = {v.x / m, v.y / m, v.z / m, v.w};
#else
    const vec r = {v.x / m, v.y / m};
#endif
    return r;
}

#ifdef HAVE_3D
const vec cross(const vec &a, const vec &b)
{
    const vec c = {a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x,
                   0.f};
    return c;
}
#endif

const unsigned int numberOfDigits(const unsigned int number)
{
    return number > 0 ? (int)log10((double)number) + 1 : 1;
}

}   // namespace
