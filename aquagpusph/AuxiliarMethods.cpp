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
#include <iomanip>
#include <numeric>
#include <mutex>

#include "AuxiliarMethods.hpp"
#include "ProblemSetup.hpp"
#include "InputOutput/Logger.hpp"


#include <iostream>

namespace Aqua {

int
isKeyPressed()
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

	if (ch != EOF) {
		ungetc(ch, stdin);
		return 1;
	}

	return 0;
}

bool
hasPrefix(const std::string& str, const std::string& prefix)
{
	return str.size() >= prefix.size() &&
	       str.compare(0, prefix.size(), prefix) == 0;
}

bool
hasSuffix(const std::string& str, const std::string& suffix)
{
	return str.size() >= suffix.size() &&
	       str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void
replaceAll(std::string& str,
           const std::string& search,
           const std::string& replace)
{
	size_t pos = 0;
	while ((pos = str.find(search, pos)) != std::string::npos) {
		str.erase(pos, search.size());
		str.insert(pos, replace);
	}
}

std::string
replaceAllCopy(std::string str, std::string search, std::string replace)
{
	replaceAll(str, search, replace);
	return str;
}

void
ltrim(std::string& s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		        return !std::isspace(ch);
	        }));
}

void
rtrim(std::string& s)
{
	s.erase(std::find_if(
	            s.rbegin(), s.rend(), [](int ch) { return !std::isspace(ch); })
	            .base(),
	        s.end());
}

void
trim(std::string& s)
{
	ltrim(s);
	rtrim(s);
}

std::string
ltrimCopy(std::string s)
{
	ltrim(s);
	return s;
}

std::string
rtrimCopy(std::string s)
{
	rtrim(s);
	return s;
}

std::string
trimCopy(std::string s)
{
	trim(s);
	return s;
}

std::string
xxd2string(unsigned char* arr, unsigned int len)
{
	char txt[len + 1];
	strncpy(txt, (const char*)arr, len);
	txt[len] = '\0';
	std::string xxd_str(txt);
	return xxd_str;
}

void
toLower(std::string& str)
{
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

std::string
toLowerCopy(std::string str)
{
	toLower(str);
	return str;
}

void
setStrConstants(std::string& str)
{
	std::ostringstream number_str;

	int mpi_rank = 0;
#ifdef HAVE_MPI
	mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
#endif
	number_str << mpi_rank;
	replaceAll(str, "{mpi_rank}", number_str.str());
	replaceAll(str, "{version}", PACKAGE_VERSION);
}

std::string
setStrConstantsCopy(std::string str)
{
	setStrConstants(str);
	return str;
}

std::vector<std::string>
split(std::string str, char chr)
{
	std::vector<std::string> substrs;
	std::istringstream istr(str);
	std::string substr;
	while (getline(istr, substr, chr)) {
		substrs.push_back(substr);
	}
	return substrs;
}

std::vector<std::string>
split_formulae(std::string str)
{
	// Replace all the commas outside functions by semicolons, to be taken into
	// account as separators
	std::string edited_str = str;
	int parenthesis_counter = 0;
	for (auto it = edited_str.begin(); it != edited_str.end(); ++it) {
		// We does not care about unbalanced parenthesis, muparser will do it
		if (*it == '(')
			parenthesis_counter++;
		else if (*it == ')')
			parenthesis_counter--;
		else if ((*it == ',') && (parenthesis_counter == 0)) {
			*it = ';';
		}
	}

	return split(edited_str, ';');
}

std::string
newFilePath(const std::string& basename, unsigned int& i, unsigned int digits)
{
	FILE* f;
	std::string filepath;
	std::ostringstream number_str;

	// Start replacing all the old-school formatting string instances by the new
	// one, based on a more intelligible variable name
	filepath = replaceAllCopy(basename, "%d", "{index}");

	// Set the constants
	setStrConstants(filepath);

	if (filepath.find("{index}") == std::string::npos) {
		// We cannot insert the file index anywhere, so just test if the file
		// does not exist
		f = fopen(filepath.c_str(), "r");
		if (f) {
			fclose(f);
			throw std::invalid_argument("Invalid file name pattern");
		}

		return filepath;
	}

	while (true) {
		number_str.str("");
		number_str << std::setfill('0') << std::setw(digits) << i;
		std::string newfilepath =
		    replaceAllCopy(filepath, "{index}", number_str.str());

		f = fopen(newfilepath.c_str(), "r");
		if (!f) {
			filepath = newfilepath;
			break;
		}
		fclose(f);
		i++;
	}

	return filepath;
}

unsigned int
nextPowerOf2(unsigned int n)
{
	if (n & !(n & (n - 1)))
		return n;

	unsigned int p = 1;
	while (p < n) {
		p <<= 1;
	}
	return p;
}

unsigned int
isPowerOf2(unsigned int n)
{
	return ((n & (n - 1)) == 0);
}

static std::string folder;

const std::string
getFolderFromFilePath(const std::string file_path)
{
	std::ostringstream str;
	if (file_path[0] != '/')
		str << "./";
	std::size_t last_sep = file_path.find_last_of("/\\");
	if (last_sep != std::string::npos)
		str << file_path.substr(0, last_sep);
	folder = str.str();
	return folder;
}

static std::string filename;

const std::string
getFileNameFromFilePath(const std::string file_path)
{
	std::size_t last_sep = file_path.find_last_of("/\\");
	if (last_sep != std::string::npos)
		filename = file_path.substr(last_sep + 1);
	else
		filename = file_path;
	return filename;
}

static std::string extension;

const std::string
getExtensionFromFilePath(const std::string file_path)
{
	std::size_t last_sep = file_path.find_last_of(".");
	if (last_sep != std::string::npos)
		extension = file_path.substr(last_sep + 1);
	else
		extension = "";
	return extension;
}

bool
isFile(const std::string file_name)
{
	std::ifstream f(file_name);
	return f.good();
}

bool
isRelativePath(const std::string path)
{
	if (trimCopy(path).front() == '/')
		return false;
	return true;
}

size_t
getLocalWorkSize(cl_uint n, cl_command_queue queue)
{
	cl_int flag;
	cl_device_id d;
	flag = clGetCommandQueueInfo(
	    queue, CL_QUEUE_DEVICE, sizeof(cl_device_id), &d, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	// Start trying maximum local work size per dimension
	cl_uint dims;
	flag = clGetDeviceInfo(
	    d, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &dims, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	size_t l[dims];
	flag = clGetDeviceInfo(
	    d, CL_DEVICE_MAX_WORK_ITEM_SIZES, dims * sizeof(size_t), l, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	// Correct it with maximum local size
	size_t max_l;
	flag = clGetDeviceInfo(
	    d, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_l, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	if (max_l < l[0])
		l[0] = max_l;
	return l[0];
}

size_t
getGlobalWorkSize(cl_uint n, size_t local_work_size)
{
	return roundUp(n, local_work_size);
}

vec
Vzero()
{
	vec r;
	r.x = 0.f;
	r.y = 0.f;
#ifdef HAVE_3D
	r.z = 0.f;
	r.w = 0.f;
#endif
	return r;
}

vec
Vx()
{
	vec r;
	r.x = 1.f;
	r.y = 0.f;
#ifdef HAVE_3D
	r.z = 0.f;
	r.w = 0.f;
#endif
	return r;
}

vec
Vy()
{
	vec r;
	r.x = 0.f;
	r.y = 1.f;
#ifdef HAVE_3D
	r.z = 0.f;
	r.w = 0.f;
#endif
	return r;
}

#ifdef HAVE_3D
vec
Vz()
{
	vec r;
	r.x = 0.f;
	r.y = 0.f;
	r.z = 1.f;
	r.w = 0.f;
	return r;
}
#endif

vec
mult(float n, vec v)
{
	vec r;
	r.x = n * v.x;
	r.y = n * v.y;
#ifdef HAVE_3D
	r.z = n + v.z;
	r.w = n + v.w;
#endif
	return r;
}

vec
add(vec a, vec b)
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

vec
sub(vec a, vec b)
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

float
dot(vec a, vec b)
{
	float d = a.x * b.x + a.y * b.y;
#ifdef HAVE_3D
	d += a.z * b.z;
	d += a.w * b.w;
#endif
	return d;
}

float
length(vec v)
{
	float m = v.x * v.x + v.y * v.y;
#ifdef HAVE_3D
	m += v.z * v.z;
#endif
	return sqrt(m);
}

vec
normalize(vec v)
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
vec
cross(vec a, vec b)
{
	vec c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	c.w = 0.f;
	return c;
}
#endif

unsigned int
numberOfDigits(unsigned int number)
{
	return number > 0 ? (int)log10((double)number) + 1 : 1;
}

#ifdef HAVE_MPI

namespace MPI {

void error_handler(MPI_Comm *comm, int *err, ...) {
	LOG(L_ERROR, "MPI reported an error: " + error_str(*err) + "\n");
}

std::string error_str(int errorcode)
{
	char str[MPI_MAX_ERROR_STRING + 1];
	int l;
	MPI_Error_string(errorcode, str, &l);
	if (l > MPI_MAX_ERROR_STRING)
		l = MPI_MAX_ERROR_STRING;
	str[l] = '\0';
	return std::string(str);
}

void init(int *argc, char ***argv)
{
	int err = MPI_Init(argc, argv);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Init() failed: " + error_str(err) + "\n");
		throw std::runtime_error("MPI error");
	}
	MPI_Errhandler mpi_err_handler;
	err = MPI_Comm_create_errhandler(&Aqua::MPI::error_handler,
	                                 &mpi_err_handler);
	if (err != MPI_SUCCESS) {
		LOG(L_WARNING, "Failure creating a MPI error handler: " +
			error_str(err) + "\n");
	}
    err = MPI_Comm_set_errhandler(MPI_COMM_WORLD, mpi_err_handler);
	if (err != MPI_SUCCESS) {
		LOG(L_WARNING, "Failure installing the MPI error handler: " +
			error_str(err) + "\n");
	}
}

void finalize()
{
	const int err = MPI_Finalize();
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Finalize() failed: " + error_str(err) + "\n");
		throw std::runtime_error("MPI error");
	}
}

int rank(MPI_Comm comm)
{
	int rank;
	const int err = MPI_Comm_rank(comm, &rank);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Comm_rank() failed: " + error_str(err) + "\n");
		throw std::runtime_error("MPI error");
	}
	return rank;
}

int size(MPI_Comm comm)
{
	int size;
	const int err = MPI_Comm_size(comm, &size);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Comm_size() failed: " + error_str(err) + "\n");
		throw std::runtime_error("MPI error");
	}
	return size;
}

void barrier(MPI_Comm comm)
{
	const int err = MPI_Barrier(comm);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Barrier() failed: " + error_str(err) + "\n");
		MPI_Abort(comm, err);
		throw std::runtime_error("MPI error");
	}
}

void send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
          MPI_Comm comm)
{
	const int err = MPI_Send(buf, count, datatype, dest, tag, comm);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Isend() failed: " + error_str(err) + "\n");
		MPI_Abort(comm, err);
		throw std::runtime_error("MPI error");
	}
}

MPI_Request isend(const void *buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm)
{
	MPI_Request req;
	const int err = MPI_Isend(buf, count, datatype, dest, tag, comm, &req);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Isend() failed: " + error_str(err) + "\n");
		MPI_Abort(comm, err);
		throw std::runtime_error("MPI error");
	}
	return req;
}

MPI_Status recv(void *buf, int count, MPI_Datatype datatype, int source,
                int tag, MPI_Comm comm)
{
	MPI_Status status;
	const int err = MPI_Recv(buf, count, datatype, source, tag, comm, &status);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Irecv() failed: " + error_str(err) + "\n");
		MPI_Abort(comm, err);
		throw std::runtime_error("MPI error");
	}
	return status;
}

MPI_Request irecv(void *buf, int count, MPI_Datatype datatype,
                  int source, int tag, MPI_Comm comm)
{
	MPI_Request req;
	const int err = MPI_Irecv(buf, count, datatype, source, tag, comm, &req);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Irecv() failed: " + error_str(err) + "\n");
		MPI_Abort(comm, err);
		throw std::runtime_error("MPI error");
	}
	return req;
}

MPI_Status wait(MPI_Request *request)
{
	MPI_Status status;
	const int err = MPI_Wait(request, &status);
	if (err != MPI_SUCCESS) {
		LOG(L_ERROR, "MPI_Wait() failed: " + error_str(err) + "\n");
		status.MPI_ERROR = err;
	}
	return status;
}

} // ::MPI

#endif


} // ::Aqua
