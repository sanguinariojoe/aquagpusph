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
#include <string>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <mutex>

#ifdef _WIN32
#include <windows.h>
#else
#include <limits.h>
#include <unistd.h>
#endif

#include "AuxiliarMethods.hpp"
#include "ProblemSetup.hpp"
#include "InputOutput/Logger.hpp"


#include <iostream>

namespace Aqua {

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
	char* txt = (char*)malloc((len + 1) * sizeof(char));
	if (!txt) {
		LOG(L_ERROR, std::string("Failure allocating ") +
		             std::to_string((len + 1) * sizeof(char)) + "bytes\n");
		throw std::bad_alloc();
	}
	strncpy(txt, (const char*)arr, len);
	txt[len] = '\0';
	std::string xxd_str(txt);
	free(txt);
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

const std::string
getExePath()
{
	#ifdef _WIN32
		wchar_t path[MAX_PATH] = { 0 };
		GetModuleFileNameW(NULL, path, MAX_PATH);
		std::filesystem::path fp(path);
		return fp.string();
	#else
		char result[PATH_MAX];
		ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
		return std::string(result, (count > 0) ? count : 0);
	#endif
}

const std::string
getRootPath()
{
	std::filesystem::path ep(getFolderFromFilePath(getExePath()));
	if (std::filesystem::exists(ep / "../resources")) {
		return std::filesystem::canonical(ep / "../").string();
	}
	return std::filesystem::canonical(ep / "../share/aquagpusph/").string();
}

const std::string
getFolderFromFilePath(const std::string file_path)
{
	std::filesystem::path fp(file_path);
	return std::filesystem::canonical(fp.parent_path()).string();
}

const std::string
getFileNameFromFilePath(const std::string file_path)
{
	std::filesystem::path fp(file_path);
	return fp.filename().string();
}

const std::string
getExtensionFromFilePath(const std::string file_path)
{
	std::filesystem::path fp(file_path);
	return fp.extension().string();
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
getLocalWorkSize(cl_command_queue queue)
{
	cl_int flag;
	cl_device_id d;
	flag = clGetCommandQueueInfo(
	    queue, CL_QUEUE_DEVICE, sizeof(cl_device_id), &d, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	size_t max_l;
	flag = clGetDeviceInfo(
	    d, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_l, NULL);
	if (flag != CL_SUCCESS) {
		return 0;
	}
	return max_l;
}

size_t
getGlobalWorkSize(size_t n, size_t local_work_size)
{
	return roundUp<size_t>(n, local_work_size);
}

unsigned int
numberOfDigits(unsigned int number)
{
	return number > 0 ? (int)log10((double)number) + 1 : 1;
}

#ifdef HAVE_MPI

namespace MPI {

void error_handler(MPI_Comm UNUSED_PARAM *comm, int *err, ...) {
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
