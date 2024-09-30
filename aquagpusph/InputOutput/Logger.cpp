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
 * @brief Terminal output, with Log automatic copying.
 * (See Aqua::InputOutput::Logger for details)
 */

#include <sys/time.h>
#include <unistd.h>
#include <sstream>
#include <iostream>

#include "aquagpusph/AuxiliarMethods.hpp"
#include "Logger.hpp"
#include "aquagpusph/ProblemSetup.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace Aqua {
namespace InputOutput {


Logger::Logger()
	: _level(L_DEBUG)
{
#ifdef HAVE_MPI
	auto mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
	if (mpi_rank > 0)
		_level = L_ERROR;
#endif
	open();
	gettimeofday(&_start_time, NULL);
}

Logger::~Logger()
{
	close();
}

void
Logger::writeReport(std::string input)
{
	const std::lock_guard<std::recursive_mutex> lock(_mutex);

	if (!input.size())
		return;

#ifdef HAVE_MPI
	auto mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
	std::cout << "[Proc " << mpi_rank << "] ";
#endif
	std::cout << input;
	if (!hasSuffix(input, "\n")) {
		std::cout << std::endl;
	}
}

void
Logger::addMessage(TLogLevel level, std::string log, std::string func)
{
	const std::lock_guard<std::recursive_mutex> lock(_mutex);

	std::ostringstream fname;
	if (func != "")
		fname << "(" << func << "): ";

	// Send the info to the log file (if possible)
	if (_log_file.is_open()) {
		if (level == L_INFO)
			_log_file << "<b><font color=\"#000000\">[INFO] " << fname.str()
			          << log << "</font></b><br>";
		else if (level == L_WARNING)
			_log_file << "<b><font color=\"#ff9900\">[WARNING] " << fname.str()
			          << log << "</font></b><br>";
		else if (level == L_ERROR)
			_log_file << "<b><font color=\"#dd0000\">[ERROR] " << fname.str()
			          << log << "</font></b><br>";
		else {
			_log_file << "<font color=\"#000000\">" << fname.str() << log
			          << "</font></b>";
			if (hasSuffix(log, "\n"))
				_log_file << "<br>";
		}
		_log_file.flush();
	}

	if (level < _level)
		return;

#ifdef HAVE_MPI
	auto mpi_rank = Aqua::MPI::rank(MPI_COMM_WORLD);
	std::cout << "[Proc " << mpi_rank << "] ";
#endif
	// Just in case the Logger has been destroyed
	if (level == L_INFO)
		std::cout << "INFO ";
	else if (level == L_WARNING)
		std::cout << "WARNING ";
	else if (level == L_ERROR)
		std::cout << "ERROR ";
	std::cout << fname.str() << log << std::flush;
}

void
Logger::printDate(TLogLevel level)
{
	std::ostringstream msg;
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	msg << ctime(&seconds) << std::endl;
	addMessage(level, msg.str());
}

void
Logger::printOpenCLError(cl_int error, TLogLevel level)
{
	switch (error) {
		case CL_DEVICE_NOT_FOUND:
			addMessage(level, "\tCL_DEVICE_NOT_FOUND\n");
			break;
		case CL_DEVICE_NOT_AVAILABLE:
			addMessage(level, "\tCL_DEVICE_NOT_AVAILABLE\n");
			break;
		case CL_COMPILER_NOT_AVAILABLE:
			addMessage(level, "\tCL_COMPILER_NOT_AVAILABLE\n");
			break;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			addMessage(level, "\tCL_MEM_OBJECT_ALLOCATION_FAILURE\n");
			break;
		case CL_OUT_OF_RESOURCES:
			addMessage(level, "\tCL_OUT_OF_RESOURCES\n");
			break;
		case CL_OUT_OF_HOST_MEMORY:
			addMessage(level, "\tCL_OUT_OF_HOST_MEMORY\n");
			break;
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			addMessage(level, "\tCL_PROFILING_INFO_NOT_AVAILABLE\n");
			break;
		case CL_MEM_COPY_OVERLAP:
			addMessage(level, "\tCL_MEM_COPY_OVERLAP\n");
			break;
		case CL_IMAGE_FORMAT_MISMATCH:
			addMessage(level, "\tCL_IMAGE_FORMAT_MISMATCH\n");
			break;
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			addMessage(level, "\tCL_IMAGE_FORMAT_NOT_SUPPORTED\n");
			break;
		case CL_BUILD_PROGRAM_FAILURE:
			addMessage(level, "\tCL_BUILD_PROGRAM_FAILURE\n");
			break;
		case CL_MAP_FAILURE:
			addMessage(level, "\tCL_MAP_FAILURE\n");
			break;
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			addMessage(level, "\tCL_MISALIGNED_SUB_BUFFER_OFFSET\n");
			break;
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			addMessage(level,
			           "\tCL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST\n");
			break;
		case CL_COMPILE_PROGRAM_FAILURE:
			addMessage(level, "\tCL_COMPILE_PROGRAM_FAILURE\n");
			break;
		case CL_LINKER_NOT_AVAILABLE:
			addMessage(level, "\tCL_LINKER_NOT_AVAILABLE\n");
			break;
		case CL_LINK_PROGRAM_FAILURE:
			addMessage(level, "\tCL_LINK_PROGRAM_FAILURE\n");
			break;
		case CL_DEVICE_PARTITION_FAILED:
			addMessage(level, "\tCL_DEVICE_PARTITION_FAILED\n");
			break;
		case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:
			addMessage(level, "\tCL_KERNEL_ARG_INFO_NOT_AVAILABLE\n");
			break;
		case CL_INVALID_VALUE:
			addMessage(level, "\tCL_INVALID_VALUE\n");
			break;
		case CL_INVALID_DEVICE_TYPE:
			addMessage(level, "\tCL_INVALID_DEVICE_TYPE\n");
			break;
		case CL_INVALID_PLATFORM:
			addMessage(level, "\tCL_INVALID_PLATFORM\n");
			break;
		case CL_INVALID_DEVICE:
			addMessage(level, "\tCL_INVALID_DEVICE\n");
			break;
		case CL_INVALID_CONTEXT:
			addMessage(level, "\tCL_INVALID_CONTEXT\n");
			break;
		case CL_INVALID_QUEUE_PROPERTIES:
			addMessage(level, "\tCL_INVALID_QUEUE_PROPERTIES\n");
			break;
		case CL_INVALID_COMMAND_QUEUE:
			addMessage(level, "\tCL_INVALID_COMMAND_QUEUE\n");
			break;
		case CL_INVALID_HOST_PTR:
			addMessage(level, "\tCL_INVALID_HOST_PTR\n");
			break;
		case CL_INVALID_MEM_OBJECT:
			addMessage(level, "\tCL_INVALID_MEM_OBJECT\n");
			break;
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			addMessage(level, "\tCL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
			break;
		case CL_INVALID_IMAGE_SIZE:
			addMessage(level, "\tCL_INVALID_IMAGE_SIZE\n");
			break;
		case CL_INVALID_SAMPLER:
			addMessage(level, "\tCL_INVALID_SAMPLER\n");
			break;
		case CL_INVALID_BINARY:
			addMessage(level, "\tCL_INVALID_BINARY\n");
			break;
		case CL_INVALID_BUILD_OPTIONS:
			addMessage(level, "\tCL_INVALID_BUILD_OPTIONS\n");
			break;
		case CL_INVALID_PROGRAM:
			addMessage(level, "\tCL_INVALID_PROGRAM\n");
			break;
		case CL_INVALID_PROGRAM_EXECUTABLE:
			addMessage(level, "\tCL_INVALID_PROGRAM_EXECUTABLE\n");
			break;
		case CL_INVALID_KERNEL_NAME:
			addMessage(level, "\tCL_INVALID_KERNEL_NAME\n");
			break;
		case CL_INVALID_KERNEL_DEFINITION:
			addMessage(level, "\tCL_INVALID_KERNEL_DEFINITION\n");
			break;
		case CL_INVALID_KERNEL:
			addMessage(level, "\tCL_INVALID_KERNEL\n");
			break;
		case CL_INVALID_ARG_INDEX:
			addMessage(level, "\tCL_INVALID_ARG_INDEX\n");
			break;
		case CL_INVALID_ARG_VALUE:
			addMessage(level, "\tCL_INVALID_ARG_VALUE\n");
			break;
		case CL_INVALID_ARG_SIZE:
			addMessage(level, "\tCL_INVALID_ARG_SIZE\n");
			break;
		case CL_INVALID_KERNEL_ARGS:
			addMessage(level, "\tCL_INVALID_KERNEL_ARGS\n");
			break;
		case CL_INVALID_WORK_DIMENSION:
			addMessage(level, "\tCL_INVALID_WORK_DIMENSION\n");
			break;
		case CL_INVALID_WORK_GROUP_SIZE:
			addMessage(level, "\tCL_INVALID_WORK_GROUP_SIZE\n");
			break;
		case CL_INVALID_WORK_ITEM_SIZE:
			addMessage(level, "\tCL_INVALID_WORK_ITEM_SIZE\n");
			break;
		case CL_INVALID_GLOBAL_OFFSET:
			addMessage(level, "\tCL_INVALID_GLOBAL_OFFSET\n");
			break;
		case CL_INVALID_EVENT_WAIT_LIST:
			addMessage(level, "\tCL_INVALID_EVENT_WAIT_LIST\n");
			break;
		case CL_INVALID_EVENT:
			addMessage(level, "\tCL_INVALID_EVENT\n");
			break;
		case CL_INVALID_OPERATION:
			addMessage(level, "\tCL_INVALID_OPERATION\n");
			break;
		case CL_INVALID_GL_OBJECT:
			addMessage(level, "\tCL_INVALID_GL_OBJECT\n");
			break;
		case CL_INVALID_BUFFER_SIZE:
			addMessage(level, "\tCL_INVALID_BUFFER_SIZE\n");
			break;
		case CL_INVALID_MIP_LEVEL:
			addMessage(level, "\tCL_INVALID_MIP_LEVEL\n");
			break;
		case CL_INVALID_GLOBAL_WORK_SIZE:
			addMessage(level, "\tCL_INVALID_GLOBAL_WORK_SIZE\n");
			break;
		case CL_INVALID_PROPERTY:
			addMessage(level, "\tCL_INVALID_PROPERTY\n");
			break;
		case CL_INVALID_IMAGE_DESCRIPTOR:
			addMessage(level, "\tCL_INVALID_IMAGE_DESCRIPTOR\n");
			break;
		case CL_INVALID_COMPILER_OPTIONS:
			addMessage(level, "\tCL_INVALID_COMPILER_OPTIONS\n");
			break;
		case CL_INVALID_LINKER_OPTIONS:
			addMessage(level, "\tCL_INVALID_LINKER_OPTIONS\n");
			break;
		case CL_INVALID_DEVICE_PARTITION_COUNT:
			addMessage(level, "\tCL_INVALID_DEVICE_PARTITION_COUNT\n");
			break;
		default:
			addMessage(level, "\tUnhandled exception\n");
			break;
	}
}

void
Logger::open()
{
	file("log.proc{mpi_rank}.{index}.html", 0);
	_log_file.open(file().c_str());
	if (!_log_file.is_open()) {
		std::ostringstream msg;
		msg << "Failure creating the log file \"" << file() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Failure creating the Log file");
	}

	_log_file << "<html>" << std::endl;
	_log_file << "<head><title>AQUAgpusph log file.</title></head>"
	          << std::endl;
	_log_file << "<body bgcolor=\"#f0ffff\">" << std::endl;
	_log_file << "<h1 align=\"center\">AQUAgpusph log file.</h1>" << std::endl;
	// Starting data
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	_log_file << "<p align=\"left\">" << ctime(&seconds) << "</p>" << std::endl;
	_log_file << "<hr><br>" << std::endl;
	_log_file.flush();
}

void
Logger::close()
{
	if (!_log_file.is_open()) {
		return;
	}

	_log_file << "<br><hr>" << std::endl;
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	_log_file << "<b><font color=\"#000000\">End of simulation</font></b><br>"
	          << std::endl;
	_log_file << "<p align=\"left\">" << ctime(&seconds) << "</p>" << std::endl;
	_log_file << "</body>" << std::endl;
	_log_file << "</html>" << std::endl;
	_log_file.close();
}

}
} // namespace
