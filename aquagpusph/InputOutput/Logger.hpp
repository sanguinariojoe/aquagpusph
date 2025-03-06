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

#ifndef LOGGER_H_INCLUDED
#define LOGGER_H_INCLUDED

#include "aquagpusph/sphPrerequisites.hpp"

#include <string>
#include <fstream>
#include <vector>
#include <mutex>
#if __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "Report.hpp"

#ifndef addMessageF
/** @def addMessageF
 * Overloaded version of Aqua::InputOutput::Logger::addMessage()
 * method, which is automatically setting the class and function names.
 */
#define addMessageF(level, log) addMessage(level, log, __METHOD_CLASS_NAME__)
#endif

#ifndef LOG
/** @def LOG
 * Overloaded version of
 * Aqua::InputOutput::Logger::singleton()->addMessageF(), such that
 * this macro can ve conveniently called to fill the log file.
 */
#define LOG(level, log)                                                        \
	Aqua::InputOutput::Logger::singleton()->addMessageF(level, log)
#endif
#ifndef LOG0
/** @def LOG0
 * Overloaded version of
 * Aqua::InputOutput::Logger::singleton()->addMessage(), such that
 * this macro can ve conveniently called to fill the log file.
 */
#define LOG0(level, log)                                                       \
	Aqua::InputOutput::Logger::singleton()->addMessage(level, log)
#endif

namespace Aqua {

enum TLogLevel
{
	L_DEBUG,
	L_INFO,
	L_WARNING,
	L_ERROR,
	L_NONE
};

namespace InputOutput {

/** @class Logger Logger.h Logger.h
 * @brief On screen and log file output manager.
 *
 * AQUAgpusph is generating, during runtime, an HTML log file, placed in the
 * execution folder, and named log.X.html, where X is replaced by the first
 * unsigned integer which generates a non-existing file.
 */
class DECLDIR Logger
  : public Aqua::InputOutput::Report
{
public:
	/// Destructor.
	~Logger();

    // Deleting the copy constructor to prevent copies
    Logger(const Logger& obj) = delete;

	/** @brief Get --creating it the first time-- the logger instance
	 *
	 * Just one logger can simultaneously exist
	 * @return The logger
	 */
	static Logger* singleton();

	/** @brief Set the log level
	 * @param level The minimum message log level to get printed
	 */
	inline void setLevel(int level)
	{
		_level = (std::min)((std::max)(level, (int)L_DEBUG), (int)L_NONE);
	}

	/** @brief Write a new message in the terminal output.
	 *
	 * This method is not redirecting the data to the log file.
	 * A line break '\\n' is appended if it is not detected at the end
	 * @param msg Message to print in the screen.
	 */
	void writeReport(std::string msg);

	/** @brief Add a new log record message.
	 *
	 * The old messages may be removed from the terminal if no more space left.
	 *
	 * @param level Message classification (L_DEBUG, L_INFO, L_WARNING, L_ERROR)
	 * @param log Log message.
	 * @param func Function name to print, NULL if it should not be printed.
	 * @note In order to append the class and the method name before the
	 * message use #addMessageF instead of this one.
	 */
	void addMessage(TLogLevel level, std::string log, std::string func = "");

	/** @brief Print a time stamp in the screen and the log file.
	 * @param level Message classification (L_DEBUG, L_INFO, L_WARNING, L_ERROR)
	 */
	void printDate(TLogLevel level = L_DEBUG);

	/** @brief Print an OpenCL error.
	 * @param error Error code returned by OpenCL.
	 * @param level Message classification (L_DEBUG, L_INFO, L_WARNING, L_ERROR)
	 */
	void printOpenCLError(cl_int error, TLogLevel level = L_DEBUG);

	/** @brief Do nothing.
	 *
	 * @param t Simulation time
	 */
	void save(float UNUSED_PARAM t) {};

  protected:
	/// Create the log file
	void open();

	/// Close the log file
	void close();

  private:
	/** @brief Constructor
	 *
	 * This calss is a singleton, and thus we keep the constructor away
	 */
	Logger();

	/// Minimum log level to get the message printed
	int _level;

	/// Output log file
	std::ofstream _log_file;

	/// Mutex to avoid several threads printing at the same time
	std::recursive_mutex _mutex;
};

}
} // namespace

#endif // LOGGER_H_INCLUDED
