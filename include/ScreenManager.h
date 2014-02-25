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

#ifndef SCREENMANAGER_H_INCLUDED
#define SCREENMANAGER_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

// ----------------------------------------------------------------------------
// Include the ncurses library
// ----------------------------------------------------------------------------
#include <ncurses.h>

// ----------------------------------------------------------------------------
// Include Singleton abstract class
// ----------------------------------------------------------------------------
#include <Singleton.h>

// ----------------------------------------------------------------------------
// Include OpenCL libraries
// ----------------------------------------------------------------------------
#include <CL/cl.h>

#ifndef addMessageF
    #define addMessageF(level, log) addMessage(level, log, __METHOD_CLASS_NAME__)
#endif

namespace Aqua{ namespace InputOutput{

/** @class ScreenManager ScreenManager.h ScreenManager.h
 *  On screen output manager.
 */
struct ScreenManager : public Aqua::Singleton<Aqua::InputOutput::ScreenManager>
{
public:
	/** Constructor.
	 */
	ScreenManager();

	/** Destructor.
	 */
	~ScreenManager();

	/** Some relevant data show will printed at terminal.
	 */
	void update();

	/** Add a message to the log registry, and erase old messages (only in terminal).
	 * @param Level Califier of message:
	 *   - 0 = Empty message.
	 *   - 1 = Info message.
	 *   - 2 = Warning message.
	 *   - 3 = Error message.
	 * @param log Log message.
	 * @param print_func_name Set the function name printing:
	 *   - 0 = Function name must not be printed.
	 *   - 1 = Function name must be printed.
	 *   - 2 = Function name will be printed if level > 0.
	 * @note In order to append the class and the method name before the
	 * messaage use addMessageF instead of this one.
	 */
	void addMessage(int Level, const char *log, const char *func=NULL);

	/** Print time stamp in the screen and log file.
	 */
    void printDate();

	/** Print an OpenCL error.
	 * @param error Error code returned by OpenCL
	 * @param level Califier of message:
	 *   - 0 = Empty message.
	 *   - 1 = Info message.
	 *   - 2 = Warning message.
	 *   - 3 = Error message.
	 */
    void printOpenCLError(int error, int level=0);
private:
	/// Start time
	struct timeval _start_time;
	/// Actual time
	struct timeval _actual_time;
	/// Maximum number of log messages
	long _n_log;
	/** Califier of message: \n
	 * 0 = Empty message. \n
	 * 1 = Info message. \n
	 * 2 = Warning message. \n
	 * 3 = Error message.
	 */
	int *_c_log;
	/// Log messages.
	char **_m_log;

	/// Previous Frame
	int _old_frame;
};

}}  // namespace

#endif // SCREENMANAGER_H_INCLUDED
