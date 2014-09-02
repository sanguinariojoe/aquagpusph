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
 * (See Aqua::InputOutput::ScreenManager for details)
 */

#ifndef SCREENMANAGER_H_INCLUDED
#define SCREENMANAGER_H_INCLUDED

#include <sphPrerequisites.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <ncurses.h>

#include <Singleton.h>

#include <CL/cl.h>

#ifndef addMessageF
    /** @def addMessageF
     * Overloaded version of Aqua::InputOutput::ScreenManager::addMessage()
     * method, which is automatically setting the class and function names.
     */
    #define addMessageF(level, log) addMessage(level, log, __METHOD_CLASS_NAME__)
#endif

namespace Aqua{ namespace InputOutput{

/** @class ScreenManager ScreenManager.h ScreenManager.h
 * @brief On screen output manager.
 *
 * This class is able to conveniently redirect the data to the log file.
 */
struct ScreenManager : public Aqua::Singleton<Aqua::InputOutput::ScreenManager>
{
public:
    /// Constructor.
    ScreenManager();

    /// Destructor.
    ~ScreenManager();

    /** @brief Call to update the terminal output.
     *
     * It may compute some additional data to print by terminal, depending on
     * the selected verbosity.
     *
     * @see Aqua::InputOutput::ProblemSetup::sphSettings::verbose_level
     */
    void update();

    /** @brief Add a new log record message.
     *
     * The old messages may be removed from the terminal if no more space left.
     *
     * @param Level Califier of message:
     *   - 0 = Empty message.
     *   - 1 = Info message.
     *   - 2 = Warning message.
     *   - 3 = Error message.
     * @param log Log message.
     * @param func Function name to print, NULL if it should not be printed.
     * @note In order to append the class and the method name before the
     * message use #addMessageF instead of this one.
     */
    void addMessage(int Level, const char *log, const char *func=NULL);

    /// Print a time stamp in the screen and the log file.
    void printDate();

    /** @brief Print an OpenCL error.
     * @param error Error code returned by OpenCL.
     * @param level Qualifier of message:
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
    /** @brief Qualifier of the messages
     *
     *    - 0 = Empty message.
     *    - 1 = Info message.
     *    - 2 = Warning message.
     *    - 3 = Error message.
     */
    int *_c_log;
    /// Log messages.
    char **_m_log;

    /// Previous Frame index.
    int _old_frame;
};

}}  // namespace

#endif // SCREENMANAGER_H_INCLUDED
