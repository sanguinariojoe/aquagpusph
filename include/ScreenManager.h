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
#include <deque>

#ifdef HAVE_NCURSES
    #include <ncurses.h>
#endif

#include <Singleton.h>

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

    /** @brief Call to setup a new terminal frame.
     *
     * This method should be called at the start of every time step.
     */
    void initFrame();

    /** @brief Call to refresh the terminal frame.
     *
     * This method should be called at the end of every time step.
     */
    void endFrame();

    /** @brief Write a new message in the terminal output.
     *
     * This method is not redirecting the data to the log file.
     * In case that ncurses is active:
     *    - Tabulators '\t' are interpreted as 1 blank space
     *    - Line breaks '\n' are intepreted as line breaks
     *    - It is fitting the message replacing spaces by lines break
     * Otherwise:
     *    - stdout will be used
     *    - A line break '\n' is appended if it is not detected at the end
     *
     * @param msg Message to print in the screen.
     * @param color Color name. Valid colors are:
     *    -# white
     *    -# green
     *    -# blue
     *    -# yellow
     *    -# red
     *    -# magenta
     *    -# cyan
     * @param bold true if bold font should be used, false otherwise
     */
    void writeReport(const char *msg,
                     const char *color="white",
                     bool bold=false);

    /** @brief Add a new log record message.
     *
     * The old messages may be removed from the terminal if no more space left.
     *
     * @param level Califier of message:
     *   - 0 = Empty message.
     *   - 1 = Info message.
     *   - 2 = Warning message.
     *   - 3 = Error message.
     * @param log Log message.
     * @param func Function name to print, NULL if it should not be printed.
     * @note In order to append the class and the method name before the
     * message use #addMessageF instead of this one.
     */
    void addMessage(int level, const char *log, const char *func=NULL);

    /** @brief Print a time stamp in the screen and the log file.
     * @param level Qualifier of message:
     *   - 0 = Empty message.
     *   - 1 = Info message.
     *   - 2 = Warning message.
     *   - 3 = Error message.
     */
    void printDate(int level=0);

    /** @brief Print an OpenCL error.
     * @param error Error code returned by OpenCL.
     * @param level Qualifier of message:
     *   - 0 = Empty message.
     *   - 1 = Info message.
     *   - 2 = Warning message.
     *   - 3 = Error message.
     */
    void printOpenCLError(int error, int level=0);
protected:
    /** @brief Print the log record
     *
     * This function will compute automatically where the log record should be
     * placed
     */
    void printLog();

    /** Call to refresh the output screen, useless if ncurses is not active.
     */
    void refreshAll();
private:

    /// Last row where datas was printed (used to locate the registry position)
    int _last_row;

    /// Start time
    struct timeval _start_time;
    /// Actual time
    struct timeval _actual_time;

    /// List of log messages level
    std::deque<int> _log_level;
    /// List of log messages
    std::deque<char*> _log;
};

}}  // namespace

#endif // SCREENMANAGER_H_INCLUDED
