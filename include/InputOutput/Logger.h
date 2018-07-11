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

#include <sphPrerequisites.h>

#include <string>
#include <fstream>
#include <vector>
#include <CL/cl.h>

#ifdef HAVE_NCURSES
    #include <ncurses.h>
#endif

#include <InputOutput/Report.h>
#include <Singleton.h>

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
    #define LOG(level, log) Aqua::InputOutput::Logger::singleton()->addMessageF(level, log)
#endif
#ifndef LOG0
    /** @def LOG0
     * Overloaded version of
     * Aqua::InputOutput::Logger::singleton()->addMessage(), such that
     * this macro can ve conveniently called to fill the log file.
     */
    #define LOG0(level, log) Aqua::InputOutput::Logger::singleton()->addMessage(level, log)
#endif

namespace Aqua{

enum TLogLevel {L_DEBUG, L_INFO, L_WARNING, L_ERROR};

namespace InputOutput{

/** @class Logger Logger.h Logger.h
 * @brief On screen and log file output manager.
 *
 * AQUAgpusph is generating, during runtime, an HTML log file, placed in the
 * execution folder, and named log.X.html, where X is replaced by the first
 * unsigned integer which generates a non-existing file.
 */
struct Logger : public Aqua::Singleton<Aqua::InputOutput::Logger>
              , public Aqua::InputOutput::Report
{
public:
    /// @brief Constructor.
    Logger();

    /// Destructor.
    ~Logger();

    /** @brief Transform the terminal into an ncurses one.
     *
     * This method should be called ONCE before starting the simulation main
     * loop
     */
    void initNCurses();

    /** @brief Transform the terminal into an ncurses one.
     *
     * This method should be called after finishing the simulation main loop
     */
    void endNCurses();

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
     *    - Tabulators '\\t' are interpreted as 1 blank space
     *    - Line breaks '\\n' are intepreted as line breaks
     *    - It is fitting the message replacing spaces by lines break
     * Otherwise:
     *    - stdout will be used
     *    - A line break '\\n' is appended if it is not detected at the end
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
    void writeReport(std::string msg,
                     std::string color="white",
                     bool bold=false);

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
    void addMessage(TLogLevel level, std::string log, std::string func="");

    /** @brief Print a time stamp in the screen and the log file.
     * @param level Message classification (L_DEBUG, L_INFO, L_WARNING, L_ERROR)
     */
    void printDate(TLogLevel level=L_DEBUG);

    /** @brief Print an OpenCL error.
     * @param error Error code returned by OpenCL.
     * @param level Message classification (L_DEBUG, L_INFO, L_WARNING, L_ERROR)
     */
    void printOpenCLError(cl_int error, TLogLevel level=L_DEBUG);

    /** @brief Do nothing.
     *
     * @param t Simulation time
     */
    void save(float t) {};
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

    /// Create the log file
    void open();

    /// Close the log file
    void close();
private:
    /// Last row where datas was printed (used to locate the registry position)
    int _last_row;

    /// Start time
    struct timeval _start_time;
    /// Actual time
    struct timeval _actual_time;

    /// List of log messages level
    std::vector<int> _log_level;
    /// List of log messages
    std::vector<std::string> _log;
    /// Output log file
    std::ofstream _log_file;
};

}}  // namespace

#endif // LOGGER_H_INCLUDED
