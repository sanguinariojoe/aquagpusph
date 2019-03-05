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
 * @brief Python script execution tool.
 * (See Aqua::CalcServer::Python for details)
 */

#ifndef PYTHON_H_INCLUDED
#define PYTHON_H_INCLUDED

#include <Python.h>
#include <CalcServer.h>
#include <CalcServer/Tool.h>

namespace Aqua{ namespace CalcServer{

/** @class Python Python.h CalcServer/Python.h
 * @brief Execute a Python script.
 *
 * The Python script should has a function main which will be called each time
 * step, when execute is called.
 *
 * AQUAgpusph is providing a module called aquagpusph which allows the Python
 * script to get and set variable values.
 */
class Python : public Aqua::CalcServer::Tool
{
public:
    /** @brief Constructor.
     * @param tool_name Tool name.
     * @param script Python script path.
     * @param once Run this tool just once. Useful to make initializations.
     */
    Python(const std::string tool_name,
           const std::string script,
           bool once=false);

    /// Destructor.
    ~Python();

    /** @brief Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accesing the dependencies
     */
    cl_event _execute(const std::vector<cl_event> events);

    /** @brief Initialize the Python interpreter.
     *
     * This method is safely creating the Python environment just one time,
     * adding the execution folder to the system path (in order to import
     * modules).
     */
    void initPython();

    /** @brief Load the script and extract the callable function.
     */
    void load();

private:
    /// Script path
    std::string _script;

    /// Python module object
    PyObject *_module;
    /// Python function to be called
    PyObject *_func;
};

}}  // namespace

#endif // PYTHON_H_INCLUDED
