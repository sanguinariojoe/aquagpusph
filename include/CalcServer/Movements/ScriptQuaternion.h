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
 * @brief Scripted quaternion based motion.
 * (See Aqua::CalcServer::Movement::ScriptQuaternion for details)
 */

#ifndef SCRIPTQUATERNION_H_INCLUDED
#define SCRIPTQUATERNION_H_INCLUDED

#include <Python.h>
#include <CalcServer/Movements/Quaternion.h>
#include <CalcServer/Torque.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class ScriptQuaternion ScriptQuaternion.h CalcServer/Movements/ScriptQuaternion.h
 * @brief Scripted quaternion based motion.
 *
 * In this model the motion is not prescribed in a data file, but it is obtained
 * in runtime from a Python script.
 *
 * The Python script must have 2 methods:
 *   - `init()`: It should return the quaternion at \f$ t = 0 \f$ seconds.
 *   - `perform(COR, X, Y, Z, Torque, Force, t, dt)`: It should return the
 *     quaternion at each time instant.
 *
 * @see Aqua::CalcServer::Movement::Quaternion
 * @see Aqua::CalcServer::Torque
 * perezrojas_etal_stab_2012 example
 */
class ScriptQuaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
    /// Constructor.
    ScriptQuaternion();

    /// Destructor.
    ~ScriptQuaternion();

    /** @brief Execute the motion.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

protected:
    /** @brief Get the Python script file from the XML definition file.
     * @param root Input node of the parser.
     * @return false if all gone right, true otherwise.
     */
    bool _parse(xercesc::DOMElement *root);

private:
    /** @brief Initialize the quaternion.
     *
     * `init()` script method is called in order to get it.
     *
     * @return false if all gone right, true otherwise.
     */
    bool initQuaternion();

    /** @brief Setup Python, loading the script.
     * @return false if all gone right, true otherwise.
     */
    bool initPython();

    /** @brief Check a result value validity.
     * @param quat Resultant quaternion.
     * @return true if valid script output provided, false otherwise.
     */
    bool isValidOutput(PyObject *quat);

    /** Create the arguments tuple.
     * @param torque Fluid moment measured.
     * @param force Fluid force measured.
     * @return The arguments object.
     * @see Aqua::CalcServer::Torque
     */
    PyObject* args(vec torque, vec force);

    /// Python script path
    char* _script;
    /// Forces and moments calculator
    Torque *_torque;
    /// Python module object
    PyObject *_module;
    /// Python function
    PyObject *_func;
};

}}} // namespace

#endif // SCRIPTQUATERNION_H_INCLUDED
