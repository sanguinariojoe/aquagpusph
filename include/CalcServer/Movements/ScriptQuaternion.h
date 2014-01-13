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

#ifndef SCRIPTQUATERNION_H_INCLUDED
#define SCRIPTQUATERNION_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Python script embedding
// ----------------------------------------------------------------------------
#include <Python.h>

// ----------------------------------------------------------------------------
// Include the Movement
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/Quaternion.h>

// ----------------------------------------------------------------------------
// Include the Torque computation kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Torque.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class ScriptQuaternion ScriptQuaternion.h CalcServer/Movements/ScriptQuaternion.h
 * @brief Scripted solid quaternion movement. Quaternion, that is specified using an
 * external Python script, defines the solid position at any time.
 */
class ScriptQuaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
	/** Constructor.
	 */
	ScriptQuaternion();

	/** Destructor.
	 */
	~ScriptQuaternion();

	/** Executes the movement.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

protected:
	/** Parse input definition file.
	 * @param root Input node of reader.
	 * @return false if all gone right, true otherwise.
	 */
	bool _parse(xercesc::DOMElement *root);

private:
	/** Init quaternion. init() script method will called
	 * in order to get it.
	 * @return false if all gone right, true otherwise.
	 */
	bool initQuaternion();

	/** Initialize perform() script method to future usage.
	 * @return false if all gone right, true otherwise.
	 */
	bool initPython();

	/** Check an result value in order to know if is a valid
	 * returned value.
	 * @param quat Resultant quaternion.
	 * @return true if valid script output provided, \n false
	 * otherwise.
	 */
	bool isValidOutput(PyObject *quat);

	/** Create arguments tuple.
	 * @param torque Fluid torque measured.
	 * @return Arguments object.
	 */
	PyObject* args(vec torque, vec force);

	/// Python script path
	char* mScript;
	/// Torque calculator
	Torque *_torque;
	/// Python module object
	PyObject *mModule;
	/// Python function
	PyObject *mFunc;
};

}}} // namespace

#endif // SCRIPTQUATERNION_H_INCLUDED
