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
 * @brief Python casting helpers
 * (see Aqua::InpuOutput::Variable and Aqua::InpuOutput::Variables)
 */

#ifndef PYCAST_H_INCLUDED
#define PYCAST_H_INCLUDED

#include <string>
#include <limits>
#include <sphPrerequisites.h>
#include <InputOutput/Logger.h>

#include <Python.h>

namespace Aqua{ namespace InputOutput{

/** @class PyCast Variable.h Variable.h
 * @brief A class to carry out conversions between Python objects and raw
 * variables, and viceversa
 */
template <class T>
class PyCast
{
public:
    /** @brief Constructor
     */
    PyCast() {};

    /** @brief Destructor
     */
    virtual ~PyCast() {};

protected:
    /** @brief Convert a raw variable into a Python object
     * @param value Value to convert
     * @return Python object, NULL upon failure conversion
     */
    PyObject* valToPython(const T& value) const {
        PyErr_Format(PyExc_TypeError,
                     "PyCast::valToPython cannot handle \"%s\" variables",
                     typeid(value).name());
        return NULL;
    }

    /** @brief Convert a Python object to a raw variable pointer
     * @param obj Python object
     * @return Pointer to value memory address, NULL upon failure conversion
     * @warning This is thread unsafe
     * @note The returned pointer will be overwritten next time this function is
     * called
     */
    void* PythonToPtr(PyObject* obj) {
        PyErr_Format(PyExc_TypeError,
                     "PyCast::PythonToPtr cannot handle \"%s\" variables",
                     typeid(_py_ptr_val).name());
        return NULL;
    }

private:
    /// Variable value, use to return pointers to already allocated memory
    T _py_ptr_val;
};

}}  // namespace

#endif // PYCAST_H_INCLUDED
