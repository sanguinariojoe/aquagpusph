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

#include <PyCast.h>


namespace Aqua{ namespace InputOutput{


template<>
PyObject* PyCast<int>::valToPython(const int& value) const
{
    return PyLong_FromLong(value);
}

template<>
void* PyCast<int>::PythonToPtr(PyObject* obj){
    const long val = PyLong_AsLong(obj);
    if(val == -1 && PyErr_Occurred()) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError,
                     "PyCast<int>::PythonToPtr got \"%s\" variable",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    if((val < std::numeric_limits<int>::min()) ||
       (val > std::numeric_limits<int>::max())) {
        PyErr_Format(PyExc_OverflowError,
                     "\"%s\" variable overflows \"int\" type",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    _py_ptr_val = (int)val;
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<unsigned int>::valToPython(const unsigned int& value) const
{
    return PyLong_FromUnsignedLong(value);
}

template<>
void* PyCast<unsigned int>::PythonToPtr(PyObject* obj){
    const long val = PyLong_AsLong(obj);
    if(val == -1 && PyErr_Occurred()) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError,
                     "PyCast<unsigned int>::PythonToPtr got \"%s\" variable",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    if((val < 0) || (val > std::numeric_limits< int>::max())) {
        PyErr_Format(PyExc_OverflowError,
                     "\"%s\" variable overflows \"unsigned int\" type",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    _py_ptr_val = (unsigned int)val;
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<float>::valToPython(const float& value) const
{
    return PyFloat_FromDouble(value);
}

template<>
void* PyCast<float>::PythonToPtr(PyObject* obj){
    const double val = PyFloat_AsDouble(obj);
    if(val == -1.0 && PyErr_Occurred()) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError,
                     "PyCast<float>::PythonToPtr got \"%s\" variable",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    if((val < std::numeric_limits<float>::min()) ||
       (val > std::numeric_limits<float>::max())) {
        PyErr_Format(PyExc_OverflowError,
                     "\"%s\" variable overflows \"float\" type",
                     Py_TYPE(obj)->tp_name);
        return NULL;
    }
    _py_ptr_val = (float)val;
    return (void*)&_py_ptr_val;
}

}}  // namespace
