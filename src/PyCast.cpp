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
/** @def PY_ARRAY_UNIQUE_SYMBOL
 * @brief Define the extension module which this Python stuff should be linked
 * to.
 *
 * In AQUAgpusph all the Python stuff is linked in the same group AQUA_ARRAY_API
 * @see http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define PY_ARRAY_UNIQUE_SYMBOL AQUA_ARRAY_API
/** @def NO_IMPORT_ARRAY
 * @brief Set this file as a helper of the group AQUA_ARRAY_API.
 * @see http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define NO_IMPORT_ARRAY
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

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

static char pyarray_type_name[32];

/** @brief Get the Python array inner type name
 * @param pyarray_type_id Type identifier
 * @return Type name
 * @warning This not a thread safe function
 */
const char* pyArrayTypeName(const int& pyarray_type_id)
{
    switch(pyarray_type_id)
    {
        case PyArray_FLOAT:
            strcpy(pyarray_type_name, "PyArray_FLOAT");
            break;
        case PyArray_INT:
            strcpy(pyarray_type_name, "PyArray_INT");
            break;
        case PyArray_UINT:
            strcpy(pyarray_type_name, "PyArray_UINT");
            break;
        default:
            strcpy(pyarray_type_name, "PyArray_UNKNOWN");
    }
    return pyarray_type_name;
}

/** @brief Check that a Python object is an array with the desired features
 *
 * This function is checking the the Python object is an array with specific
 * dimensions and type
 *
 * @param obj Python object
 * @param nd Number of dimensions
 * @param dimensions Dimensions array, with \p nd components
 * @param innertype Inner array type (PyArray_FLOAT / PyArray_INT / PyArray_UINT)
 * @return true if the array has the expected specifications, false otherwise
 * @note This function is registering Python exceptions
 */
const bool checkPyhonObjectArray(PyObject* obj,
                                 const int& nd,
                                 const npy_intp* dimensions,
                                 const int innertype){
    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        PyErr_Format(PyExc_TypeError,
                     "PyArrayObject expected, but \"%s\" has been received",
                     Py_TYPE(obj)->tp_name);
        return false;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;

    const int typ = PyArray_TYPE(array_obj);
    if (PyArray_TYPE(array_obj) != innertype) {
        char tname[32], visitor_tname[32];
        strncpy(tname, pyArrayTypeName(innertype), 32);
        strncpy(visitor_tname, pyArrayTypeName(PyArray_TYPE(array_obj)), 32);
        PyErr_Format(PyExc_TypeError,
                     "Array of type \"%s\" expected, but got \"%s\"",
                     tname, visitor_tname);
        return false;
    }

    if(array_obj->nd != nd){
        PyErr_Format(PyExc_ValueError,
                     "PyArrayObject of %d dims expected, but got %d",
                     nd, array_obj->nd);
        return false;
    }

    for (int i = 0; i < nd; i++) {
        if(array_obj->dimensions[i] != dimensions[i]) {
            PyErr_Format(PyExc_ValueError,
                        "length = %d expected along axis %d, but got %d",
                        dimensions[i], i, array_obj->dimensions[i]);
            return false;
        }
    }

    return true;
}

// NOTE: vec3 is just a typedef of vec4, inherited from the following
// CL/cl_platform.h line:
//
// typedef  cl_float4  cl_float3;
//
// The same apply to ivec and uivec types

template<>
PyObject* PyCast<vec2>::valToPython(const vec2& value) const
{
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, (void*)value.s);
}

template<>
void* PyCast<vec2>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {2};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_FLOAT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(vec2));
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<vec4>::valToPython(const vec4& value) const
{
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, (void*)value.s);
}

template<>
void* PyCast<vec4>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {4};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_FLOAT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(vec4));
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<ivec2>::valToPython(const ivec2& value) const
{
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_INT, (void*)value.s);
}

template<>
void* PyCast<ivec2>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {2};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_INT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(ivec2));
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<ivec4>::valToPython(const ivec4& value) const
{
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_INT, (void*)value.s);
}

template<>
void* PyCast<ivec4>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {4};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_INT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(ivec4));
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<uivec2>::valToPython(const uivec2& value) const
{
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_UINT, (void*)value.s);
}

template<>
void* PyCast<uivec2>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {2};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_UINT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(uivec2));
    return (void*)&_py_ptr_val;
}

template<>
PyObject* PyCast<uivec4>::valToPython(const uivec4& value) const
{
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_UINT, (void*)value.s);
}

template<>
void* PyCast<uivec4>::PythonToPtr(PyObject* obj){
    const npy_intp dims[] = {4};
    if(!checkPyhonObjectArray(obj, 1, dims, PyArray_UINT)) {
        return NULL;
    }
    PyArrayObject* array_obj = (PyArrayObject*) obj;
    memcpy(_py_ptr_val.s, array_obj->data, sizeof(uivec4));
    return (void*)&_py_ptr_val;
}

}}  // namespace
