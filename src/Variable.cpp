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

#include <Variable.h>
#include <AuxiliarMethods.h>
#include <ScreenManager.h>
#include <CalcServer.h>

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

Variable::Variable(const char *varname, const char *vartype)
    : _name(NULL)
    , _typename(NULL)
{
    unsigned int len;

    len = strlen(varname) + 1;
    _name = new char[len];
    strcpy(_name, varname);
    len = strlen(vartype) + 1;
    _typename = new char[len];
    strcpy(_typename, vartype);
}

Variable::~Variable()
{
    delete[] _name; _name = NULL;
    delete[] _typename; _typename = NULL;
}

IntVariable::IntVariable(const char *varname)
    : Variable(varname, "int")
    , _value(0)
{
}

IntVariable::~IntVariable()
{
}

PyObject* IntVariable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }

    long val = *(int*)get();
    return PyLong_FromLong(val);
}

bool IntVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyLong_Check(obj)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyLongObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    int value = (int) PyLong_AsLong(obj);
    set(&value);

    return false;
}

UIntVariable::UIntVariable(const char *varname)
    : Variable(varname, "unsigned int")
    , _value(0)
{
}

UIntVariable::~UIntVariable()
{
}

PyObject* UIntVariable::getPythonObject(int i0, int n)
{
    unsigned long val = *(unsigned int*)get();
    return PyLong_FromUnsignedLong(val);
}

bool UIntVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyLong_Check(obj)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyLongObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    unsigned int value = (unsigned int) PyLong_AsLong(obj);
    set(&value);

    return false;
}

FloatVariable::FloatVariable(const char *varname)
    : Variable(varname, "float")
    , _value(0.f)
{
}

FloatVariable::~FloatVariable()
{
}

PyObject* FloatVariable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    double val = *(float*)get();
    return PyFloat_FromDouble(val);
}

bool FloatVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyFloat_Check(obj)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyFloatObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    float value = (float) PyFloat_AsDouble(obj);
    set(&value);

    return false;
}

Vec2Variable::Vec2Variable(const char *varname)
    : Variable(varname, "vec2")
{
    _value.x = 0.f;
    _value.y = 0.f;
}

Vec2Variable::~Vec2Variable()
{
}

PyObject* Vec2Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    vec2 *vv = (vec2*)get();
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, vv->s);
}

bool Vec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 2){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 2 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    vec2 *vv = (vec2*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(vec2));

    return false;
}

Vec3Variable::Vec3Variable(const char *varname)
    : Variable(varname, "vec3")
{
    _value.x = 0.f;
    _value.y = 0.f;
    _value.z = 0.f;
}

Vec3Variable::~Vec3Variable()
{
}

PyObject* Vec3Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    vec3 *vv = (vec3*)get();
    npy_intp dims[] = {3};
    return PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, vv->s);
}

bool Vec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 3){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 3 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    vec3 *vv = (vec3*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(vec3));

    return false;
}

Vec4Variable::Vec4Variable(const char *varname)
    : Variable(varname, "vec4")
{
    _value.x = 0.f;
    _value.y = 0.f;
    _value.z = 0.f;
    _value.w = 0.f;
}

Vec4Variable::~Vec4Variable()
{
}

PyObject* Vec4Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    vec4 *vv = (vec4*)get();
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, vv->s);
}

bool Vec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 4){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 4 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    vec4 *vv = (vec4*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(vec4));

    return false;
}

IVec2Variable::IVec2Variable(const char *varname)
    : Variable(varname, "ivec2")
{
    _value.x = 0;
    _value.y = 0;
}

IVec2Variable::~IVec2Variable()
{
}

PyObject* IVec2Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    ivec2 *vv = (ivec2*)get();
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_INT, vv->s);
}

bool IVec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 2){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 2 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    ivec2 *vv = (ivec2*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(ivec2));

    return false;
}

IVec3Variable::IVec3Variable(const char *varname)
    : Variable(varname, "ivec3")
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
}

IVec3Variable::~IVec3Variable()
{
}

PyObject* IVec3Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    ivec3 *vv = (ivec3*)get();
    npy_intp dims[] = {3};
    return PyArray_SimpleNewFromData(1, dims, PyArray_INT, vv->s);
}

bool IVec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 3){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 3 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    ivec3 *vv = (ivec3*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(ivec3));

    return false;
}

IVec4Variable::IVec4Variable(const char *varname)
    : Variable(varname, "ivec4")
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
    _value.w = 0;
}

IVec4Variable::~IVec4Variable()
{
}

PyObject* IVec4Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    ivec4 *vv = (ivec4*)get();
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_INT, vv->s);
}

bool IVec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 4){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 4 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    ivec4 *vv = (ivec4*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(ivec4));

    return false;
}

UIVec2Variable::UIVec2Variable(const char *varname)
    : Variable(varname, "uivec2")
{
    _value.x = 0;
    _value.y = 0;
}

UIVec2Variable::~UIVec2Variable()
{
}

PyObject* UIVec2Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    uivec2 *vv = (uivec2*)get();
    npy_intp dims[] = {2};
    return PyArray_SimpleNewFromData(1, dims, PyArray_UINT, vv->s);
}

bool UIVec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 2){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 2 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    uivec2 *vv = (uivec2*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(uivec2));

    return false;
}

UIVec3Variable::UIVec3Variable(const char *varname)
    : Variable(varname, "uivec3")
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
}

UIVec3Variable::~UIVec3Variable()
{
}

PyObject* UIVec3Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    uivec3 *vv = (uivec3*)get();
    npy_intp dims[] = {3};
    return PyArray_SimpleNewFromData(1, dims, PyArray_UINT, vv->s);
}

bool UIVec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 3){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 3 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    uivec3 *vv = (uivec3*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(uivec3));

    return false;
}

UIVec4Variable::UIVec4Variable(const char *varname)
    : Variable(varname, "uivec4")
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
    _value.w = 0;
}

UIVec4Variable::~UIVec4Variable()
{
}

PyObject* UIVec4Variable::getPythonObject(int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    uivec4 *vv = (uivec4*)get();
    npy_intp dims[] = {4};
    return PyArray_SimpleNewFromData(1, dims, PyArray_UINT, vv->s);
}

bool UIVec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"offset\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n != 0){
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", but \"n\" different from 0 has been received",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 1){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp dim = array_obj->dimensions[0];
    if(dim != 4){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected a 4 components array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    uivec4 *vv = (uivec4*)get();
    void *data = array_obj->data;
    memcpy(vv->s, data, sizeof(uivec4));

    return false;
}

ArrayVariable::ArrayVariable(const char *varname, const char *vartype)
    : Variable(varname, vartype)
    , _value(NULL)
{
}

ArrayVariable::~ArrayVariable()
{
    unsigned int i;
    for(i = 0; i < _objects.size(); i++){
        if(_objects.at(i))
            Py_DECREF(_objects.at(i));
        _objects.at(i) = NULL;
    }
    _objects.clear();
    for(i = 0; i < _data.size(); i++){
        if(_data.at(i))
            free(_data.at(i));
        _data.at(i) = NULL;
    }
    _data.clear();
    if(_value) clReleaseMemObject(_value); _value=NULL;
}

size_t ArrayVariable::size() const
{
    if(!_value)
        return 0;

    size_t memsize=0;
    cl_int status = clGetMemObjectInfo(_value,
                                       CL_MEM_SIZE,
                                       sizeof(size_t),
                                       &memsize,
                                       NULL);
    if(status != CL_SUCCESS){
        char msg[256];
        ScreenManager *S = ScreenManager::singleton();
        sprintf(msg,
                "Failure getting allocated memory from variable \"%s\"\n",
                name());
        S->addMessageF(3, msg);
        S->printOpenCLError(status);
    }
    return memsize;
}

PyObject* ArrayVariable::getPythonObject(int i0, int n)
{
    if(i0 < 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" cannot handle \"offset\" lower than 0",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if(n < 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" cannot handle \"n\" lower than 0",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	Variables *vars = C->variables();
	cl_int err_code;
    // Clear outdated references
    cleanMem();
    // Get the dimensions
    unsigned components = vars->typeToN(type());
    size_t typesize = vars->typeToBytes(type());
    size_t memsize = size();
    size_t offset = i0;
    if(offset * typesize > memsize){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure reading variable \"%s\" out of bounds",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    size_t len = memsize / typesize - offset;
    if(n != 0){
        len = n;
    }
    if(len == 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "0 bytes asked to be read from variable \"%s\"",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    if((offset + len) * typesize > memsize){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure reading variable \"%s\" out of bounds",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    npy_intp dims[] = {len, components};
    // Get the appropiate type
    int pytype = PyArray_FLOAT;
    if(strstr(type(), "unsigned int") ||
       strstr(type(), "uivec")){
       pytype = PyArray_UINT;
    }
    else if(strstr(type(), "int") ||
       strstr(type(), "ivec")){
       pytype = PyArray_INT;
    }
    else if(strstr(type(), "float") ||
       strstr(type(), "vec")){
       pytype = PyArray_FLOAT;
    }
    else{
        char errstr[128 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" is of type \"%s\", which is not handled by Python",
                name(),
                type());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    // Reallocate memory
    void *data = malloc(len * typesize);
    if(!data){
        char errstr[128 + strlen(name())];
        sprintf(errstr,
                "Failure allocating %lu bytes for variable \"%s\"",
                len * typesize,
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    _data.push_back(data);
    // Download the data
    err_code = clEnqueueReadBuffer(C->command_queue(),
                                   _value,
                                   CL_TRUE,
                                   offset * typesize,
                                   len * typesize,
                                   data,
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure downloading variable \"%s\"",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    // Build and return the Python object
    PyObject *obj = PyArray_SimpleNewFromData(2, dims, pytype, data);
    if(!obj){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure creating a Python object for variable \"%s\"",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }
    _objects.push_back(obj);
    Py_INCREF(obj);
    return obj;
}

bool ArrayVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 < 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" cannot handle \"offset\" lower than 0",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if(n < 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" cannot handle \"n\" lower than 0",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	Variables *vars = C->variables();
	cl_int err_code;
    // Clear outdated references
    cleanMem();
    // Get the dimensions
    unsigned components = vars->typeToN(type());
    size_t typesize = vars->typeToBytes(type());
    size_t memsize = size();
    size_t offset = i0;
    if(offset * typesize > memsize){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure writing variable \"%s\" out of bounds",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    size_t len = memsize / typesize - offset;
    if(n != 0){
        len = n;
    }
    if(len == 0){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "0 bytes asked to be written to variable \"%s\"",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if((offset + len) * typesize > memsize){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure writing variable \"%s\" out of bounds",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        char errstr[64 + strlen(name()) + strlen(type())];
        sprintf(errstr,
                "Variable \"%s\" expected a PyArrayObject",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 2){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Variable \"%s\" expected an one dimensional array",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    npy_intp *dims = array_obj->dimensions;
    if((size_t)dims[0] != len){
        char errstr[128 + strlen(name())];
        sprintf(errstr,
                "%lu elements have been asked to be written in variable \"%s\" but %lu have been provided",
                len,
                name(),
                (size_t)dims[0]);
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }
    if((size_t)dims[1] != components){
        char errstr[128 + strlen(name())];
        sprintf(errstr,
                "%lu components per elements are expected by variable \"%s\" but %lu have been provided",
                components,
                name(),
                (size_t)dims[1]);
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    void *data = array_obj->data;

    err_code =  clEnqueueWriteBuffer(C->command_queue(),
                                     _value,
                                     CL_TRUE,
                                     offset * typesize,
                                     len * typesize,
                                     data,
                                     0,
                                     NULL,
                                     NULL);
    if(err_code != CL_SUCCESS){
        char errstr[64 + strlen(name())];
        sprintf(errstr,
                "Failure uploading variable \"%s\"",
                name());
        PyErr_SetString(PyExc_ValueError, errstr);
        return true;
    }

    return false;
}

void ArrayVariable::cleanMem()
{
    int i;  // unsigned cannot be used here
    for(i = _objects.size() - 1; i >= 0; i--){
        if(_objects.at(i)->ob_refcnt == 1){
            Py_DECREF(_objects.at(i));
            free(_data.at(i));
            _data.at(i) = NULL;
            _data.erase(_data.begin() + i);
            _objects.erase(_objects.begin() + i);
        }
    }
}

// ---------------------------------------------------------------------------
// Variables manager
// ---------------------------------------------------------------------------

Variables::Variables()
{
}

Variables::~Variables()
{
    unsigned int i;
    for(i = 0; i < _vars.size(); i++){
        delete _vars.at(i);
    }
    _vars.clear();
}

bool Variables::registerVariable(const char* name,
                                 const char* type,
                                 const char* length,
                                 const char* value)
{
    // Look for an already existing variable with the same name
    unsigned int i;
    for(i = 0; i < _vars.size(); i++){
        if(!strcmp(_vars.at(i)->name(), name)){
            delete _vars.at(i);
            _vars.erase(_vars.begin() + i);
        }
    }

    // Discriminate scalar vs. array
    if(strstr(type, "*")){
        return registerClMem(name, type, length);
    }
    else{
        return registerScalar(name, type, value);
    }
    return false;
}

Variable* Variables::get(unsigned int index)
{
    if(index >= _vars.size()){
        return NULL;
    }
    return _vars.at(index);
}

Variable* Variables::get(const char* name)
{
    unsigned int i;
    for(i = 0; i < _vars.size(); i++){
        if(!strcmp(name, _vars.at(i)->name())){
            return _vars.at(i);
        }
    }
    return NULL;
}

size_t Variables::allocatedMemory(){
    unsigned int i;
    size_t allocated_mem = 0;
    for(i = 0; i < size(); i++){
        if(!strchr(_vars.at(i)->type(), '*')){
            continue;
        }
        ArrayVariable *var = (ArrayVariable *)_vars.at(i);
        allocated_mem += var->size();
    }
    return allocated_mem;
}

size_t Variables::typeToBytes(const char* type) const
{
    unsigned int n = typeToN(type);
    size_t type_size = 0;

    if(strstr(type, "unsigned int") ||
       strstr(type, "uivec")){
        type_size = sizeof(unsigned int);
    }
    else if(strstr(type, "int") ||
            strstr(type, "ivec")){
        type_size = sizeof(int);
    }
    else if(strstr(type, "float") ||
            strstr(type, "vec")){
        type_size = sizeof(float);
    }
    else{
        char msg[256];
        ScreenManager *S = ScreenManager::singleton();
        sprintf(msg,
                "Unvalid type \"%s\"\n",
                type);
        S->addMessageF(3, msg);
        return 0;
    }
    return n * type_size;
}

unsigned int Variables::typeToN(const char* type) const
{
    unsigned int n = 1;
    if(strstr(type, "vec2")) {
        n = 2;
    }
    else if(strstr(type, "vec3")) {
        n = 3;
    }
    else if(strstr(type, "vec4")) {
        n = 4;
    }
    else if(strstr(type, "vec")) {
        #ifdef HAVE_3D
            n = 4;
        #else
            n = 2;
        #endif // HAVE_3D
    }
    return n;
}

bool Variables::isSameType(const char* type_a,
                           const char* type_b,
                           bool ignore_asterisk)
{
    if(typeToN(type_a) != typeToN(type_b)){
        return false;
    }

    if(!ignore_asterisk){
        if(strchr(type_a, '*') && !strchr(type_b, '*')){
            return false;
        }
        else if(!strchr(type_a, '*') && strchr(type_b, '*')){
            return false;
        }
    }

    size_t len_a = strlen(type_a);
    if((type_a[len_a - 1] == '*')){
        len_a--;
    }
    if((type_a[len_a - 1] == '2') ||
       (type_a[len_a - 1] == '3') ||
       (type_a[len_a - 1] == '4')){
        len_a--;
    }
    size_t len_b = strlen(type_b);
    if((type_b[len_b - 1] == '*')){
        len_b--;
    }
    if((type_b[len_b - 1] == '2') ||
       (type_b[len_b - 1] == '3') ||
       (type_b[len_b - 1] == '4')){
        len_b--;
    }

    if(len_a != len_b){
        return false;
    }

    if(strncmp(type_a, type_b, len_a)){
        return false;
    }

    return true;
}

bool Variables::solve(const char *type_name,
                      const char *value,
                      void *data,
                      const char* name)
{
    char msg[256];
    ScreenManager *S = ScreenManager::singleton();
    size_t typesize = typeToBytes(type_name);
    if(!typesize){
        return true;
    }
    if(!strcmp(value, "")){
        S->addMessageF(3, "Empty value received\n");
        return 0;
    }

    char *type = new char[strlen(type_name) + 1];
    strcpy(type, type_name);
    if(strchr(type, '*'))
        strcpy(strchr(type, '*'), "");

    if(!strcmp(type, "int")){
        int val;
        float auxval;
        if(readComponents(name, value, 1, &auxval))
            return true;
        val = round(auxval);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "unsigned int")){
        unsigned int val;
        float auxval;
        if(readComponents(name, value, 1, &auxval))
            return true;
        val = (unsigned int)round(auxval);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "float")){
        float val;
        if(readComponents(name, value, 1, &val))
            return true;
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "vec")){
        vec val;
        #ifdef HAVE_3D
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            val.w = auxval[3];
            memcpy(data, &val, typesize);
        #else
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
        #endif
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "vec2")){
        vec2 val;
        float auxval[2];
        if(readComponents(name, value, 2, auxval))
            return true;
        val.x = auxval[0];
        val.y = auxval[1];
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "vec3")){
        vec3 val;
        float auxval[3];
        if(readComponents(name, value, 3, auxval))
            return true;
        val.x = auxval[0];
        val.y = auxval[1];
        val.z = auxval[2];
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "vec4")){
        vec4 val;
        float auxval[4];
        if(readComponents(name, value, 4, auxval))
            return true;
        val.x = auxval[0];
        val.y = auxval[1];
        val.z = auxval[2];
        val.w = auxval[3];
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "ivec")){
        ivec val;
        #ifdef HAVE_3D
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            val.z = round(auxval[2]);
            val.w = round(auxval[3]);
            memcpy(data, &val, typesize);
        #else
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
        #endif
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "ivec2")){
        ivec2 val;
        float auxval[2];
        if(readComponents(name, value, 2, auxval))
            return true;
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "ivec3")){
        ivec3 val;
        float auxval[3];
        if(readComponents(name, value, 3, auxval))
            return true;
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        val.z = round(auxval[2]);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "ivec4")){
        ivec4 val;
        float auxval[4];
        if(readComponents(name, value, 4, auxval))
            return true;
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        val.z = round(auxval[2]);
        val.w = round(auxval[3]);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "uivec")){
        uivec val;
        #ifdef HAVE_3D
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            val.z = (unsigned int)round(auxval[2]);
            val.w = (unsigned int)round(auxval[3]);
            memcpy(data, &val, typesize);
        #else
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
        #endif
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "uivec2")){
        uivec2 val;
        float auxval[2];
        if(readComponents(name, value, 2, auxval))
            return true;
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "uivec3")){
        uivec3 val;
        float auxval[3];
        if(readComponents(name, value, 3, auxval))
            return true;
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        val.z = (unsigned int)round(auxval[2]);
        memcpy(data, &val, typesize);
    }
    else if(!strcmp(type, "uivec4")){
        uivec4 val;
        float auxval[4];
        if(readComponents(name, value, 4, auxval))
            return true;
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        val.z = (unsigned int)round(auxval[2]);
        val.w = (unsigned int)round(auxval[3]);
        memcpy(data, &val, typesize);
    }
    else{
        return true;
    }

    delete[] type; type = NULL;
    return false;
}

bool Variables::registerScalar(const char* name,
                               const char* type,
                               const char* value)
{
    if(!strcmp(type, "int")){
        IntVariable *var = new IntVariable(name);
        if(strcmp(value, "")){
            int val = round(tok.solve(value));
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "unsigned int")){
        UIntVariable *var = new UIntVariable(name);
        if(strcmp(value, "")){
            unsigned int val = (unsigned int)round(tok.solve(value));
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "float")){
        FloatVariable *var = new FloatVariable(name);
        if(strcmp(value, "")){
            float val = tok.solve(value);
            tok.registerVariable(name, val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec")){
        VecVariable *var = new VecVariable(name);
        if(strcmp(value, "")){
            vec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = auxval[0];
                val.y = auxval[1];
                val.z = auxval[2];
                val.w = auxval[3];
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = auxval[0];
                val.y = auxval[1];
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec2")){
        Vec2Variable *var = new Vec2Variable(name);
        if(strcmp(value, "")){
            vec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec3")){
        Vec3Variable *var = new Vec3Variable(name);
        if(strcmp(value, "")){
            vec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec4")){
        Vec4Variable *var = new Vec4Variable(name);
        if(strcmp(value, "")){
            vec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            val.w = auxval[3];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec")){
        IVecVariable *var = new IVecVariable(name);
        if(strcmp(value, "")){
            ivec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = round(auxval[0]);
                val.y = round(auxval[1]);
                val.z = round(auxval[2]);
                val.w = round(auxval[3]);
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = round(auxval[0]);
                val.y = round(auxval[1]);
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec2")){
        IVec2Variable *var = new IVec2Variable(name);
        if(strcmp(value, "")){
            ivec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec3")){
        IVec3Variable *var = new IVec3Variable(name);
        if(strcmp(value, "")){
            ivec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            val.z = round(auxval[2]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec4")){
        IVec4Variable *var = new IVec4Variable(name);
        if(strcmp(value, "")){
            ivec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            val.z = round(auxval[2]);
            val.w = round(auxval[3]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec")){
        UIVecVariable *var = new UIVecVariable(name);
        if(strcmp(value, "")){
            uivec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = (unsigned int)round(auxval[0]);
                val.y = (unsigned int)round(auxval[1]);
                val.z = (unsigned int)round(auxval[2]);
                val.w = (unsigned int)round(auxval[3]);
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = (unsigned int)round(auxval[0]);
                val.y = (unsigned int)round(auxval[1]);
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec2")){
        UIVec2Variable *var = new UIVec2Variable(name);
        if(strcmp(value, "")){
            uivec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec3")){
        UIVec3Variable *var = new UIVec3Variable(name);
        if(strcmp(value, "")){
            uivec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            val.z = (unsigned int)round(auxval[2]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec4")){
        UIVec4Variable *var = new UIVec4Variable(name);
        if(strcmp(value, "")){
            uivec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            val.z = (unsigned int)round(auxval[2]);
            val.w = (unsigned int)round(auxval[3]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else{
        char msg[256];
        ScreenManager *S = ScreenManager::singleton();
        sprintf(msg,
                "\"%s\" declared as \"%s\", which is not a valid scalar type.\n",
                name,
                type);
        S->addMessageF(3, msg);
        S->addMessageF(0, "Valid types are:\n");
        S->addMessageF(0, "\tint\n");
        S->addMessageF(0, "\tunsigned int\n");
        S->addMessageF(0, "\tfloat\n");
        S->addMessageF(0, "\tvec\n");
        S->addMessageF(0, "\tvec2\n");
        S->addMessageF(0, "\tvec3\n");
        S->addMessageF(0, "\tvec4\n");
        S->addMessageF(0, "\tivec\n");
        S->addMessageF(0, "\tivec2\n");
        S->addMessageF(0, "\tivec3\n");
        S->addMessageF(0, "\tivec4\n");
        S->addMessageF(0, "\tuivec\n");
        S->addMessageF(0, "\tuivec2\n");
        S->addMessageF(0, "\tuivec3\n");
        S->addMessageF(0, "\tuivec4\n");
        return true;
    }
    return false;
}

bool Variables::registerClMem(const char* name,
                              const char* type,
                              const char* length)
{
    size_t typesize;
    unsigned int n;
    char msg[256];
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    ScreenManager *S = ScreenManager::singleton();
    // Get the type size
    char auxtype[strlen(type) + 1];
    strcpy(auxtype, type);
    strcpy(strchr(auxtype, '*'), "");
    if(!strcmp(auxtype, "int")){
        typesize = sizeof(int);
    }
    else if(!strcmp(auxtype, "unsigned int")){
        typesize = sizeof(unsigned int);
    }
    else if(!strcmp(auxtype, "float")){
        typesize = sizeof(float);
    }
    else if(!strcmp(auxtype, "vec")){
        typesize = sizeof(vec);
    }
    else if(!strcmp(auxtype, "vec2")){
        typesize = sizeof(vec2);
    }
    else if(!strcmp(auxtype, "vec3")){
        typesize = sizeof(vec3);
    }
    else if(!strcmp(auxtype, "vec4")){
        typesize = sizeof(vec4);
    }
    else if(!strcmp(auxtype, "ivec")){
        typesize = sizeof(ivec);
    }
    else if(!strcmp(auxtype, "ivec2")){
        typesize = sizeof(ivec2);
    }
    else if(!strcmp(auxtype, "ivec3")){
        typesize = sizeof(ivec3);
    }
    else if(!strcmp(auxtype, "ivec4")){
        typesize = sizeof(ivec4);
    }
    else if(!strcmp(auxtype, "uivec")){
        typesize = sizeof(uivec);
    }
    else if(!strcmp(auxtype, "uivec2")){
        typesize = sizeof(uivec2);
    }
    else if(!strcmp(auxtype, "uivec3")){
        typesize = sizeof(uivec3);
    }
    else if(!strcmp(auxtype, "uivec4")){
        typesize = sizeof(uivec4);
    }
    else{
        sprintf(msg,
                "\"%s\" declared as \"%s\", which is not a valid array type.\n",
                name,
                type);
        S->addMessageF(3, msg);
        S->addMessageF(0, "Valid types are:\n");
        S->addMessageF(0, "\tint*\n");
        S->addMessageF(0, "\tunsigned int*\n");
        S->addMessageF(0, "\tfloat*\n");
        S->addMessageF(0, "\tvec*\n");
        S->addMessageF(0, "\tvec2*\n");
        S->addMessageF(0, "\tvec3*\n");
        S->addMessageF(0, "\tvec4*\n");
        S->addMessageF(0, "\tivec*\n");
        S->addMessageF(0, "\tivec2*\n");
        S->addMessageF(0, "\tivec3*\n");
        S->addMessageF(0, "\tivec4*\n");
        S->addMessageF(0, "\tuivec*\n");
        S->addMessageF(0, "\tuivec2*\n");
        S->addMessageF(0, "\tuivec3*\n");
        S->addMessageF(0, "\tuivec4*\n");
        return true;
    }

    // Get the length
    n = 0;
    if(strcmp(length, ""))
        n = (unsigned int)round(tok.solve(length));

    // Generate the variable
    ArrayVariable *var = new ArrayVariable(name, type);
    if(n){
        // Allocate memory on device
        cl_int status;
        cl_mem mem;
        mem = clCreateBuffer(C->context(),
                             CL_MEM_READ_WRITE,
                             n * typesize,
                             NULL,
                             &status);
        if(status != CL_SUCCESS) {
            S->addMessageF(3, "Allocation failure.\n");
            S->printOpenCLError(status);
            return true;
        }
        var->set(&mem);
    }
    _vars.push_back(var);

    return false;
}

bool Variables::readComponents(const char* name,
                               const char* value,
                               unsigned int n,
                               float* v)
{
    float val;
    unsigned int i;
    char msg[256];
    bool error;
    ScreenManager *S = ScreenManager::singleton();
    if(n == 0){
        sprintf(msg,
                "%u components required for the variable \"%s\".\n",
                n,
                name);
        S->addMessageF(3, msg);
    }
    if(n > 4){
        S->addMessageF(3, "No more than 4 components can be required\n");
        sprintf(msg,
                "%u components required for the variable \"%s\".\n",
                n,
                name);
        S->addMessageF(0, msg);
    }
    char* remain = (char *)value;
    char aux[strlen(value) + 1];
    char nameaux[strlen(name) + 3];
    const char* extensions[4] = {"_x", "_y", "_z", "_w"};
    for(i = 0; i < n - 1; i++){
        strcpy(aux, remain);
        if(!strchr(aux, ',')){
            sprintf(msg, "Failure reading the variable \"%s\" value", name);
            S->addMessageF(3, msg);
            sprintf(msg, "%u fields expected, %u received.\n", n, i);
            S->addMessageF(0, msg);
            return true;
        }
        strcpy(strchr(aux, ','), "");
        remain = strchr(remain, ',') + 1;
        val = tok.solve(aux, &error);
        if(error)
            return true;
        strcpy(nameaux, name);
        strcat(nameaux, extensions[i]);
        tok.registerVariable(nameaux, val);
        v[i] = val;
    }

    strcpy(aux, remain);
    if(strchr(aux, ','))
        strcpy(strchr(aux, ','), "");
    val = tok.solve(aux, &error);
    if(error)
        return true;
    strcpy(nameaux, name);
    if(n != 1)
        strcat(nameaux, extensions[n - 1]);
    tok.registerVariable(nameaux, val);
    v[n - 1] = val;
    return false;
}

}}  // namespace
