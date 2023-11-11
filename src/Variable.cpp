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
 * @brief Virtual variables environment to allow the user define/manipulate
 * the variables used in the simulation externally.
 * (see Aqua::InpuOutput::Variable and Aqua::InpuOutput::Variables)
 */

#include <Variable.h>
#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <algorithm>

/** @def PY_ARRAY_UNIQUE_SYMBOL
 * @brief Define the extension module which this Python stuff should be linked
 * to.
 *
 * In AQUAgpusph all the Python stuff is linked in the same group AQUA_ARRAY_API
 * @see
 * http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define PY_ARRAY_UNIQUE_SYMBOL AQUA_ARRAY_API
/** @def NO_IMPORT_ARRAY
 * @brief Set this file as a helper of the group AQUA_ARRAY_API.
 * @see
 * http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define NO_IMPORT_ARRAY
#include <numpy/npy_no_deprecated_api.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

namespace Aqua {
namespace InputOutput {

static std::string str_val;
static std::ostringstream pyerr;

Variable::Variable(const std::string varname, const std::string vartype)
  : _name(varname)
  , _typename(vartype)
  , _event(NULL)
  , _synced(true)
  , _synced_for_read(true)
{
	cl_int err_code;
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	// Create a dummy starting event for the queue
	_event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure creating user event for \"" << varname << "\" variable."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clSetUserEventStatus(_event, CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting user event status for \"" << varname
		    << "\" variable." << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
}

void
Variable::setEvent(cl_event event)
{
	cl_int err_code;
	if (_event != NULL) {
		// Forgive the former/predecessor event
		err_code = clReleaseEvent(_event);
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure releasing the predecessor event for \"" << name()
			    << "\" variable." << std::endl;
			LOG(L_ERROR, msg.str());
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL execution error");
		}
	}
	// And get retained the current event
	err_code = clRetainEvent(event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure reteaning the writing event for \"" << name()
		    << "\" variable." << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	_event = event;

	// Output events are blocking all reading events, so we can forget about
	// them
	for (auto e : _reading_events) {
		err_code = clReleaseEvent(e);
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure releasing a reading event for \"" << name()
			    << "\" variable." << std::endl;
			LOG(L_ERROR, msg.str());
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL execution error");
		}
	}
	_reading_events.clear();
	_synced = false;
	_synced_for_read = false;
}

void
Variable::addReadingEvent(cl_event event)
{
	cl_int err_code;
	// Tidy up the former events
	cleanReadingEvents();
	// And get retained the current event
	err_code = clRetainEvent(event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure reteaning the reading event for \"" << name()
		    << "\" variable." << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	_reading_events.push_back(event);
	_synced = false;
}

void
Variable::sync(bool readonly)
{
	if (_synced || (readonly && _synced_for_read))
		return;

	cl_int err_code;
	std::vector<cl_event> events;
	if (!readonly)
		events = _reading_events;
	if (_event != NULL)
		events.push_back(_event);

	if (!events.size()) {
		if (readonly) {
			_synced_for_read = true;
			if (!_reading_events.size())
				_synced = true;
		} else {
			_synced = true;
		}
		return;
	}
	err_code = clWaitForEvents(events.size(), events.data());
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure syncing variable \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	if (readonly) {
		_synced_for_read = true;
		if (!_reading_events.size())
			_synced = true;
	} else {
		_synced = true;
	}
}

void
Variable::cleanReadingEvents()
{
	cl_int err_code, event_status;
	for (auto it = _reading_events.begin(); it != _reading_events.end();) {
		cl_event event = *it;
		err_code = clGetEventInfo(event,
		                          CL_EVENT_COMMAND_EXECUTION_STATUS,
		                          sizeof(cl_int),
		                          &event_status,
		                          NULL);
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure querying a reading event status for variable \""
			    << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL execution error");
		}

		if (event_status == CL_COMPLETE) {
			err_code = clReleaseEvent(event);
			if (err_code != CL_SUCCESS) {
				std::stringstream msg;
				msg << "Failure releasing a reading event for \"" << name()
				    << "\" variable." << std::endl;
				LOG(L_ERROR, msg.str());
				Logger::singleton()->printOpenCLError(err_code);
				throw std::runtime_error("OpenCL execution error");
			}
			it = _reading_events.erase(it);
		} else
			++it;
	}
}

template<class T>
const std::string
ScalarNumberVariable<T>::asString()
{
	std::ostringstream msg;
	msg << *((T*)this->get());
	str_val = msg.str();
	return str_val;
}

IntVariable::IntVariable(const std::string varname)
  : ScalarNumberVariable<int>(varname, "int")
{
	_value = 0;
}

PyObject*
IntVariable::getPythonObject(int i0, int n)
{
	long val = *((int*)get());
	return PyLong_FromLong(val);
}

bool
IntVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyLong_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyLongObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	int value = (int)PyLong_AsLong(obj);
	set(&value);

	return false;
}

UIntVariable::UIntVariable(const std::string varname)
  : ScalarNumberVariable<unsigned int>(varname, "unsigned int")
{
	_value = 0u;
}

PyObject*
UIntVariable::getPythonObject(int i0, int n)
{
	unsigned long val = *(unsigned int*)get();
	return PyLong_FromUnsignedLong(val);
}

bool
UIntVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyLong_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyLongObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	unsigned int value = (unsigned int)PyLong_AsLong(obj);
	set(&value);

	return false;
}

FloatVariable::FloatVariable(const std::string varname)
  : ScalarNumberVariable<float>(varname, "float")
{
	_value = 0.f;
}

PyObject*
FloatVariable::getPythonObject(int i0, int n)
{
	double val = *((float*)get());
	return PyFloat_FromDouble(val);
}

bool
FloatVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyFloat_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyFloatObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	float value = (float)PyFloat_AsDouble(obj);
	set(&value);

	return false;
}

template<class T>
ScalarVecVariable<T>::ScalarVecVariable(const std::string varname,
                                        const std::string vartype,
                                        const unsigned int dims)
  : ScalarVariable<T>(varname, vartype)
  , _dims(dims)
{
}

template<class T>
bool
ScalarVecVariable<T>::checkPyhonObjectDims(PyObject* obj)
{
	if (!PyObject_TypeCheck(obj, &PyArray_Type)) {
		pyerr.str("");
		pyerr << "Variable \"" << this->name() << "\" expected a PyArrayObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	PyArrayObject* array_obj = (PyArrayObject*)obj;
	if (PyArray_NDIM(array_obj) != 1) {
		pyerr.str("");
		pyerr << "Variable \"" << this->name()
		      << "\" expected an one dimensional array" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	npy_intp dim = PyArray_DIMS(array_obj)[0];
	if (dim != _dims) {
		pyerr.str("");
		pyerr << "Variable \"" << this->name() << "\" expected an array with "
		      << _dims << " components" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	return false;
}

template<class T>
const std::string
ScalarVecVariable<T>::asString()
{
	T* val = (T*)(this->get());
	std::ostringstream msg;
	msg << "(";
	for (unsigned int i = 0; i < _dims; i++) {
		msg << val->s[i] << ",";
	}
	str_val = msg.str();
	str_val.back() = ')';
	return str_val;
}

Vec2Variable::Vec2Variable(const std::string varname)
  : ScalarVecVariable(varname, "vec2", 2)
{
	_value.x = 0.f;
	_value.y = 0.f;
}

PyObject*
Vec2Variable::getPythonObject(int i0, int n)
{
	vec2* vv = (vec2*)get();
	npy_intp dims[] = { 2 };
	return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, vv->s);
}

bool
Vec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_FLOAT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_FLOAT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	vec2* vv = (vec2*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(vec2));

	return false;
}

Vec3Variable::Vec3Variable(const std::string varname)
  : ScalarVecVariable(varname, "vec3", 3)
{
	_value.x = 0.f;
	_value.y = 0.f;
	_value.z = 0.f;
}

PyObject*
Vec3Variable::getPythonObject(int i0, int n)
{
	vec3* vv = (vec3*)get();
	npy_intp dims[] = { 3 };
	return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, vv->s);
}

bool
Vec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_FLOAT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_FLOAT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	vec3* vv = (vec3*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(vec3));

	return false;
}

Vec4Variable::Vec4Variable(const std::string varname)
  : ScalarVecVariable(varname, "vec4", 4)
{
	_value.x = 0.f;
	_value.y = 0.f;
	_value.z = 0.f;
	_value.w = 0.f;
}

PyObject*
Vec4Variable::getPythonObject(int i0, int n)
{
	vec4* vv = (vec4*)get();
	npy_intp dims[] = { 4 };
	return PyArray_SimpleNewFromData(1, dims, NPY_FLOAT32, vv->s);
}

bool
Vec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_FLOAT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_FLOAT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	vec4* vv = (vec4*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(vec4));

	return false;
}

IVec2Variable::IVec2Variable(const std::string varname)
  : ScalarVecVariable(varname, "ivec2", 2)
{
	_value.x = 0;
	_value.y = 0;
}

PyObject*
IVec2Variable::getPythonObject(int i0, int n)
{
	ivec2* vv = (ivec2*)get();
	npy_intp dims[] = { 2 };
	return PyArray_SimpleNewFromData(1, dims, NPY_INT32, vv->s);
}

bool
IVec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_INT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_INT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	ivec2* vv = (ivec2*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(ivec2));

	return false;
}

IVec3Variable::IVec3Variable(const std::string varname)
  : ScalarVecVariable(varname, "ivec3", 3)
{
	_value.x = 0;
	_value.y = 0;
	_value.z = 0;
}

PyObject*
IVec3Variable::getPythonObject(int i0, int n)
{
	ivec3* vv = (ivec3*)get();
	npy_intp dims[] = { 3 };
	return PyArray_SimpleNewFromData(1, dims, NPY_INT32, vv->s);
}

bool
IVec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_INT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_INT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	ivec3* vv = (ivec3*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(ivec3));

	return false;
}

IVec4Variable::IVec4Variable(const std::string varname)
  : ScalarVecVariable(varname, "ivec4", 4)
{
	_value.x = 0;
	_value.y = 0;
	_value.z = 0;
	_value.w = 0;
}

PyObject*
IVec4Variable::getPythonObject(int i0, int n)
{
	ivec4* vv = (ivec4*)get();
	npy_intp dims[] = { 4 };
	return PyArray_SimpleNewFromData(1, dims, NPY_INT32, vv->s);
}

bool
IVec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_INT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_INT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	ivec4* vv = (ivec4*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(ivec4));

	return false;
}

UIVec2Variable::UIVec2Variable(const std::string varname)
  : ScalarVecVariable(varname, "uivec2", 2)
{
	_value.x = 0;
	_value.y = 0;
}

PyObject*
UIVec2Variable::getPythonObject(int i0, int n)
{
	uivec2* vv = (uivec2*)get();
	npy_intp dims[] = { 2 };
	return PyArray_SimpleNewFromData(1, dims, NPY_UINT32, vv->s);
}

bool
UIVec2Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_UINT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_UINT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	uivec2* vv = (uivec2*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(uivec2));

	return false;
}

UIVec3Variable::UIVec3Variable(const std::string varname)
  : ScalarVecVariable(varname, "uivec3", 3)
{
	_value.x = 0;
	_value.y = 0;
	_value.z = 0;
}

PyObject*
UIVec3Variable::getPythonObject(int i0, int n)
{
	uivec3* vv = (uivec3*)get();
	npy_intp dims[] = { 3 };
	return PyArray_SimpleNewFromData(1, dims, NPY_UINT32, vv->s);
}

bool
UIVec3Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_UINT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_UINT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	uivec3* vv = (uivec3*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(uivec3));

	return false;
}

UIVec4Variable::UIVec4Variable(const std::string varname)
  : ScalarVecVariable(varname, "uivec4", 4)
{
	_value.x = 0;
	_value.y = 0;
	_value.z = 0;
	_value.w = 0;
}

PyObject*
UIVec4Variable::getPythonObject(int i0, int n)
{
	uivec4* vv = (uivec4*)get();
	npy_intp dims[] = { 4 };
	return PyArray_SimpleNewFromData(1, dims, NPY_UINT32, vv->s);
}

bool
UIVec4Variable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	int typ = PyArray_TYPE(array_obj);
	if (typ != NPY_UINT32) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a NPY_UINT32 array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	uivec4* vv = (uivec4*)get();
	void* data = PyArray_DATA(array_obj);
	memcpy(vv->s, data, sizeof(uivec4));

	return false;
}

ArrayVariable::ArrayVariable(const std::string varname,
                             const std::string vartype)
  : Variable(varname, vartype)
  , _value(NULL)
{
}

ArrayVariable::~ArrayVariable()
{
	for (auto object : _objects) {
		if (object)
			Py_DECREF(object);
	}
	_objects.clear();
	for (auto data : _data) {
		if (data)
			free(data);
	}
	_data.clear();
	if (_value)
		clReleaseMemObject(_value);
	_value = NULL;
}

inline bool
ArrayVariable::isArray()
{
	return true;
}

size_t
ArrayVariable::size() const
{
	if (!_value)
		return 0;

	size_t memsize = 0;
	cl_int status =
	    clGetMemObjectInfo(_value, CL_MEM_SIZE, sizeof(size_t), &memsize, NULL);
	if (status != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure getting allocated memory from variable \"" << name()
		    << "\"." << std::endl,
		    LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(status);
	}
	return memsize;
}

PyObject*
ArrayVariable::getPythonObject(int i0, int n)
{
	if (i0 < 0) {
		pyerr.str("");
		pyerr << "Variable \"" << name()
		      << "\" cannot handle \"offset\" lower than 0" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	if (n < 0) {
		pyerr.str("");
		pyerr << "Variable \"" << name()
		      << "\" cannot handle \"n\" lower than 0" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	Variables* vars = C->variables();
	cl_int err_code;
	// Clear outdated references
	cleanMem();
	// Get the dimensions
	unsigned components = vars->typeToN(type());
	size_t typesize = vars->typeToBytes(type());
	size_t memsize = size();
	size_t offset = i0;
	if (offset * typesize > memsize) {
		pyerr.str("");
		pyerr << "Failure reading variable \"" << name() << "\" out of bounds"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	size_t len = memsize / typesize - offset;
	if (n != 0) {
		len = n;
	}
	if (len == 0) {
		pyerr.str("");
		pyerr << "0 bytes asked to be read from variable \"" << name() << "\""
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	if ((offset + len) * typesize > memsize) {
		pyerr.str("");
		pyerr << "Failure reading variable \"" << name() << "\" out of bounds"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	npy_intp dims[] = { static_cast<npy_intp>(len), components };
	// Get the appropiate type
	int pytype = NPY_FLOAT32;
	if (!type().compare("unsigned int") || !type().compare("unsigned int*") ||
	    !type().compare("uivec") || !type().compare("uivec*")) {
		pytype = NPY_UINT32;
	} else if (!type().compare("int") || !type().compare("int*") ||
	           !type().compare("ivec") || !type().compare("ivec*")) {
		pytype = NPY_INT32;
	} else if (!type().compare("float") || !type().compare("float*") ||
	           !type().compare("vec") || !type().compare("vec*")) {
		pytype = NPY_FLOAT32;
	} else {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" is of type \"" << type()
		      << "\", which can't be handled by Python" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	// Reallocate memory
	void* data = malloc(len * typesize);
	if (!data) {
		pyerr.str("");
		pyerr << "Failure allocating " << len * typesize
		      << " bytes for variable \"" << name() << "\"" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	_data.push_back(data);
	// Download the data
	cl_event event_wait = getEvent();
	err_code = clEnqueueReadBuffer(C->command_queue(),
	                               _value,
	                               CL_TRUE,
	                               offset * typesize,
	                               len * typesize,
	                               data,
	                               1,
	                               &event_wait,
	                               NULL);
	if (err_code != CL_SUCCESS) {
		pyerr.str("");
		pyerr << "Failure downloading variable \"" << name() << "\""
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	// Build and return the Python object
	PyObject* obj = PyArray_SimpleNewFromData(2, dims, pytype, data);
	if (!obj) {
		pyerr.str("");
		pyerr << "Failure creating a Python object for variable \"" << name()
		      << "\"" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return NULL;
	}
	_objects.push_back(obj);
	Py_INCREF(obj);
	return obj;
}

bool
ArrayVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (i0 < 0) {
		pyerr.str("");
		pyerr << "Variable \"" << name()
		      << "\" cannot handle \"offset\" lower than 0" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	if (n < 0) {
		pyerr.str("");
		pyerr << "Variable \"" << name()
		      << "\" cannot handle \"n\" lower than 0" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	Variables* vars = C->variables();
	cl_int err_code;
	// Clear outdated references
	cleanMem();
	// Get the dimensions
	unsigned components = vars->typeToN(type());
	size_t typesize = vars->typeToBytes(type());
	size_t memsize = size();
	size_t offset = i0;
	if (offset * typesize > memsize) {
		pyerr.str("");
		pyerr << "Failure writing variable \"" << name() << "\" out of bounds"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	size_t len = memsize / typesize - offset;
	if (n != 0) {
		len = n;
	}
	if (len == 0) {
		pyerr.str("");
		pyerr << "0 bytes asked to be written into variable \"" << name()
		      << "\"" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	if ((offset + len) * typesize > memsize) {
		pyerr.str("");
		pyerr << "Failure writing variable \"" << name() << "\" out of bounds"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	if (!PyObject_TypeCheck(obj, &PyArray_Type)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyArrayObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	PyArrayObject* array_obj = (PyArrayObject*)obj;
	if (PyArray_NDIM(array_obj) != 2) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a 2D array"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	npy_intp* dims = PyArray_DIMS(array_obj);
	if ((size_t)dims[0] != len) {
		pyerr.str("");
		pyerr << len << " elements have been asked to be written in variable \""
		      << name() << "\" but " << (size_t)dims[0] << " have been provided"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	if ((size_t)dims[1] != components) {
		pyerr.str("");
		pyerr << components
		      << " components per elements are expected by variable \""
		      << name() << "\" but " << (size_t)dims[1] << " have been provided"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	void* data = PyArray_DATA(array_obj);
	cl_event event, event_wait = getEvent();
	err_code = clEnqueueWriteBuffer(C->command_queue(),
	                                _value,
	                                CL_FALSE,
	                                offset * typesize,
	                                len * typesize,
	                                data,
	                                1,
	                                &event_wait,
	                                &event);
	if (err_code != CL_SUCCESS) {
		pyerr.str("");
		pyerr << "Failure uploading variable \"" << name() << "\"" << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	setEvent(event);
	err_code = clReleaseEvent(event);
	if (err_code != CL_SUCCESS) {
		pyerr.str("");
		pyerr << "Failure releasing variable \"" << name() << "\" event"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	return false;
}

const std::string
ArrayVariable::asString()
{
	cl_mem* val = (cl_mem*)get();
	std::ostringstream msg;
	msg << val;
	str_val = msg.str();
	return str_val;
}

const std::string
ArrayVariable::asString(size_t i)
{
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	size_t type_size = Variables::typeToBytes(type());
	size_t length = size() / type_size;
	if (i > length) {
		std::ostringstream msg;
		msg << "Out of bounds (length = " << length << ")" << std::endl;
		LOG(L_ERROR, msg.str());
		return NULL;
	}
	void* ptr = malloc(type_size);
	if (!ptr) {
		std::ostringstream msg;
		msg << "Failure allocating memory to download the variable \"" << name()
		    << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << type_size << "bytes requested" << std::endl;
		LOG0(L_DEBUG, msg.str());
		return NULL;
	}
	cl_event event_wait = getEvent();
	cl_int err_code = clEnqueueReadBuffer(C->command_queue(),
	                                      _value,
	                                      CL_TRUE,
	                                      i * type_size,
	                                      type_size,
	                                      ptr,
	                                      1,
	                                      &event_wait,
	                                      NULL);
	if (err_code != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure downloading the variable \"" << name() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		free(ptr);
		return NULL;
	}

	std::ostringstream str_stream;

	if (!type().compare("unsigned int*")) {
		str_stream << ((unsigned int*)ptr)[0];
	} else if (!type().compare("uivec2*")) {
		str_stream << "(" << ((unsigned int*)ptr)[0] << ","
		           << ((unsigned int*)ptr)[1] << ")";
	} else if (!type().compare("uivec3*")) {
		str_stream << "(" << ((unsigned int*)ptr)[0] << ","
		           << ((unsigned int*)ptr)[1] << "," << ((unsigned int*)ptr)[2]
		           << ")";
	} else if (!type().compare("uivec4*")) {
		str_stream << "(" << ((unsigned int*)ptr)[0] << ","
		           << ((unsigned int*)ptr)[1] << "," << ((unsigned int*)ptr)[2]
		           << "," << ((unsigned int*)ptr)[3] << ")";
	} else if (!type().compare("uivec*")) {
#ifdef HAVE_3D
		str_stream << "(" << ((unsigned int*)ptr)[0] << ","
		           << ((unsigned int*)ptr)[1] << "," << ((unsigned int*)ptr)[2]
		           << "," << ((unsigned int*)ptr)[3] << ")";
#else
		str_stream << "(" << ((unsigned int*)ptr)[0] << ","
		           << ((unsigned int*)ptr)[1] << ")";
#endif
	} else if (!type().compare("int*")) {
		str_stream << ((int*)ptr)[0];
	} else if (!type().compare("ivec2*")) {
		str_stream << "(" << ((int*)ptr)[0] << "," << ((int*)ptr)[1] << ")";
	} else if (!type().compare("ivec3*")) {
		str_stream << "(" << ((int*)ptr)[0] << "," << ((int*)ptr)[1] << ","
		           << ((int*)ptr)[2] << ")";
	} else if (!type().compare("ivec4*")) {
		str_stream << "(" << ((int*)ptr)[0] << "," << ((int*)ptr)[1] << ","
		           << ((int*)ptr)[2] << "," << ((int*)ptr)[3] << ")";
	} else if (!type().compare("ivec*")) {
#ifdef HAVE_3D
		str_stream << "(" << ((int*)ptr)[0] << "," << ((int*)ptr)[1] << ","
		           << ((int*)ptr)[2] << "," << ((int*)ptr)[3] << ")";
#else
		str_stream << "(" << ((int*)ptr)[0] << "," << ((int*)ptr)[1] << ")";
#endif
	} else if (!type().compare("float*")) {
		str_stream << ((float*)ptr)[0];
	} else if (!type().compare("vec2*")) {
		str_stream << "(" << ((float*)ptr)[0] << "," << ((float*)ptr)[1] << ")";
	} else if (!type().compare("vec3*")) {
		str_stream << "(" << ((float*)ptr)[0] << "," << ((float*)ptr)[1] << ","
		           << ((float*)ptr)[2] << ")";
	} else if (!type().compare("vec4*")) {
		str_stream << "(" << ((float*)ptr)[0] << "," << ((float*)ptr)[1] << ","
		           << ((float*)ptr)[2] << "," << ((float*)ptr)[3] << ")";
	} else if (!type().compare("vec*")) {
#ifdef HAVE_3D
		str_stream << "(" << ((float*)ptr)[0] << "," << ((float*)ptr)[1] << ","
		           << ((float*)ptr)[2] << "," << ((float*)ptr)[3] << ")";
#else
		str_stream << "(" << ((float*)ptr)[0] << "," << ((float*)ptr)[1] << ")";
#endif
	} else {
		std::ostringstream msg;
		msg << "Variable \"" << name() << "\" has unknown type \"" << type()
		    << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		free(ptr);
		return NULL;
	}
	free(ptr);

	str_val = str_stream.str();
	return str_val;
}

void
ArrayVariable::cleanMem()
{
	for (int i = _objects.size() - 1; i >= 0; i--) {
		if (_objects.at(i)->ob_refcnt == 1) {
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

Variables::Variables() {}

Variables::~Variables()
{
	unsigned int i;
	for (auto var : _vars) {
		delete var;
	}
	_vars.clear();
}

void
Variables::registerVariable(const std::string name,
                            const std::string type,
                            const std::string length,
                            const std::string value)
{
	// Look for an already existing variable with the same name
	for (unsigned int i = 0; i < _vars.size(); i++) {
		if (!_vars.at(i)->name().compare(name)) {
			delete _vars.at(i);
			_vars.erase(_vars.begin() + i);
		}
	}

	// Discriminate scalar vs. array
	if (type.find('*') != std::string::npos) {
		registerClMem(name, type, length);
	} else {
		registerScalar(name, type, value);
	}
}

Variable*
Variables::get(unsigned int index)
{
	if (index >= _vars.size()) {
		return NULL;
	}
	return _vars.at(index);
}

Variable*
Variables::get(const std::string name)
{
	for (auto var : _vars) {
		if (!name.compare(var->name())) {
			return var;
		}
	}
	return NULL;
}

size_t
Variables::allocatedMemory()
{
	size_t allocated_mem = 0;
	for (auto var : _vars) {
		if (var->type().find('*') == std::string::npos) {
			continue;
		}
		allocated_mem += var->size();
	}
	return allocated_mem;
}

size_t
Variables::typeToBytes(const std::string type)
{
	unsigned int n = typeToN(type);
	size_t type_size = 0;

	if (type.find("unsigned int") != std::string::npos ||
	    type.find("uivec") != std::string::npos) {
		type_size = sizeof(unsigned int);
	} else if (type.find("int") != std::string::npos ||
	           type.find("ivec") != std::string::npos) {
		type_size = sizeof(int);
	} else if (type.find("float") != std::string::npos ||
	           type.find("vec") != std::string::npos ||
	           type.find("matrix") != std::string::npos) {
		type_size = sizeof(float);
	} else {
		std::ostringstream msg;
		msg << "Unvalid type \"" << type << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		return 0;
	}
	return n * type_size;
}

unsigned int
Variables::typeToN(const std::string type)
{
	unsigned int n = 1;
	if (type.find("vec2") != std::string::npos) {
		n = 2;
	} else if (type.find("vec3") != std::string::npos) {
		n = 3;
	} else if (type.find("vec4") != std::string::npos) {
		n = 4;
	} else if (type.find("vec") != std::string::npos) {
#ifdef HAVE_3D
		n = 4;
#else
		n = 2;
#endif // HAVE_3D
	} else if (type.find("matrix") != std::string::npos) {
#ifdef HAVE_3D
		n = 16;
#else
		n = 4;
#endif // HAVE_3D
	}
	return n;
}

bool
Variables::isSameType(const std::string type_a,
                      const std::string type_b,
                      bool ignore_asterisk)
{
	if (typeToN(type_a) != typeToN(type_b)) {
		return false;
	}

	if (!ignore_asterisk) {
		if ((type_a.find('*') != std::string::npos) &&
		    (type_b.find('*') == std::string::npos)) {
			return false;
		} else if ((type_a.find('*') == std::string::npos) &&
		           (type_b.find('*') != std::string::npos)) {
			return false;
		}
	}

	std::string ta = trimCopy(type_a);
	if (ta.back() == '*') {
		ta.pop_back();
	}
	if ((ta.back() == '2') || (ta.back() == '3') || (ta.back() == '4')) {
		ta.pop_back();
	}
	std::string tb = trimCopy(type_b);
	if (tb.back() == '*') {
		tb.pop_back();
	}
	if ((tb.back() == '2') || (tb.back() == '3') || (tb.back() == '4')) {
		tb.pop_back();
	}

	if (ta.compare(tb)) {
		return false;
	}

	return true;
}

/** @brief Convert the names populated at the tokenizer to variable names
 * @param name The name populated on the tokenizer
 * @return The variable name
 */
std::string
tokNameToVarName(const std::string& name)
{
	std::string var_name = name;
	for (auto suffix : { "_x", "_y", "_z", "_w" }) {
		if (endswith(var_name, suffix)) {
			var_name.erase(var_name.end() - 2, var_name.end());
			break;
		}
	}
	return var_name;
}

std::vector<Variable*>
Variables::exprVariables(const std::string& expr)
{
	std::vector<Variable*> vars;
	auto var_names = tok.exprVariables(expr);
	for (auto var_name : var_names) {
		Variable* var = get(tokNameToVarName(var_name));
		if (std::find(vars.begin(), vars.end(), var) != vars.end())
			continue;
		if (!var) {
			std::ostringstream msg;
			msg << "Variable \"" << var_name
			    << "\", referenced on the expression " << expr
			    << ", cannot be found" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("No such variable");
		}
		vars.push_back(var);
	}
	return vars;
}

void
Variables::solve(const std::string type_name,
                 const std::string value,
                 void* data,
                 const std::string name)
{
	size_t typesize = typeToBytes(type_name);
	if (!typesize) {
		throw std::runtime_error("0 bytes size");
	}
	if (!value.compare("")) {
		LOG(L_ERROR, "Empty value received\n");
		throw std::runtime_error("Empty value string");
	}

	// Ignore whether it is an array or a scalar
	std::string type = trimCopy(type_name);
	if (type.back() == '*')
		type.pop_back();

	if (!type.compare("int")) {
		int val;
		float auxval;
		readComponents(name, value, 1, &auxval);
		val = round(auxval);
		memcpy(data, &val, typesize);
	} else if (!type.compare("unsigned int")) {
		unsigned int val;
		float auxval;
		readComponents(name, value, 1, &auxval);
		val = (unsigned int)round(auxval);
		memcpy(data, &val, typesize);
	} else if (!type.compare("float")) {
		float val;
		readComponents(name, value, 1, &val);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec")) {
		vec val;
#ifdef HAVE_3D
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = auxval[0];
		val.y = auxval[1];
		val.z = auxval[2];
		val.w = auxval[3];
		memcpy(data, &val, typesize);
#else
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = auxval[0];
		val.y = auxval[1];
#endif
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec2")) {
		vec2 val;
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = auxval[0];
		val.y = auxval[1];
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec3")) {
		vec3 val;
		float auxval[3];
		readComponents(name, value, 3, auxval);
		val.x = auxval[0];
		val.y = auxval[1];
		val.z = auxval[2];
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec4")) {
		vec4 val;
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = auxval[0];
		val.y = auxval[1];
		val.z = auxval[2];
		val.w = auxval[3];
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec")) {
		ivec val;
#ifdef HAVE_3D
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = round(auxval[0]);
		val.y = round(auxval[1]);
		val.z = round(auxval[2]);
		val.w = round(auxval[3]);
		memcpy(data, &val, typesize);
#else
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = round(auxval[0]);
		val.y = round(auxval[1]);
#endif
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec2")) {
		ivec2 val;
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = round(auxval[0]);
		val.y = round(auxval[1]);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec3")) {
		ivec3 val;
		float auxval[3];
		readComponents(name, value, 3, auxval);
		val.x = round(auxval[0]);
		val.y = round(auxval[1]);
		val.z = round(auxval[2]);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec4")) {
		ivec4 val;
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = round(auxval[0]);
		val.y = round(auxval[1]);
		val.z = round(auxval[2]);
		val.w = round(auxval[3]);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec")) {
		uivec val;
#ifdef HAVE_3D
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = (unsigned int)round(auxval[0]);
		val.y = (unsigned int)round(auxval[1]);
		val.z = (unsigned int)round(auxval[2]);
		val.w = (unsigned int)round(auxval[3]);
		memcpy(data, &val, typesize);
#else
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = (unsigned int)round(auxval[0]);
		val.y = (unsigned int)round(auxval[1]);
#endif
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec2")) {
		uivec2 val;
		float auxval[2];
		readComponents(name, value, 2, auxval);
		val.x = (unsigned int)round(auxval[0]);
		val.y = (unsigned int)round(auxval[1]);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec3")) {
		uivec3 val;
		float auxval[3];
		readComponents(name, value, 3, auxval);
		val.x = (unsigned int)round(auxval[0]);
		val.y = (unsigned int)round(auxval[1]);
		val.z = (unsigned int)round(auxval[2]);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec4")) {
		uivec4 val;
		float auxval[4];
		readComponents(name, value, 4, auxval);
		val.x = (unsigned int)round(auxval[0]);
		val.y = (unsigned int)round(auxval[1]);
		val.z = (unsigned int)round(auxval[2]);
		val.w = (unsigned int)round(auxval[3]);
		memcpy(data, &val, typesize);
	} else {
		throw std::runtime_error("Invalid variable type");
	}
}

void
Variables::populate(const std::string name)
{
	unsigned int i;
	Variable* var = NULL;
	if (name.compare("")) {
		var = get(name);
		if (!var) {
			std::ostringstream msg;
			msg << "Variable \"" << name << "\" cannot be found" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("No such variable");
		}
		populate(var);
		return;
	}

	for (auto var : _vars) {
		populate(var);
	}
}

void
Variables::populate(Variable* var)
{
	std::ostringstream name;
	const std::string type = trimCopy(var->type());
	if (!type.compare("int")) {
		int val = *(int*)var->get_async();
		tok.registerVariable(var->name(), (float)val);
	} else if (!type.compare("unsigned int")) {
		unsigned int val = *(unsigned int*)var->get_async();
		tok.registerVariable(var->name(), (float)val);
	} else if (!type.compare("float")) {
		float val = *(float*)var->get_async();
		tok.registerVariable(var->name(), (float)val);
	} else if (!type.compare("vec")) {
		vec val = *(vec*)var->get_async();
#ifdef HAVE_3D
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
#else
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
#endif // HAVE_3D
	} else if (!type.compare("vec2")) {
		vec2 val = *(vec2*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
	} else if (!type.compare("vec3")) {
		vec3 val = *(vec3*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
	} else if (!type.compare("vec4")) {
		vec4 val = *(vec4*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
	} else if (!type.compare("ivec")) {
		ivec val = *(ivec*)var->get_async();
#ifdef HAVE_3D
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
#else
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
#endif // HAVE_3D
	} else if (!type.compare("ivec2")) {
		ivec2 val = *(ivec2*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
	} else if (!type.compare("ivec3")) {
		ivec3 val = *(ivec3*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
	} else if (!type.compare("ivec4")) {
		ivec4 val = *(ivec4*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
	} else if (!type.compare("uivec")) {
		uivec val = *(uivec*)var->get_async();
#ifdef HAVE_3D
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
#else
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
#endif // HAVE_3D
	} else if (!type.compare("uivec2")) {
		uivec2 val = *(uivec2*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
	} else if (!type.compare("uivec3")) {
		uivec3 val = *(uivec3*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
	} else if (!type.compare("uivec4")) {
		uivec4 val = *(uivec4*)var->get_async();
		name.str("");
		name << var->name() << "_x";
		tok.registerVariable(name.str(), (float)(val.x));
		name.str("");
		name << var->name() << "_y";
		tok.registerVariable(name.str(), (float)(val.y));
		name.str("");
		name << var->name() << "_z";
		tok.registerVariable(name.str(), (float)(val.z));
		name.str("");
		name << var->name() << "_w";
		tok.registerVariable(name.str(), (float)(val.w));
	} else {
		std::ostringstream msg;
		msg << "\"" << var->name() << "\" declared as \"" << type
		    << "\", which is not a scalar type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint\n");
		LOG0(L_DEBUG, "\tunsigned int\n");
		LOG0(L_DEBUG, "\tfloat\n");
		LOG0(L_DEBUG, "\tvec\n");
		LOG0(L_DEBUG, "\tvec2\n");
		LOG0(L_DEBUG, "\tvec3\n");
		LOG0(L_DEBUG, "\tvec4\n");
		LOG0(L_DEBUG, "\tivec\n");
		LOG0(L_DEBUG, "\tivec2\n");
		LOG0(L_DEBUG, "\tivec3\n");
		LOG0(L_DEBUG, "\tivec4\n");
		LOG0(L_DEBUG, "\tuivec\n");
		LOG0(L_DEBUG, "\tuivec2\n");
		LOG0(L_DEBUG, "\tuivec3\n");
		LOG0(L_DEBUG, "\tuivec4\n");
		throw std::runtime_error("Invalid variable type");
	}
}

void
Variables::registerScalar(const std::string name,
                          const std::string type_name,
                          const std::string value)
{
	std::string type = trimCopy(type_name);
	if (!type.compare("int")) {
		IntVariable* var = new IntVariable(name);
		if (value.compare("")) {
			int val = round(tok.solve(value));
			tok.registerVariable(name, (float)val);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("unsigned int")) {
		UIntVariable* var = new UIntVariable(name);
		if (value.compare("")) {
			unsigned int val = (unsigned int)round(tok.solve(value));
			tok.registerVariable(name, (float)val);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("float")) {
		FloatVariable* var = new FloatVariable(name);
		if (value.compare("")) {
			float val = tok.solve(value);
			tok.registerVariable(name, val);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("vec")) {
		VecVariable* var = new VecVariable(name);
		if (value.compare("")) {
			vec val;
#ifdef HAVE_3D
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = auxval[0];
			val.y = auxval[1];
			val.z = auxval[2];
			val.w = auxval[3];
#else
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = auxval[0];
			val.y = auxval[1];
#endif // HAVE_3D
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("vec2")) {
		Vec2Variable* var = new Vec2Variable(name);
		if (value.compare("")) {
			vec2 val;
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = auxval[0];
			val.y = auxval[1];
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("vec3")) {
		Vec3Variable* var = new Vec3Variable(name);
		if (value.compare("")) {
			vec3 val;
			float auxval[3];
			readComponents(name, value, 3, auxval);
			val.x = auxval[0];
			val.y = auxval[1];
			val.z = auxval[2];
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("vec4")) {
		Vec4Variable* var = new Vec4Variable(name);
		if (value.compare("")) {
			vec4 val;
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = auxval[0];
			val.y = auxval[1];
			val.z = auxval[2];
			val.w = auxval[3];
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("ivec")) {
		IVecVariable* var = new IVecVariable(name);
		if (value.compare("")) {
			ivec val;
#ifdef HAVE_3D
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = round(auxval[0]);
			val.y = round(auxval[1]);
			val.z = round(auxval[2]);
			val.w = round(auxval[3]);
#else
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = round(auxval[0]);
			val.y = round(auxval[1]);
#endif // HAVE_3D
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("ivec2")) {
		IVec2Variable* var = new IVec2Variable(name);
		if (value.compare("")) {
			ivec2 val;
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = round(auxval[0]);
			val.y = round(auxval[1]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("ivec3")) {
		IVec3Variable* var = new IVec3Variable(name);
		if (value.compare("")) {
			ivec3 val;
			float auxval[3];
			readComponents(name, value, 3, auxval);
			val.x = round(auxval[0]);
			val.y = round(auxval[1]);
			val.z = round(auxval[2]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("ivec4")) {
		IVec4Variable* var = new IVec4Variable(name);
		if (value.compare("")) {
			ivec4 val;
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = round(auxval[0]);
			val.y = round(auxval[1]);
			val.z = round(auxval[2]);
			val.w = round(auxval[3]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("uivec")) {
		UIVecVariable* var = new UIVecVariable(name);
		if (value.compare("")) {
			uivec val;
#ifdef HAVE_3D
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = (unsigned int)round(auxval[0]);
			val.y = (unsigned int)round(auxval[1]);
			val.z = (unsigned int)round(auxval[2]);
			val.w = (unsigned int)round(auxval[3]);
#else
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = (unsigned int)round(auxval[0]);
			val.y = (unsigned int)round(auxval[1]);
#endif // HAVE_3D
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("uivec2")) {
		UIVec2Variable* var = new UIVec2Variable(name);
		if (value.compare("")) {
			uivec2 val;
			float auxval[2];
			readComponents(name, value, 2, auxval);
			val.x = (unsigned int)round(auxval[0]);
			val.y = (unsigned int)round(auxval[1]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("uivec3")) {
		UIVec3Variable* var = new UIVec3Variable(name);
		if (value.compare("")) {
			uivec3 val;
			float auxval[3];
			readComponents(name, value, 3, auxval);
			val.x = (unsigned int)round(auxval[0]);
			val.y = (unsigned int)round(auxval[1]);
			val.z = (unsigned int)round(auxval[2]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("uivec4")) {
		UIVec4Variable* var = new UIVec4Variable(name);
		if (value.compare("")) {
			uivec4 val;
			float auxval[4];
			readComponents(name, value, 4, auxval);
			val.x = (unsigned int)round(auxval[0]);
			val.y = (unsigned int)round(auxval[1]);
			val.z = (unsigned int)round(auxval[2]);
			val.w = (unsigned int)round(auxval[3]);
			var->set(&val);
		}
		_vars.push_back(var);
	} else {
		std::ostringstream msg;
		msg << "\"" << name << "\" declared as \"" << type
		    << "\", which is not a valid scalar type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint\n");
		LOG0(L_DEBUG, "\tunsigned int\n");
		LOG0(L_DEBUG, "\tfloat\n");
		LOG0(L_DEBUG, "\tvec\n");
		LOG0(L_DEBUG, "\tvec2\n");
		LOG0(L_DEBUG, "\tvec3\n");
		LOG0(L_DEBUG, "\tvec4\n");
		LOG0(L_DEBUG, "\tivec\n");
		LOG0(L_DEBUG, "\tivec2\n");
		LOG0(L_DEBUG, "\tivec3\n");
		LOG0(L_DEBUG, "\tivec4\n");
		LOG0(L_DEBUG, "\tuivec\n");
		LOG0(L_DEBUG, "\tuivec2\n");
		LOG0(L_DEBUG, "\tuivec3\n");
		LOG0(L_DEBUG, "\tuivec4\n");
		throw std::runtime_error("Invalid scalar variable type");
	}
}

void
Variables::registerClMem(const std::string name,
                         const std::string type_name,
                         const std::string length)
{
	unsigned int n;
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	// Get the type size
	std::string type = trimCopy(type_name);
	type.pop_back(); // Remove the asterisk
	const size_t typesize = typeToBytes(type);
	if (!typesize) {
		std::ostringstream msg;
		msg << "\"" << name << "\" declared as \"" << type
		    << "*\", which is not a valid array type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint*\n");
		LOG0(L_DEBUG, "\tunsigned int*\n");
		LOG0(L_DEBUG, "\tfloat*\n");
		LOG0(L_DEBUG, "\tvec*\n");
		LOG0(L_DEBUG, "\tvec2*\n");
		LOG0(L_DEBUG, "\tvec3*\n");
		LOG0(L_DEBUG, "\tvec4*\n");
		LOG0(L_DEBUG, "\tivec*\n");
		LOG0(L_DEBUG, "\tivec2*\n");
		LOG0(L_DEBUG, "\tivec3*\n");
		LOG0(L_DEBUG, "\tivec4*\n");
		LOG0(L_DEBUG, "\tuivec*\n");
		LOG0(L_DEBUG, "\tuivec2*\n");
		LOG0(L_DEBUG, "\tuivec3*\n");
		LOG0(L_DEBUG, "\tuivec4*\n");
		LOG0(L_DEBUG, "\tmatrix*\n");
		throw std::runtime_error("Invalid array variable type");
	}

	// Get the length
	if (!length.compare("")) {
		std::ostringstream msg;
		msg << "No length specified for variable \"" << name
		    << "\", declared as array (" << type << "*)" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid array length");
	}

	try {
		n = (unsigned int)round(tok.solve(length));
	} catch (...) {
		std::ostringstream msg;
		msg << "Failure evaluating variable \"" << name << "\" length"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid array length");
	}

	if (!n) {
		std::ostringstream msg;
		msg << "Variable \"" << name << "\" has null length" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid array length");
	}

	// Generate the variable
	ArrayVariable* var = new ArrayVariable(name, trimCopy(type_name));

	// Allocate memory on device
	cl_int status;
	cl_mem mem;
	mem = clCreateBuffer(
	    C->context(), CL_MEM_READ_WRITE, n * typesize, NULL, &status);
	if (status != CL_SUCCESS) {
		LOG(L_ERROR, "Allocation failure.\n");
		Aqua::InputOutput::Logger::singleton()->printOpenCLError(status);
		throw std::bad_alloc();
	}
	var->set(&mem);

	_vars.push_back(var);
}

void
Variables::readComponents(const std::string name,
                          const std::string value,
                          unsigned int n,
                          float* v)
{
	float val;
	unsigned int i, j;
	if (n == 0) {
		std::ostringstream msg;
		msg << n << " components required for the variable \"" << name << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("0 components variable registration");
	}
	if (n > 4) {
		LOG(L_ERROR, "No more than 4 components can be required\n");
		std::ostringstream msg;
		msg << n << " components required for the variable \"" << name << "\"."
		    << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("5+ components variable registration");
	}

	// Split and parse subexpressions
	std::vector<std::string> value_strs = split_formulae(value);
	const char* extensions[4] = { "_x", "_y", "_z", "_w" };
	i = 0;
	for (auto s : value_strs) {
		try {
			val = tok.solve(s);
		} catch (...) {
			std::ostringstream msg;
			msg << "parsing variable \"" << name << "\"" << std::endl;
			LOG0(L_DEBUG, msg.str());
			throw std::runtime_error("Variable evaluation error");
		}
		std::ostringstream tok_name;
		if (n != 1) {
			tok_name << name << extensions[i];
		} else {
			tok_name << name;
		}
		tok.registerVariable(tok_name.str(), val);
		v[i] = val;
		if (++i == n) {
			break;
		}
	}
	if (i != n) {
		std::ostringstream msg;
		msg << "Failure parsing the variable \"" << name << "\" value"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << n << " fields expected, " << i << " received" << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid number of fields");
	}
}

}
} // namespace
