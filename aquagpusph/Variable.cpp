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

#include <algorithm>
#include <sstream>
#include "Variable.hpp"
#include "AuxiliarMethods.hpp"
#include "InputOutput/Logger.hpp"
#include "CalcServer/CalcServer.hpp"

namespace Aqua {
namespace InputOutput {

static std::ostringstream pyerr;

Variable::Variable(const std::string varname, const std::string vartype)
  : _name(varname)
  , _typename(vartype)
  , _event(NULL)
  , _synced(true)
  , _synced_for_read(true)
{
	_event = Variable::createDummyEvent();
}

void
Variable::set(void* ptr, bool synced)
{
	if(!synced)
		return;

	cl_int err_code;
	cl_event e = Variable::createDummyEvent();
	setWritingEvent(e);
	err_code = clReleaseEvent(e);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing the dummy writing event for \"" << name()
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

cl_event Variable::createDummyEvent()
{
	cl_int err_code;
	auto C = CalcServer::CalcServer::singleton();
	cl_event event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure creating an user event" << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clSetUserEventStatus(event, CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting the CL_COMPLETE status" << std::endl;
		LOG(L_ERROR, msg.str());
		Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	return event;
}

template<class T>
const std::string
ScalarNumberVariable<T>::asString(bool synced)
{
	std::ostringstream msg;
	msg << this->value(synced);
	return msg.str();
}

IntVariable::IntVariable(const std::string varname)
  : ScalarNumberVariable<icl>(varname, "int")
{
	_value = 0;
}

PyObject*
IntVariable::getPythonObject(int i0, int n)
{
	long val = value();
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

	icl value = (icl)PyLong_AsLong(obj);
	set(&value);

	return false;
}

LongVariable::LongVariable(const std::string varname)
  : ScalarNumberVariable<lcl>(varname, "long")
{
	_value = 0;
}

PyObject*
LongVariable::getPythonObject(int i0, int n)
{
	long val = value();
	return PyLong_FromLong(val);
}

bool
LongVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyLong_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyLongObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	lcl value = (lcl)PyLong_AsLong(obj);
	this->value(value);

	return false;
}

UIntVariable::UIntVariable(const std::string varname)
  : ScalarNumberVariable<uicl>(varname, "unsigned int")
{
	_value = 0u;
}

PyObject*
UIntVariable::getPythonObject(int i0, int n)
{
	unsigned long val = value();
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

	uicl value = (uicl)PyLong_AsUnsignedLong(obj);
	this->value(value);

	return false;
}

ULongVariable::ULongVariable(const std::string varname)
  : ScalarNumberVariable<ulcl>(varname, "unsigned long")
{
	_value = 0u;
}

PyObject*
ULongVariable::getPythonObject(int i0, int n)
{
	unsigned long val = value();
	return PyLong_FromUnsignedLong(val);
}

bool
ULongVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyLong_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyLongObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	ulcl value = (ulcl)PyLong_AsUnsignedLong(obj);
	set(&value);

	return false;
}

FloatVariable::FloatVariable(const std::string varname)
  : ScalarNumberVariable<fcl>(varname, "float")
{
	_value = 0.f;
}

PyObject*
FloatVariable::getPythonObject(int i0, int n)
{
	double val = value();
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

	fcl value = (fcl)PyFloat_AsDouble(obj);
	set(&value);

	return false;
}

DoubleVariable::DoubleVariable(const std::string varname)
  : ScalarNumberVariable<dcl>(varname, "double")
{
	_value = 0.0;
}

PyObject*
DoubleVariable::getPythonObject(int i0, int n)
{
	double val = value();
	return PyFloat_FromDouble(val);
}

bool
DoubleVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (!PyFloat_Check(obj)) {
		pyerr.str("");
		pyerr << "Variable \"" << name() << "\" expected a PyFloatObject"
		      << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}

	fcl value = (fcl)PyFloat_AsDouble(obj);
	set(&value);

	return false;
}

template<class T>
ScalarVecVariable<T>::ScalarVecVariable(const std::string varname,
                                        const std::string vartype,
                                        const unsigned int dims,
                                        int np_type)
  : ScalarVariable<T>(varname, vartype)
  , _dims(dims)
  , _np_type(np_type)
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
PyObject*
ScalarVecVariable<T>::getPythonObject(int i0, int n)
{
	T v = this->value();
	npy_intp dims[] = { _dims };
	return PyArray_SimpleNewFromData(1, dims, _np_type, v.s);
}

std::string npy_type_name(int np_type)
{
	std::string name("NPY_NOTYPE");
	switch(np_type) {
		case NPY_BOOL:
			name = "NPY_BOOL";
			break;
		case NPY_INT8:
			name = "NPY_INT8";
			break;
		case NPY_INT16:
			name = "NPY_INT16";
			break;
		case NPY_INT32:
			name = "NPY_INT32";
			break;
		case NPY_INT64:
			name = "NPY_INT64";
			break;
		case NPY_UINT8:
			name = "NPY_UINT8";
			break;
		case NPY_UINT16:
			name = "NPY_UINT16";
			break;
		case NPY_UINT32:
			name = "NPY_UINT32";
			break;
		case NPY_UINT64:
			name = "NPY_UINT64";
			break;
		case NPY_FLOAT16:
			name = "NPY_FLOAT16";
			break;
		case NPY_FLOAT32:
			name = "NPY_FLOAT32";
			break;
		case NPY_FLOAT64:
			name = "NPY_FLOAT64";
			break;
		default:
			name = "NPY_NOTYPE";
	}
	return name;
}

template<class T>
bool
ScalarVecVariable<T>::setFromPythonObject(PyObject* obj, int i0, int n)
{
	if (checkPyhonObjectDims(obj))
		return true;
	PyArrayObject* array_obj = (PyArrayObject*)obj;
	if (PyArray_TYPE(array_obj) != _np_type) {
		pyerr.str("");
		pyerr << "Variable \"" << this->name() << "\" expected a "
		      << npy_type_name(_np_type) << " but got a "
		      << npy_type_name(PyArray_TYPE(array_obj)) << std::endl;
		PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
		return true;
	}
	this->set(PyArray_DATA(array_obj));
	return false;
}

template<class T>
const std::string
ScalarVecVariable<T>::asString(bool synced)
{
	T val = this->value(synced);
	std::ostringstream msg;
	msg << "(";
	for (unsigned int i = 0; i < _dims; i++) {
		msg << val.s[i] << ",";
	}
	std::string str = msg.str();
	str.back() = ')';
	return str;
}

#define __DEFINE_AQUA_VEC(NAME, TYPE, DIMS, NPTYPE)                            \
	NAME::NAME(const std::string varname)                                      \
		: ScalarVecVariable(varname, #TYPE, DIMS, NPTYPE)                      \
	{                                                                          \
		for (unsigned int i = 0; i < DIMS; i++)                                \
			_value.s[i] = 0;                                                   \
	}

__DEFINE_AQUA_VEC(IVec2Variable, ivec2, 2, NPY_INT32)
__DEFINE_AQUA_VEC(IVec3Variable, ivec3, 3, NPY_INT32)
__DEFINE_AQUA_VEC(IVec4Variable, ivec4, 4, NPY_INT32)
__DEFINE_AQUA_VEC(IVec8Variable, ivec8, 8, NPY_INT32)

__DEFINE_AQUA_VEC(LVec2Variable, lvec2, 2, NPY_INT64)
__DEFINE_AQUA_VEC(LVec3Variable, lvec3, 3, NPY_INT64)
__DEFINE_AQUA_VEC(LVec4Variable, lvec4, 4, NPY_INT64)
__DEFINE_AQUA_VEC(LVec8Variable, lvec8, 8, NPY_INT64)

__DEFINE_AQUA_VEC(UIVec2Variable, uivec2, 2, NPY_UINT32)
__DEFINE_AQUA_VEC(UIVec3Variable, uivec3, 3, NPY_UINT32)
__DEFINE_AQUA_VEC(UIVec4Variable, uivec4, 4, NPY_UINT32)
__DEFINE_AQUA_VEC(UIVec8Variable, uivec8, 8, NPY_UINT32)

__DEFINE_AQUA_VEC(ULVec2Variable, ulvec2, 2, NPY_UINT64)
__DEFINE_AQUA_VEC(ULVec3Variable, ulvec3, 3, NPY_UINT64)
__DEFINE_AQUA_VEC(ULVec4Variable, ulvec4, 4, NPY_UINT64)
__DEFINE_AQUA_VEC(ULVec8Variable, ulvec8, 8, NPY_UINT64)

__DEFINE_AQUA_VEC(Vec2Variable, vec2, 2, NPY_FLOAT32)
__DEFINE_AQUA_VEC(Vec3Variable, vec3, 3, NPY_FLOAT32)
__DEFINE_AQUA_VEC(Vec4Variable, vec4, 4, NPY_FLOAT32)
__DEFINE_AQUA_VEC(Vec8Variable, vec8, 8, NPY_FLOAT32)

__DEFINE_AQUA_VEC(DVec2Variable, dvec2, 2, NPY_FLOAT64)
__DEFINE_AQUA_VEC(DVec3Variable, dvec3, 3, NPY_FLOAT64)
__DEFINE_AQUA_VEC(DVec4Variable, dvec4, 4, NPY_FLOAT64)
__DEFINE_AQUA_VEC(DVec8Variable, dvec8, 8, NPY_FLOAT64)

ArrayVariable::ArrayVariable(const std::string varname,
                             const std::string vartype)
  : Variable(varname, vartype)
  , _value(NULL)
  , _reallocatable(false)
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

void
ArrayVariable::set(void* ptr, bool synced)
{
	if (_value && !_reallocatable) {
		std::ostringstream msg;
		msg << "Array variable \"" << name()
		    << "\", cannot be set, because it has not been marked as "
		    << "reallocatable" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("No reallocatable variable");
	}
	_value = *(cl_mem*)ptr;
	this->Variable::set(ptr, synced);
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
	int pytype = NPY_NOTYPE;
	if (startswith(type(), "unsigned int") ||
	    startswith(type(), "uivec")) {
		pytype = NPY_UINT32;
	} else if (startswith(type(), "unsigned long") ||
	           startswith(type(), "ulvec")) {
		pytype = NPY_UINT64;
	} else if (startswith(type(), "int") ||
	           startswith(type(), "ivec")) {
		pytype = NPY_INT32;
	} else if (startswith(type(), "long") ||
	           startswith(type(), "lvec")) {
		pytype = NPY_INT64;
	} else if (startswith(type(), "float") ||
	           startswith(type(), "vec")) {
		pytype = NPY_FLOAT32;
	} else if (startswith(type(), "double") ||
	           startswith(type(), "dvec")) {
		pytype = NPY_FLOAT64;
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
ArrayVariable::asString(bool synced)
{
	cl_mem* val = (cl_mem*)get(synced);
	std::ostringstream msg;
	msg << val;
	return msg.str();
}

/** @brief Stringify a value pointer
 * @param ptr Pointer to the value
 * @param n Number of fields into the data
 */
template<typename T>
const std::string
valptr_as_string(T* ptr, unsigned int n)
{
	if (!n)
		return "";
	std::ostringstream out;
	if (n > 1)
		out << "(";
	for (unsigned int i = 0; i < n - 1; n++)
		out << ptr[i] << ",";
	out << ptr[n - 1];
	if (n > 1)
		out << ")";
	return out.str();
}

const std::string
ArrayVariable::asString(size_t i, bool synced)
{
	auto C = CalcServer::CalcServer::singleton();
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
	cl_event event_wait = getWritingEvent();
	cl_int err_code = clEnqueueReadBuffer(
	    C->command_queue(CalcServer::CalcServer::cmd_queue::cmd_queue_new),
	    _value,
	    CL_TRUE,
	    i * type_size,
	    type_size,
	    ptr,
	    synced ? 1 : 0,
	    synced ? &event_wait : NULL,
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

	const unsigned int n = Variables::typeToN(type());
	std::string res = "";
	if (startswith(type(), "int") ||
	    startswith(type(), "ivec")) {
		res = valptr_as_string((icl*)ptr, n);
	} else if (startswith(type(), "long") ||
	    startswith(type(), "lvec")) {
		res = valptr_as_string((lcl*)ptr, n);
	} else if (startswith(type(), "unsigned int") ||
	    startswith(type(), "uivec")) {
		res = valptr_as_string((uicl*)ptr, n);
	} else if (startswith(type(), "unsigned long") ||
	    startswith(type(), "ulvec")) {
		res = valptr_as_string((ulcl*)ptr, n);
	} else if (startswith(type(), "float") ||
	    startswith(type(), "vec")) {
		res = valptr_as_string((fcl*)ptr, n);
	} else if (startswith(type(), "double") ||
	    startswith(type(), "dvec")) {
		res = valptr_as_string((dcl*)ptr, n);
	} else {
		std::ostringstream msg;
		msg << "Variable \"" << name() << "\" has unknown type \"" << type()
		    << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		free(ptr);
		return NULL;
	}
	free(ptr);

	return res;
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
		type_size = sizeof(uicl);
	} else if (type.find("unigned long") != std::string::npos ||
	           type.find("ulvec") != std::string::npos) {
		type_size = sizeof(ulcl);
	} else if (type.find("int") != std::string::npos ||
	           type.find("ivec") != std::string::npos) {
		type_size = sizeof(icl);
	} else if (type.find("long") != std::string::npos ||
	           type.find("lvec") != std::string::npos) {
		type_size = sizeof(lcl);
	} else if (type.find("double") != std::string::npos ||
	           type.find("dvec") != std::string::npos) {
		type_size = sizeof(dcl);
	} else if (type.find("float") != std::string::npos ||
	           type.find("vec") != std::string::npos ||
	           type.find("matrix") != std::string::npos) {
		type_size = sizeof(fcl);
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
	} else if (type.find("vec8") != std::string::npos) {
		n = 8;
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
	ta = typeAlias(ta);
	std::string tb = trimCopy(type_b);
	if (tb.back() == '*') {
		tb.pop_back();
	}
	tb = typeAlias(tb);

	if (typeToN(ta) != typeToN(tb)) {
		return false;
	}

	if (std::isdigit(ta.back()))
		ta.pop_back();
	if (std::isdigit(tb.back()))
		tb.pop_back();

	return ta == tb;
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

#define SPLIT_EXPR_AND_CHECK(var_name, expr, exprs_out, n)                     \
	std::vector<std::string> exprs_out = split_formulae(expr);                 \
	if (exprs_out.size() < n) {                                                \
		std::ostringstream msg;                                                \
		msg << "Failure parsing the expression \"" << expr                     \
		    << "\" for variable \"" << var_name << "\"" << std::endl;          \
		LOG(L_ERROR, msg.str());                                               \
		msg.str("");                                                           \
		msg << n << " fields required instead of " << exprs_out.size()         \
		    << std::endl;                                                      \
		LOG0(L_DEBUG, msg.str());                                              \
		throw std::runtime_error("Invalid number of fields");                  \
	}

#define SCALAR_SOLVER(vartype, toktype)                                        \
	template<>                                                                 \
	vartype                                                                    \
	Variables::solve(const std::string& name, const std::string& value)        \
	{                                                                          \
		SPLIT_EXPR_AND_CHECK(name, value, value_strs, 1);                      \
		toktype val = tok.solve<toktype>(value_strs[0]);                       \
		tok.registerVariable<toktype>(name, val);                              \
		return (vartype)val;                                                   \
	}

SCALAR_SOLVER(icl, int32_t)
SCALAR_SOLVER(lcl, int64_t)
SCALAR_SOLVER(uicl, uint32_t)
SCALAR_SOLVER(ulcl, uint64_t)
SCALAR_SOLVER(fcl, float)
SCALAR_SOLVER(dcl, double)

#define VEC_SOLVER(vartype, toktype, n)                                        \
	template<>                                                                 \
	vartype                                                                    \
	Variables::solve(const std::string& name, const std::string& value)        \
	{                                                                          \
		vartype out;                                                           \
		const char* extensions[16] = { "_x",  "_y",  "_z",  "_w",              \
		                               "_yx", "_yy", "_yz", "_yw",             \
		                               "_zx", "_zy", "_zz", "_zw",             \
		                               "_wx", "_wy", "_wz", "_ww" };           \
		SPLIT_EXPR_AND_CHECK(name, value, value_strs, n);                      \
		for (unsigned int i = 0; i < n; i++) {                                 \
			toktype val = tok.solve<toktype>(value_strs[i]);                   \
			tok.registerVariable(name + extensions[i], val);                   \
			out.s[i] = val;                                                    \
		}                                                                      \
		return out;                                                            \
	}

// NOTE: 3 components OpenCL vectors are actually an alias for the 4 components
VEC_SOLVER(ivec2, int, 2)
// VEC_SOLVER(ivec3, int, 3)
VEC_SOLVER(ivec4, int, 4)
VEC_SOLVER(ivec8, int, 8)
VEC_SOLVER(lvec2, long, 2)
// VEC_SOLVER(lvec3, long, 3)
VEC_SOLVER(lvec4, long, 4)
VEC_SOLVER(lvec8, long, 8)
VEC_SOLVER(uivec2, unsigned int, 2)
// VEC_SOLVER(uivec3, unsigned int, 3)
VEC_SOLVER(uivec4, unsigned int, 4)
VEC_SOLVER(uivec8, unsigned int, 8)
VEC_SOLVER(ulvec2, unsigned long, 2)
// VEC_SOLVER(ulvec3, unsigned long, 3)
VEC_SOLVER(ulvec4, unsigned long, 4)
VEC_SOLVER(ulvec8, unsigned long, 8)
VEC_SOLVER(vec2, float, 2)
// VEC_SOLVER(vec3, float, 3)
VEC_SOLVER(vec4, float, 4)
VEC_SOLVER(vec8, float, 8)
VEC_SOLVER(dvec2, double, 2)
// VEC_SOLVER(dvec3, double, 3)
VEC_SOLVER(dvec4, double, 4)
VEC_SOLVER(dvec8, double, 8)

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
	type = typeAlias(type);

	if (!type.compare("int")) {
		icl val = solve<icl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("long")) {
		lcl val = solve<lcl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("unsigned int")) {
		uicl val = solve<uicl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("unsigned long")) {
		ulcl val = solve<ulcl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("float")) {
		fcl val = solve<fcl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("double")) {
		dcl val = solve<dcl>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec")) {
		vec val = solve<vec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec2")) {
		vec2 val = solve<vec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec3")) {
		vec3 val = solve<vec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec4")) {
		vec4 val = solve<vec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("vec8")) {
		vec8 val = solve<vec8>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("dvec")) {
		dvec val = solve<dvec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("dvec2")) {
		dvec2 val = solve<dvec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("dvec3")) {
		dvec3 val = solve<dvec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("dvec4")) {
		dvec4 val = solve<dvec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("dvec8")) {
		dvec8 val = solve<dvec8>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec")) {
		ivec val = solve<ivec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec2")) {
		ivec2 val = solve<ivec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec3")) {
		ivec3 val = solve<ivec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec4")) {
		ivec4 val = solve<ivec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ivec8")) {
		ivec8 val = solve<ivec8>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("lvec")) {
		lvec val = solve<lvec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("lvec2")) {
		lvec2 val = solve<lvec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("lvec3")) {
		lvec3 val = solve<lvec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("lvec4")) {
		lvec4 val = solve<lvec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("lvec8")) {
		lvec8 val = solve<lvec8>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec")) {
		uivec val = solve<uivec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec2")) {
		uivec2 val = solve<uivec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec3")) {
		uivec3 val = solve<uivec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec4")) {
		uivec4 val = solve<uivec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("uivec8")) {
		uivec8 val = solve<uivec8>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ulvec")) {
		ulvec val = solve<ulvec>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ulvec2")) {
		ulvec2 val = solve<ulvec2>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ulvec3")) {
		ulvec3 val = solve<ulvec3>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ulvec4")) {
		ulvec4 val = solve<ulvec4>(name, value);
		memcpy(data, &val, typesize);
	} else if (!type.compare("ulvec8")) {
		ulvec8 val = solve<ulvec8>(name, value);
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
	const char* extensions[16] = { "_x",  "_y",  "_z",  "_w",
	                               "_yx", "_yy", "_yz", "_yw",
	                               "_zx", "_zy", "_zz", "_zw",
	                               "_wx", "_wy", "_wz", "_ww" };
	const std::string type = trimCopy(var->type());
	if (!type.compare("int")) {
		icl val = *(icl*)var->get_async();
		tok.registerVariable<int>(var->name(), (int)val);
	} else if (!type.compare("long")) {
		lcl val = *(lcl*)var->get_async();
		tok.registerVariable<int64_t>(var->name(), (int64_t)val);
	} else if (!type.compare("unsigned int")) {
		uicl val = *(uicl*)var->get_async();
		tok.registerVariable<unsigned int>(var->name(), (unsigned int)val);
	} else if (!type.compare("unsigned long")) {
		ulcl val = *(ulcl*)var->get_async();
		tok.registerVariable<uint64_t>(var->name(), (uint64_t)val);
	} else if (!type.compare("float")) {
		fcl val = *(fcl*)var->get_async();
		tok.registerVariable<float>(var->name(), (float)val);
	} else if (!type.compare("double")) {
		dcl val = *(dcl*)var->get_async();
		tok.registerVariable<double>(var->name(), (double)val);
	} else if (startswith(type, "ivec")) {
		const unsigned int n = typeToN(type);
		icl* ptr = (icl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<int32_t>(var->name() + extensions[i],
			                              (int32_t)(ptr[i]));
		}
	} else if (startswith(type, "lvec")) {
		const unsigned int n = typeToN(type);
		lcl* ptr = (lcl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<int64_t>(var->name() + extensions[i],
			                              (int64_t)(ptr[i]));
		}
	} else if (startswith(type, "uivec")) {
		const unsigned int n = typeToN(type);
		uicl* ptr = (uicl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<uint32_t>(var->name() + extensions[i],
			                               (uint32_t)(ptr[i]));
		}
	} else if (startswith(type, "ulvec")) {
		const unsigned int n = typeToN(type);
		ulcl* ptr = (ulcl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<uint64_t>(var->name() + extensions[i],
			                               (uint64_t)(ptr[i]));
		}
	} else if (startswith(type, "dvec")) {
		const unsigned int n = typeToN(type);
		dcl* ptr = (dcl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<double>(var->name() + extensions[i],
			                             (double)(ptr[i]));
		}
	} else if (startswith(type, "vec")) {
		const unsigned int n = typeToN(type);
		fcl* ptr = (fcl*)(var->get_async());
		for (unsigned int i = 0; i < n; i++) {
			tok.registerVariable<float>(var->name() + extensions[i],
			                            (float)(ptr[i]));
		}
	} else {
		std::ostringstream msg;
		msg << "\"" << var->name() << "\" declared as \"" << type
		    << "\", which is not a valid scalar type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint\n");
		LOG0(L_DEBUG, "\tlong\n");
		LOG0(L_DEBUG, "\tunsigned int\n");
		LOG0(L_DEBUG, "\tunsigned long\n");
		LOG0(L_DEBUG, "\tfloat\n");
		LOG0(L_DEBUG, "\tdouble\n");
		for (auto prefix : {"vec", "dvec", "ivec", "lvec", "uivec", "ulvec"}) {
			for (auto n : {2, 3, 4, 8}) {
				LOG0(L_DEBUG, std::string("\t") + prefix + std::to_string(n) +
				              "\n");
			}
		}
		throw std::runtime_error("Invalid variable type");
	}
}

std::string
Variables::typeAlias(const std::string& t)
{
	if (t == "int32")
		return "int";
	else if (t == "int64")
		return "long";
	else if ((t == "uint32") || (t == "uint"))
		return "unsigned int";
	else if (t == "uint64")
		return "unsigned long";
	else if (t == "size_t") {
		auto C = CalcServer::CalcServer::singleton();
		if (C->device_addr_bits() == 64)
			return "unsigned long";
		return "unsigned int";
	} else if (t == "ssize_t") {
		auto C = CalcServer::CalcServer::singleton();
		if (C->device_addr_bits() == 64)
			return "long";
		return "int";
	} else if (startswith(t, "svec")) {
		auto C = CalcServer::CalcServer::singleton();
		if (C->device_addr_bits() == 64)
			return replaceAllCopy(t, "svec", "ulvec");
		return replaceAllCopy(t, "svec", "uivec");
	} else if (startswith(t, "ssvec")) {
		auto C = CalcServer::CalcServer::singleton();
		if (C->device_addr_bits() == 64)
			return replaceAllCopy(t, "ssvec", "lvec");
		return replaceAllCopy(t, "ssvec", "ivec");
	}
	// Vector types base names, i.e. int2, int3, int4, float2, float3, ...
	else if (startswith(t, "int") && (t != "int")) {
		return replaceAllCopy(t, "int", "ivec");
	} else if (startswith(t, "long") && (t != "long")) {
		return replaceAllCopy(t, "long", "lvec");
	} else if (startswith(t, "uint") && (t != "uint")) {
		return replaceAllCopy(t, "uint", "uivec");
	} else if (startswith(t, "ulong") && (t != "ulong")) {
		return replaceAllCopy(t, "ulong", "ulvec");
	} else if (startswith(t, "float") && (t != "float")) {
		return replaceAllCopy(t, "float", "vec");
	} else if (startswith(t, "double") && (t != "double")) {
		return replaceAllCopy(t, "double", "dvec");
	}
	return t;
}

#define REG_VEC_BLOCK(TYPE)                                                    \
	else if (!type.compare(#TYPE)) {                                           \
		TYPE ## Variable* var = new TYPE ## Variable(name);                    \
		if (value.compare("")) {                                               \
			TYPE val = solve<TYPE>(name, value);                               \
			var->set(&val);                                                    \
		}                                                                      \
		_vars.push_back(var);                                                  \
	}

void
Variables::registerScalar(const std::string name,
                          const std::string type_name,
                          const std::string value)
{
	LOG(L_INFO,
	    std::string("Registering the scalar \"") + name + "\" (" +
	    type_name + ") = \"" + value + "\"\n");
	std::string type = typeAlias(trimCopy(type_name));
	if (!type.compare("int")) {
		IntVariable* var = new IntVariable(name);
		if (value.compare("")) {
			icl val = solve<icl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("long")) {
		LongVariable* var = new LongVariable(name);
		if (value.compare("")) {
			lcl val = solve<lcl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("unsigned int")) {
		UIntVariable* var = new UIntVariable(name);
		if (value.compare("")) {
			uicl val = solve<uicl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("unsigned long")) {
		ULongVariable* var = new ULongVariable(name);
		if (value.compare("")) {
			ulcl val = solve<ulcl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("float")) {
		FloatVariable* var = new FloatVariable(name);
		if (value.compare("")) {
			fcl val = solve<fcl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	} else if (!type.compare("double")) {
		DoubleVariable* var = new DoubleVariable(name);
		if (value.compare("")) {
			dcl val = solve<dcl>(name, value);
			var->set(&val);
		}
		_vars.push_back(var);
	}
	REG_VEC_BLOCK(vec)
	REG_VEC_BLOCK(vec2)
	REG_VEC_BLOCK(vec3)
	REG_VEC_BLOCK(vec4)
	REG_VEC_BLOCK(vec8)
	REG_VEC_BLOCK(dvec)
	REG_VEC_BLOCK(dvec2)
	REG_VEC_BLOCK(dvec3)
	REG_VEC_BLOCK(dvec4)
	REG_VEC_BLOCK(dvec8)
	REG_VEC_BLOCK(ivec)
	REG_VEC_BLOCK(ivec2)
	REG_VEC_BLOCK(ivec3)
	REG_VEC_BLOCK(ivec4)
	REG_VEC_BLOCK(ivec8)
	REG_VEC_BLOCK(lvec)
	REG_VEC_BLOCK(lvec2)
	REG_VEC_BLOCK(lvec3)
	REG_VEC_BLOCK(lvec4)
	REG_VEC_BLOCK(lvec8)
	REG_VEC_BLOCK(uivec)
	REG_VEC_BLOCK(uivec2)
	REG_VEC_BLOCK(uivec3)
	REG_VEC_BLOCK(uivec4)
	REG_VEC_BLOCK(uivec8)
	REG_VEC_BLOCK(ulvec)
	REG_VEC_BLOCK(ulvec2)
	REG_VEC_BLOCK(ulvec3)
	REG_VEC_BLOCK(ulvec4)
	REG_VEC_BLOCK(ulvec8)
	else {
		std::ostringstream msg;
		msg << "\"" << name << "\" declared as \"" << type
		    << "\", which is not a valid scalar type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint\n");
		LOG0(L_DEBUG, "\tlong\n");
		LOG0(L_DEBUG, "\tunsigned int\n");
		LOG0(L_DEBUG, "\tunsigned long\n");
		for (auto alias : {"int32", "int64", "uint32", "uint64", "size_t"}) {
			LOG0(L_DEBUG, std::string("\t") + alias + " (same than " +
			              typeAlias(alias) + ")\n");
		}
		LOG0(L_DEBUG, "\tfloat\n");
		LOG0(L_DEBUG, "\tdouble\n");
		for (auto prefix : {"vec", "dvec", "ivec", "lvec", "uivec", "ulvec"}) {
			for (auto n : {2, 3, 4, 8}) {
				LOG0(L_DEBUG, std::string("\t") + prefix + std::to_string(n) +
				              "\n");
			}
		}
		throw std::runtime_error("Invalid scalar variable type");
	}
}

void
Variables::registerClMem(const std::string name,
                         const std::string type_name,
                         const std::string length)
{
	LOG(L_INFO,
	    std::string("Registering the array \"") + name + "\" (" +
	    type_name + ") of length \"" + length + "\"\n");

	size_t n;
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();
	// Get the type size
	std::string type = trimCopy(type_name);
	type.pop_back(); // Remove the asterisk
	type = typeAlias(type);
	const size_t typesize = typeToBytes(type);
	if (!typesize) {
		std::ostringstream msg;
		msg << "\"" << name << "\" declared as \"" << type
		    << "*\", which is not a valid array type" << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "Valid types are:\n");
		LOG0(L_DEBUG, "\tint*\n");
		LOG0(L_DEBUG, "\tlong*\n");
		LOG0(L_DEBUG, "\tunsigned int*\n");
		LOG0(L_DEBUG, "\tunsigned long*\n");
		for (auto alias : {"int32", "int64", "uint32", "uint64", "size_t"}) {
			LOG0(L_DEBUG, std::string("\t") + alias + "* (same than " +
			              typeAlias(alias) + "*)\n");
		}
		LOG0(L_DEBUG, "\tfloat*\n");
		LOG0(L_DEBUG, "\tdouble*\n");
		for (auto prefix : {"vec", "dvec", "ivec", "lvec", "uivec", "ulvec"}) {
			for (auto n : {2, 3, 4, 8}) {
				LOG0(L_DEBUG, std::string("\t") + prefix + "*" +
				              std::to_string(n) + "\n");
			}
		}
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
		n = round(tok.solve<size_t>(length));
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
	ArrayVariable* var = new ArrayVariable(name, type + '*');

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

}
} // namespace
