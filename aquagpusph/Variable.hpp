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

#ifndef VARIABLE_H_INCLUDED
#define VARIABLE_H_INCLUDED

#include <Python.h>
#include <string>
#include <vector>

/** @def PY_ARRAY_UNIQUE_SYMBOL
 * @brief Define the extension module which this Python stuff should be linked
 * to.
 *
 * In AQUAgpusph all the Python stuff is linked in the same group
 * AQUA_ARRAY_API
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

#include "sphPrerequisites.hpp"
#include "Tokenizer/Tokenizer.hpp"

#ifndef MAX_TYPE_NAME_LEN
#define MAX_TYPE_NAME_LEN 128
#endif

namespace Aqua {
namespace InputOutput {

/** @class Variable Variable.h Variable.h
 * @brief A generic variable. Almost useless, use the overloaded classes
 * instead of this one.
 */
class DECLDIR Variable
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 * @param vartype Type of the variable.
	 */
	Variable(const std::string varname, const std::string vartype);

	/** @brief Destructor.
	 */
	virtual ~Variable();

	/** @brief Let efficiently know whether the variable is an array or not
	 *
	 * @return true if the variable is an array, false if the variable is a
	 * scalar
	 */
	virtual bool isArray() = 0;

	/** @brief Let efficiently know whether the variable is a scalar or not
	 *
	 * @return true if the variable is a scalar, false if the variable is an
	 * array
	 */
	inline bool isScalar() { return !isArray(); }

	/** @brief Name of the variable
	 * @return The name of the variable
	 */
	std::string name() const { return _name; }

	/** @brief Type of the variable
	 * @return The type of the variable
	 */
	virtual std::string type() const { return _typename; }

	/** @brief Get the variable type size.
	 * @return Variable type size (in bytes)
	 */
	virtual size_t typesize() const { return 0; }

	/** @brief Get the variable type size.
	 * @return Variable type size (in bytes)
	 */
	virtual size_t size() const { return typesize(); }

	/** @brief Get variable pointer basis pointer
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 * @return Implementation pointer, NULL for this class.
	 */
	virtual inline void* get(bool UNUSED_PARAM synced = true) { return NULL; }

	/** @brief Get variable pointer basis pointer
	 *
	 * Sometimes a tool would make events micromanagement, becoming useful to
	 * can retrieve the value in an asynchronous way. That is for instance the
	 * case of events manipulated in callbacks
	 *
	 * @return Implementation pointer
	 */
	inline void* get_async() { return get(false); }

	/** @brief Set variable from memory
	 *
	 * If @p synced, this method is setting a completed user event as the
	 * variable writing event. Remeber to call this at the end of any
	 * overloaded function
	 * @param ptr Unused parameter, it is of use for overloaded functions
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 */
	virtual void set(void* ptr, bool synced = true);

	/** @brief Set variable from memory
	 *
	 * Sometimes a tool would make events micromanagement, becoming useful to
	 * can retrieve the value in an asynchronous way. That is for instance the
	 * case of events manipulated in callbacks
	 *
	 * @param ptr Memory to copy
	 */
	inline void set_async(void* ptr) { set(ptr, false); }

	/** @brief Get a Python interpretation of the variable
	 * @param i0 First component to be read, just for array variables.
	 * @param n Number of component to be read, just for array variables.
	 * @return Python object, NULL for this class.
	 */
	virtual PyObject* getPythonObject(int UNUSED_PARAM i0 = 0,
	                                  int UNUSED_PARAM n = 0)
	{
		return NULL;
	}

	/** @brief Set the variable from a Python object
	 * @param obj Python object.
	 * @param i0 First component to be set, just for array variables.
	 * @param n Number of component to be set, just for array variables.
	 * @return true, i.e. an error.
	 */
	virtual bool setFromPythonObject(PyObject UNUSED_PARAM *obj,
	                                 int UNUSED_PARAM i0 = 0,
	                                 int UNUSED_PARAM n = 0)
	{
		return true;
	}

	/** @brief Get the variable text representation
	 * @return The variable represented as a string, NULL in case of errors.
	 */
	virtual const std::string asString(bool UNUSED_PARAM synced = true)
	{
		return "";
	}

	/**
	 * \defgroup VariableEventsGroup Events tracked on variables
	 * @brief On AQUAgpusph the variables are used to track the events, and
	 * thus when each tool can be executed.
	 *
	 * The rationale is that each variable can store 2 types of events,
	 * reading events and writing events.
	 * Along this line, reading processes shall only wait for the last writing
	 * event, while writing processes shall wait for both the last writing
	 * event as well all the recorded reading events
	 * @{
	 */

	/** @brief Set the variable current writing event
	 *
	 * clRetainEvent() is called on the provided event, while clReleaseEvent()
	 * is called on the eventually previous stored event
	 *
	 * @remarks Events are used even for non-OpenCL variables
	 * @param event OpenCL event
	 */
	void setEvent(cl_event event);

	/** @brief Alias of InputOutput::Variable::setEvent()
	 */
	inline void setWritingEvent(cl_event event) { setEvent(event); }

	/** @brief Returns the last writing event associated to this variable
	 *
	 * @return OpenCL event
	 */
	inline cl_event getEvent() const { return _event; }

	/** @brief Alias of InputOutput::Variable::getEvent()
	 */
	inline cl_event getWritingEvent() const { return getEvent(); }

	/** @brief Add a new reading event to the variable
	 *
	 * clRetainEvent() is called on the provided event.
	 *
	 * On top of that, the list of reading events is traversed, calling
	 * clReleaseEvent() on those completed and droping them from the list
	 *
	 * @remarks Events are used even for non-OpenCL variables
	 * @param event OpenCL event
	 */
	void addReadingEvent(cl_event event);

	/** @brief Get the list of reading events
	 *
	 * The list of reading events is only valid until
	 * InputOutput::Variable::addReadingEvent() is called again
	 *
	 * @remarks Events are used even for non-OpenCL variables
	 * @param event OpenCL event
	 */
	inline std::vector<cl_event> getReadingEvents() const
	{
		return _reading_events;
	}

	/** @ingroup VariableEventsGroup
	 * @brief Wait for variable reading and writing events to be completed
	 *
	 * This function is tracking the syncing state, so clWaitForEvents() is not
	 * called until the previous syncing has been completed.
	 *
	 * This is obviously a blocking function
	 * @param readonly true if it is just needed to sync for reading
	 * operations, false otherwise
	 */
	void sync(bool readonly = false);

	/**
	 * @}
	 */
  protected:
	/** @ingroup VariableEventsGroup
	 * @brief Clean up the list of reading events
	 *
	 * All the events marked as completed are released (calling
	 * clReleaseEvent()) and dropped from the list
	 */
	void cleanReadingEvents();

	/** @brief Create a competed user event to be set as the writing event when
	 * synced writing operations are carried out.
	 * @return The completed event
	 * @throw std::runtime_error If the OpenCL event creation fails
	 */
	static cl_event createDummyEvent();

  private:
	/// Name of the variable
	std::string _name;

	/// Type of the variable
	std::string _typename;

	/// Last writing event
	cl_event _event;

	/// List of reading events
	std::vector<cl_event> _reading_events;

	/** @brief Shortcut to avoid calling the expensive OpenCL API
	 *
	 * This is true when all the events (reading and writing) are marked as
	 * completed
	 */
	bool _synced;

	/** @brief Shortcut to avoid calling the expensive OpenCL API
	 *
	 * This is true when all the writing events are marked as completed
	 */
	bool _synced_for_read;
};

/** @class ScalarVariable Variable.h Variable.h
 * @brief A generic Scalar variable.
 */
template<class T>
class DECLDIR ScalarVariable : public Variable
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 * @param vartype Type of the variable.
	 */
	ScalarVariable(const std::string varname, const std::string vartype)
	  : Variable(varname, vartype)
	{
	}

	/** @brief Destructor.
	 */
	~ScalarVariable() {};

	/** @brief Report that the varaible is not an array
	 *
	 * @return false
	 */
	inline bool isArray() { return false; }

	/** @brief Get the variable type size.
	 * @return Variable type size (in bytes)
	 */
	size_t typesize() const { return sizeof(T); }

	/** @brief Get variable pointer basis pointer
	 *
	 * This is a blocking operation, that will retain the program until the
	 * underlying variable event is complete. For asynchronous variable
	 * retrieval see get_async()
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 * @return Implementation pointer
	 */
	inline void* get(bool synced = true)
	{
		if (synced)
			sync(true);
		return &_value;
	}

	/** @brief Get variable value
	 *
	 * This is a blocking operation, that will retain the program until the
	 * underlying variable event is complete. For asynchronous variable
	 * retrieval see get_async()
	 * @param value Output value
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 */
	inline void get(T& value, bool synced = true)
	{
		if (synced)
			sync(true);
		value = _value;
	}

	/** @brief Set variable from memory
	 *
	 * This is a blocking operation, that will retain the program until the
	 * underlying variable event is complete.
	 * @param ptr Memory to copy
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 */
	inline void set(void* ptr, bool synced = true)
	{
		if (synced)
			sync();
		memcpy(&_value, ptr, sizeof(T));
		this->Variable::set(ptr, synced);
	}

	/** @brief Set variable value
	 *
	 * This is a blocking operation, that will retain the program until the
	 * underlying variable event is complete.
	 * @param value New value
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 */
	inline void set(T& value, bool synced = true)
	{
		if (synced)
			sync();
		_value = value;
		this->Variable::set(NULL, synced);
	}

	/** @brief Get the variable value
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 * @return Output value
	 */
	inline T value(bool synced = true)
	{
		if (synced)
			sync(true);
		return _value;
	}

	/** @brief Set the variable value
	 * @param synced true if the operation shall block the execution until the
	 * variable events are dispatched
	 * @return Output value
	 */
	inline void value(T& value, bool synced = true)
	{
		this->set(value, synced);
	}

	/** @brief Get a Python Object representation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter.
	 * @return Python object.
	 */
	virtual PyObject* getPythonObject(int i0 = 0, int n = 0) = 0;

	/** @brief Set the variable from a Python object
	 * @param obj Python object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0) = 0;

	/** @brief Get the variable text representation
	 * @param synced true if the function shall block until the last writing
	 * event is dispatched, false otherwise
	 * @return The variable represented as a string, NULL in case of errors.
	 */
	virtual const std::string asString(bool synced) = 0;

  protected:
	/// Variable value
	T _value;
};

/** @class ScalarVariable Variable.h Variable.h
 * @brief A generic Scalar variable.
 */
template<class T>
class DECLDIR ScalarNumberVariable : public ScalarVariable<T>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 * @param vartype Type of the variable.
	 */
	ScalarNumberVariable(const std::string varname, const std::string vartype)
	  : ScalarVariable<T>(varname, vartype)
	{
	}

	/** @brief Destructor.
	 */
	~ScalarNumberVariable(){};

	/** @brief Get the variable text representation
	 * @param synced true if the function shall block until the last writing
	 * event is dispatched, false otherwise
	 * @return The variable represented as a string, NULL in case of errors.
	 */
	virtual const std::string asString(bool synced = true);
};

/** @class IntVariable Variable.h Variable.h
 * @brief An integer variable.
 */
class DECLDIR IntVariable final : public ScalarNumberVariable<icl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	IntVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~IntVariable(){};

	/** @brief Get a PyLongObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter.
	 * @return PyLongObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyLongObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class LongVariable Variable.h Variable.h
 * @brief A 64bits integer variable.
 */
class DECLDIR LongVariable final : public ScalarNumberVariable<lcl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	LongVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~LongVariable(){};

	/** @brief Get a PyLongObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter.
	 * @return PyLongObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyLongObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class UIntVariable Variable.h Variable.h
 * @brief An unsigned integer variable.
 */
class DECLDIR UIntVariable final : public ScalarNumberVariable<uicl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	UIntVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~UIntVariable(){};

	/** @brief Get a PyLongObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return PyLongObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyLongObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class ULongVariable Variable.h Variable.h
 * @brief An integer variable.
 */
class DECLDIR ULongVariable final : public ScalarNumberVariable<ulcl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	ULongVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~ULongVariable(){};

	/** @brief Get a PyLongObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return PyLongObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyLongObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class FloatVariable Variable.h Variable.h
 * @brief A float variable.
 */
class DECLDIR FloatVariable final : public ScalarNumberVariable<fcl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	FloatVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~FloatVariable(){};

	/** @brief Get a PyFloatObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return PyFloatObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyFloatObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class DoubleVariable Variable.h Variable.h
 * @brief A double variable.
 */
class DECLDIR DoubleVariable final : public ScalarNumberVariable<dcl>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 */
	DoubleVariable(const std::string varname);

	/** @brief Destructor.
	 */
	~DoubleVariable(){};

	/** @brief Get a PyFloatObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return PyFloatObject Python object.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** @brief Set the variable from a Python object
	 * @param obj PyFloatObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);
};

/** @class ScalarVecVariable Variable.h Variable.h
 * @brief A generic Scalar variable, of 2 or more components.
 */
template<class T>
class DECLDIR ScalarVecVariable : public ScalarVariable<T>
{
  public:
	/** @brief Constructor.
	 * @param varname Name of the variable.
	 * @param vartype Type of the variable.
	 * @param dims Number of components of the type.
	 * @param np_type Type of the numpy associated object.
	 */
	ScalarVecVariable(const std::string varname,
	                  const std::string vartype,
	                  const unsigned int dims,
	                  int np_type);

	/** @brief Destructor.
	 */
	~ScalarVecVariable(){};

	/** Get a PyArrayObject interpretation of the variable
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return PyArrayObject Python object (PyArray_FLOAT subtype).
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** Set the variable from a Python object
	 * @param obj PyArrayObject object (PyArray_FLOAT subtype).
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);

	/** @brief Get the variable text representation
	 * @param synced true if the function shall block until the last writing
	 * event is dispatched, false otherwise
	 * @return The variable represented as a string, NULL in case of errors.
	 */
	virtual const std::string asString(bool synced = true);

  protected:
	/** @brief Check that a Python object is compatible with the variable type
	 *
	 * This method is checking that the Python object is an array with the same
	 * number of components of the variable type. The internal data type is not
	 * checked.
	 *
	 * @param obj Python object object.
	 * @return false if the Python object matchs with the variable type, true
	 * otherwise.
	 */
	bool checkPyhonObjectDims(PyObject* obj);

  private:
	/// Number of components of the type
	unsigned int _dims;
	/// Number of components of the type
	int _np_type;
};

#define __DECLARE_AQUA_VEC(NAME, TYPE)                                         \
	class DECLDIR NAME final : public ScalarVecVariable<TYPE>                  \
	{                                                                          \
	public:                                                                    \
		NAME(const std::string varname);                                       \
		~NAME(){};                                                             \
	};

/** @class IVec2Variable Variable.h Variable.h
 * @brief A ivec2 variable.
 */
__DECLARE_AQUA_VEC(IVec2Variable, ivec2)

/** @class IVec3Variable Variable.h Variable.h
 * @brief A ivec3 variable.
 */
__DECLARE_AQUA_VEC(IVec3Variable, ivec3)

/** @class IVec4Variable Variable.h Variable.h
 * @brief A ivec4 variable.
 */
__DECLARE_AQUA_VEC(IVec4Variable, ivec4)

/** @class IVec8Variable Variable.h Variable.h
 * @brief A ivec8 variable.
 */
__DECLARE_AQUA_VEC(IVec8Variable, ivec8)

/** @class LVec2Variable Variable.h Variable.h
 * @brief A lvec2 variable.
 */
__DECLARE_AQUA_VEC(LVec2Variable, lvec2)

/** @class LVec3Variable Variable.h Variable.h
 * @brief A lvec3 variable.
 */
__DECLARE_AQUA_VEC(LVec3Variable, lvec3)

/** @class LVec4Variable Variable.h Variable.h
 * @brief A lvec4 variable.
 */
__DECLARE_AQUA_VEC(LVec4Variable, lvec4)

/** @class LVec8Variable Variable.h Variable.h
 * @brief A ivec8 variable.
 */
__DECLARE_AQUA_VEC(LVec8Variable, lvec8)

/** @class UIVec2Variable Variable.h Variable.h
 * @brief A uivec2 variable.
 */
__DECLARE_AQUA_VEC(UIVec2Variable, uivec2)

/** @class UIVec3Variable Variable.h Variable.h
 * @brief A uivec3 variable.
 */
__DECLARE_AQUA_VEC(UIVec3Variable, uivec3)

/** @class UIVec4Variable Variable.h Variable.h
 * @brief A uivec4 variable.
 */
__DECLARE_AQUA_VEC(UIVec4Variable, uivec4)

/** @class UIVec8Variable Variable.h Variable.h
 * @brief A uivec8 variable.
 */
__DECLARE_AQUA_VEC(UIVec8Variable, uivec8)

/** @class ULVec2Variable Variable.h Variable.h
 * @brief A ulvec2 variable.
 */
__DECLARE_AQUA_VEC(ULVec2Variable, ulvec2)

/** @class ULVec3Variable Variable.h Variable.h
 * @brief A ulvec3 variable.
 */
__DECLARE_AQUA_VEC(ULVec3Variable, ulvec3)

/** @class ULVec4Variable Variable.h Variable.h
 * @brief A ulvec4 variable.
 */
__DECLARE_AQUA_VEC(ULVec4Variable, ulvec4)

/** @class ULVec8Variable Variable.h Variable.h
 * @brief A ulvec8 variable.
 */
__DECLARE_AQUA_VEC(ULVec8Variable, ulvec8)

/** @class Vec2Variable Variable.h Variable.h
 * @brief A vec2 variable.
 */
__DECLARE_AQUA_VEC(Vec2Variable, vec2)

/** @class Vec3Variable Variable.h Variable.h
 * @brief A vec3 variable.
 */
__DECLARE_AQUA_VEC(Vec3Variable, vec3)

/** @class Vec4Variable Variable.h Variable.h
 * @brief A vec4 variable.
 */
__DECLARE_AQUA_VEC(Vec4Variable, vec4)

/** @class Vec8Variable Variable.h Variable.h
 * @brief A vec8 variable.
 */
__DECLARE_AQUA_VEC(Vec8Variable, vec8)

/** @class DVec2Variable Variable.h Variable.h
 * @brief A dvec2 variable.
 */
__DECLARE_AQUA_VEC(DVec2Variable, dvec2)

/** @class DVec3Variable Variable.h Variable.h
 * @brief A dvec3 variable.
 */
__DECLARE_AQUA_VEC(DVec3Variable, dvec3)

/** @class DVec4Variable Variable.h Variable.h
 * @brief A dvec4 variable.
 */
__DECLARE_AQUA_VEC(DVec4Variable, dvec4)

/** @class DVec8Variable Variable.h Variable.h
 * @brief A dvec8 variable.
 */
__DECLARE_AQUA_VEC(DVec8Variable, dvec8)

/** @name Variable aliases
 * @brief Some alias to make easier to create preprocessor directives
 */
/// @{
typedef IntVariable iclVariable;
typedef LongVariable lclVariable;
typedef UIntVariable uiclVariable;
typedef ULongVariable ulclVariable;
typedef FloatVariable fclVariable;
typedef DoubleVariable dclVariable;
typedef IVec2Variable ivec2Variable;
typedef IVec3Variable ivec3Variable;
typedef IVec4Variable ivec4Variable;
typedef IVec8Variable ivec8Variable;
typedef LVec2Variable lvec2Variable;
typedef LVec3Variable lvec3Variable;
typedef LVec4Variable lvec4Variable;
typedef LVec8Variable lvec8Variable;
typedef UIVec2Variable uivec2Variable;
typedef UIVec3Variable uivec3Variable;
typedef UIVec4Variable uivec4Variable;
typedef UIVec8Variable uivec8Variable;
typedef ULVec2Variable ulvec2Variable;
typedef ULVec3Variable ulvec3Variable;
typedef ULVec4Variable ulvec4Variable;
typedef ULVec8Variable ulvec8Variable;
typedef Vec2Variable vec2Variable;
typedef Vec3Variable vec3Variable;
typedef Vec4Variable vec4Variable;
typedef Vec8Variable vec8Variable;
typedef DVec2Variable dvec2Variable;
typedef DVec3Variable dvec3Variable;
typedef DVec4Variable dvec4Variable;
typedef DVec8Variable dvec8Variable;
/// @}

/** @class ArrayVariable Variable.h Variable.h
 * @brief An array variable.
 */
class DECLDIR ArrayVariable : public Variable
{
  public:
	/** Constructor.
	 * @param varname Name of the variable.
	 * @param vartype Type of the variable.
	 */
	ArrayVariable(const std::string varname, const std::string vartype);

	/** Destructor.
	 */
	~ArrayVariable();

	/** @brief Report that the varaible is an array
	 *
	 * @return true
	 */
	bool isArray();

	/** Get the cl_mem type size.
	 * @return cl_mem type size (in bytes)
	 * @note In order to know the typesize of the components into the array you
	 * may use Variables::typeToBytes()
	 */
	size_t typesize() const { return sizeof(cl_mem); }

	/** Get the array size.
	 * @return Array allocated memory (in bytes)
	 * @note In order to get the length of the array the command
	 * size() / Variables::typeToBytes(type()) can be used
	 */
	size_t size() const;

	/** Get variable pointer basis pointer
	 * @param synced Unused parameter
	 * @return Implementation pointer.
	 */
	void* get(bool UNUSED_PARAM synced = false) { return &_value; }

	/** Set variable from memory
	 * @param ptr Memory to copy.
	 * @param synced Unused parameter
	 * @throw std::runtime_error if the variable is not reallocatable and it
	 * has been already set
	 */
	void set(void* ptr, bool synced = false);

	/** @brief Get if a variable is reallocatable
	 * @return true if the variable is marked as reallocatable, false otherwise
	 */
	inline bool reallocatable() const { return _reallocatable; }

	/** @brief Set if a variable is reallocatable
	 * @param is true if the variable is reallocatable, false otherwise
	 * @note Non-reallocatable variables are in general way more efficient
	 */
	inline void reallocatable(bool is) { _reallocatable = is; }

	/** Get a PyArrayObject interpretation of the variable
	 * @param i0 First component to be read.
	 * @param n Number of component to be read, 0 to read all available memory,
	 * i.e. All the array after i0.
	 * @return PyArrayObject Python object. NULL if the memory cannot be read.
	 */
	PyObject* getPythonObject(int i0 = 0, int n = 0);

	/** Set the variable from a Python object
	 * @param obj PyArrayObject object.
	 * @param i0 ignored parameter.
	 * @param n ignored parameter
	 * @return false if all gone right, true otherwise.
	 */
	bool setFromPythonObject(PyObject* obj, int i0 = 0, int n = 0);

	/** Get the variable text representation
	 * @param synced true if the function shall block until the last writing
	 * event is dispatched, false otherwise
	 * @return The variable represented as a string, NULL in case of errors.
	 */
	const std::string asString(bool synced = true);

	/** Get a component text representation
	 * @param i Index of the component to be extracted.
	 * @param synced true if the function shall block until the last writing
	 * event is dispatched, false otherwise
	 * @return The component represented as a string, NULL in case of errors.
	 */
	const std::string asString(size_t i, bool synced = true);

  private:
	/// Check for abandoned python objects to destroy them.
	void cleanMem();

	/// Variable value
	cl_mem _value;

	/// Mark if this memory buffer can be reallocated or not (false by default)
	bool _reallocatable;

	/** @brief List of helpers data array storages for the Python objects
	 *
	 * The memory array inside numpy objects must be dynamically allocated and
	 * preserved, otherwise wrong values will be received in the Python script.
	 *
	 * On the other hand, this memory is not automatically freed by Python when
	 * the object is destroyed, and therefore we need to call Py_INCREF after
	 * the object generation to assert that Python is not automatically
	 * detroying it, such that we can control the reference count, deleting the
	 * memory array and the Python object when 0 is reached.
	 * @see getPythonObject()
	 * @see _objects
	 */
	std::vector<void*> _data;

	/** @brief List of helpers data array storages for the Python objects
	 * @see getPythonObject()
	 * @see _data
	 */
	std::vector<PyObject*> _objects;
};

// ---------------------------------------------------------------------------
// Variables manager
// ---------------------------------------------------------------------------

/** @class Variables Variables.h Variables.h
 * @brief Variables manager, which can interpret the types on the fly.
 */
class DECLDIR Variables
{
  public:
	/** Constructor.
	 */
	Variables();

	/** Destructor.
	 */
	~Variables();

	/** Register a new variable.
	 * @param name Name of the variable.
	 * @param type Type of the variable.
	 * @param length Array length, 1 for scalars, 0 for arrays that will
	 * not be allocated at the start (for instance the heads of chains,
	 * which requires the number of cells).
	 * @param value Variable value, NULL for arrays. It is optional for
	 * scalar variables.
	 */
	void registerVariable(const std::string name,
	                      const std::string type,
	                      const std::string length,
	                      const std::string value);

	/** Get a variable.
	 * @param index Index of the variable.
	 * @return Variable, NULL if the variable cannot be found.
	 */
	Variable* get(unsigned int index);

	/** Get a variable.
	 * @param name Name of the variable.
	 * @return Variable, NULL if the variable cannot be found.
	 */
	Variable* get(const std::string name);

	/** Get all the registered variables.
	 * @return Variable, NULL if the variable cannot be found.
	 */
	std::vector<Variable*> getAll() const { return _vars; }

	/** Get the number of variables.
	 * @return Number of registered variables.
	 */
	unsigned int size() const { return _vars.size(); }

	/** Get the allocated memory.
	 * @return Allocated memory on device. Just the arrays can contribute to
	 * this value.
	 */
	size_t allocatedMemory();

	/** Convert a type name to bytes.
	 * @param type Type name.
	 * @return Type size in bytes, 0 if the type is not recognized.
	 */
	static size_t typeToBytes(const std::string type);

	/** Get the number of components of a type name.
	 * @param type Type name.
	 * @return Number of components (1 for not recognized types).
	 */
	static unsigned int typeToN(const std::string type);

	/** Get if two types strings are the same one.
	 * @param type_a First type name.
	 * @param type_b Second type name.
	 * @param ignore_asterisk true to ignore the asterisk of arrays.
	 * @return true if it is the same type, false otherwise.
	 */
	static bool isSameType(const std::string type_a,
	                       const std::string type_b,
	                       bool ignore_asterisk = true);

	/** Get the list of variables called on a expression.
	 * @param expr Expression to parse.
	 * @return The list of input variables.
	 * @see Aqua::Tokenizer::exprVariables()
	 */
	std::vector<Variable*> exprVariables(const std::string& expr);

	/** Solve a string, interpreting the variables.
	 * @param type_name Type of the output desired value.
	 * @param value Expression to evaluate.
	 * @param data Allocated memory where the result should be stored.
	 * @param name Variable name to register in the tokenizer.
	 * @note typeToBytes(type) bytes should be allocated in data.
	 */
	void solve(const std::string type_name,
	           const std::string value,
	           void* data,
	           const std::string name = "NULL");

	/** @brief Populate variables in order that the tokenizer may get the
	 * updated value.
	 * @param name Name of the variable to be populated, "" if all the
	 * variables should be populated.
	 */
	void populate(const std::string name = "");

	/** @brief Populate a variable in order that the tokenizer may know the
	 * updated value.
	 * @param var Variable to be populated.
	 */
	void populate(Variable* var);

	/** @brief Resolve some type aliases
	 *
	 * Several types have aliases to make them easier to use, which might
	 * depend on the device arch:
	 *
	 *  - "int32" = "int"
	 *  - "int64" = "long"
	 *  - "uint32" = "unsigned int"
	 *  - "uint64" = "unsigned long"
	 *  - "size_t" = Either "unsigned int" or "unsigned long"
	 *  - "ssize_t" = Either "int" or "long"
	 *  - "svec" = Either "uivec" or "ulvec"
	 *  - "svec2" = Either "uivec2" or "ulvec2"
	 *  - "svec3" = Either "uivec3" or "ulvec3"
	 *  - "svec4" = Either "uivec4" or "ulvec4"
	 *  - "svec8" = Either "uivec8" or "ulvec8"
	 *  - "ssvec" = Either "ivec" or "lvec"
	 *  - "ssvec2" = Either "ivec2" or "lvec2"
	 *  - "ssvec3" = Either "ivec3" or "lvec3"
	 *  - "ssvec4" = Either "ivec4" or "lvec4"
	 *  - "ssvec8" = Either "ivec8" or "lvec8"
	 * @param t Input aliased type
	 * @return Base type
	 * @see Aqua::CalcServer::CalcServer::device_addr_bits()
	 */
	static std::string typeAlias(const std::string& t);

  private:
	/** Register a scalar variable
	 * @param name Name of the variable.
	 * @param type Type of the variable.
	 * @param value Variable value, NULL for arrays. It is optional for
	 * scalar variables.
	 */
	void registerScalar(const std::string name,
	                    const std::string type,
	                    const std::string value);
	/** Register a cl_mem variable
	 * @param name Name of the variable.
	 * @param type Type of the variable.
	 * @param length Array length, 1 for scalars, 0 for arrays that will
	 * not be allocated at the start (for instance the heads of chains,
	 * which requires the number of cells).
	 */
	void registerClMem(const std::string name,
	                   const std::string type,
	                   const std::string length);

	/** Read a set of components from a value array.
	 * @param name Name of the variable. It is used to register variables in
	 * the tokenizer and to report errors.
	 * @param value Value string where the components are extracted from.
	 * @param n Number of components to read.
	 * @param v Allocated array where the components should be stored
	 */
	template<typename T>
	T solve(const std::string& name, const std::string& value);

	/// Set of available variables
	std::vector<Variable*> _vars;
	/// Tokenizer to evaluate variables
	Tokenizer tok;
};

}
} // namespace

#endif // VARIABLE_H_INCLUDED
