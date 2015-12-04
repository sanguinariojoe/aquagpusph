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

#ifndef VARIABLE_H_INCLUDED
#define VARIABLE_H_INCLUDED

#include <Python.h>
#include <CL/cl.h>

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <deque>
#include <sphPrerequisites.h>
#include <Tokenizer/Tokenizer.h>

#ifdef HAVE_3D
    #define VecVariable Vec4Variable
    #define IVecVariable IVec4Variable
    #define UIVecVariable UIVec4Variable
#else
    #define VecVariable Vec2Variable
    #define IVecVariable IVec2Variable
    #define UIVecVariable UIVec2Variable
#endif

namespace Aqua{ namespace InputOutput{

/** @class Variable Variable.h Variable.h
 * @brief A generic variable. Almost useless, use the overloaded classes
 * instead of this one.
 */
class Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    Variable(const char *varname, const char *vartype);

    /** Destructor.
     */
    virtual ~Variable();

    /** Name of the variable
     * @return The name of the variable
     */
    const char* name() const {return (const char*)_name;}

    /** Type of the variable
     * @return The type of the variable
     */
    virtual const char* type() const {return _typename;}

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    virtual size_t typesize() const {return sizeof(void);}

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    virtual size_t size() const {return typesize();}

    /** Get variable pointer basis pointer
     * @return Implementation pointer, NULL for this class.
     */
    virtual void* get(){return NULL;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    virtual void set(void* ptr)=0;

    /** Get a Python interpretation of the variable
     * @param i0 First component to be read, just for array variables.
     * @param n Number of component to be read, just for array variables.
     * @return Python object, NULL for this class.
     */
    virtual PyObject* getPythonObject(int i0=0, int n=0){return NULL;}

    /** Set the variable from a Python object
     * @param obj Python object.
     * @param i0 First component to be set, just for array variables.
     * @param n Number of component to be set, just for array variables.
     * @return true, i.e. an error.
     */
    virtual bool setFromPythonObject(PyObject* obj, int i0=0, int n=0)
    {
        return true;
    }

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    virtual const char* asString(){return NULL;}
private:
    /// Name of the variable
    char* _name;

    /// Type of the variable
    char* _typename;
};

/** @class IntVariable Variable.h Variable.h
 * @brief An integer variable.
 */
class IntVariable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    IntVariable(const char *varname);

    /** Destructor.
     */
    ~IntVariable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(int);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    int* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){_value = *(int*)ptr;}

    /** Get a PyLongObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyLongObject Python object.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyLongObject object.
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    int _value;
};

/** @class UIntVariable Variable.h Variable.h
 * @brief An integer variable.
 */
class UIntVariable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    UIntVariable(const char *varname);

    /** Destructor.
     */
    ~UIntVariable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(unsigned int);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    unsigned int* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){_value = *(unsigned int*)ptr;}

    /** Get a PyLongObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyLongObject Python object.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyLongObject object.
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    unsigned int _value;
};

/** @class FloatVariable Variable.h Variable.h
 * @brief A float variable.
 */
class FloatVariable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    FloatVariable(const char *varname);

    /** Destructor.
     */
    ~FloatVariable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(float);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    float* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){_value = *(float*)ptr;}

    /** Get a PyFloatObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyFloatObject Python object.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyFloatObject object.
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    float _value;
};

/** @class Vec2Variable Variable.h Variable.h
 * @brief A vec2 variable.
 */
class Vec2Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    Vec2Variable(const char *varname);

    /** Destructor.
     */
    ~Vec2Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(vec2);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    vec2* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(vec2));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    vec2 _value;
};

/** @class Vec3Variable Variable.h Variable.h
 * @brief A vec3 variable.
 */
class Vec3Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    Vec3Variable(const char *varname);

    /** Destructor.
     */
    ~Vec3Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(vec3);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    vec3* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(vec3));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    vec3 _value;
};

/** @class Vec4Variable Variable.h Variable.h
 * @brief A vec4 variable.
 */
class Vec4Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    Vec4Variable(const char *varname);

    /** Destructor.
     */
    ~Vec4Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(vec4);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    vec4* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(vec4));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    vec4 _value;
};

/** @class IVec2Variable Variable.h Variable.h
 * @brief A ivec2 variable.
 */
class IVec2Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    IVec2Variable(const char *varname);

    /** Destructor.
     */
    ~IVec2Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(ivec2);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    ivec2* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(ivec2));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    ivec2 _value;
};

/** @class IVec3Variable Variable.h Variable.h
 * @brief A ivec3 variable.
 */
class IVec3Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    IVec3Variable(const char *varname);

    /** Destructor.
     */
    ~IVec3Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(ivec3);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    ivec3* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(ivec3));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    ivec3 _value;
};

/** @class IVec4Variable Variable.h Variable.h
 * @brief A ivec4 variable.
 */
class IVec4Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    IVec4Variable(const char *varname);

    /** Destructor.
     */
    ~IVec4Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(ivec4);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    ivec4* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(ivec4));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    ivec4 _value;
};

/** @class UIVec2Variable Variable.h Variable.h
 * @brief A uivec2 variable.
 */
class UIVec2Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    UIVec2Variable(const char *varname);

    /** Destructor.
     */
    ~UIVec2Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(uivec2);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    uivec2* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(uivec2));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    uivec2 _value;
};

/** @class IVec3Variable Variable.h Variable.h
 * @brief A uivec3 variable.
 */
class UIVec3Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    UIVec3Variable(const char *varname);

    /** Destructor.
     */
    ~UIVec3Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(uivec3);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    uivec3* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(uivec3));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    uivec3 _value;
};

/** @class UIVec4Variable Variable.h Variable.h
 * @brief A uivec4 variable.
 */
class UIVec4Variable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     */
    UIVec4Variable(const char *varname);

    /** Destructor.
     */
    ~UIVec4Variable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(uivec4);}

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    uivec4* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){memcpy(&_value, ptr, sizeof(uivec4));}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();
private:
    /// Variable value
    uivec4 _value;
};

/** @class FloatVariable Variable.h Variable.h
 * @brief A float variable.
 */
class ArrayVariable : public Variable
{
public:
    /** Constructor.
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    ArrayVariable(const char *varname, const char *vartype);

    /** Destructor.
     */
    ~ArrayVariable();

    /** Get the cl_mem type size.
     * @return cl_mem type size (in bytes)
     * @note In order to know the typesize of the components into the array you
     * may use Variables::typeToBytes()
     */
    size_t typesize() const {return sizeof(cl_mem);}

    /** Get the array size.
     * @return Array allocated memory (in bytes)
     * @note In order to get the length of the array the command
     * size() / typesize() can be used
     */
    size_t size() const;

    /** Get variable pointer basis pointer
     * @return Implementation pointer.
     */
    cl_mem* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy.
     */
    void set(void* ptr){_value = *(cl_mem*)ptr;}

    /** Get a PyArrayObject interpretation of the variable
     * @param i0 First component to be read.
     * @param n Number of component to be read, 0 to read all available memory,
     * i.e. All the array after i0.
     * @return PyArrayObject Python object. NULL if the memory cannot be read.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** Set the variable from a Python object
     * @param obj PyArrayObject object.
     * @param i0 0, otherwise an error will be returned.
     * @param n 0, otherwise an error will be returned.
     * @return false if all gone right, true otherwise.
     */
    bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const char* asString();

    /** Get a component text representation
     * @param i Index of the component to be extracted.
     * @return The component represented as a string, NULL in case of errors.
     */
    const char* asString(size_t i);
private:
    /// Check for abandoned python objects to destroy them.
    void cleanMem();

    /// Variable value
    cl_mem _value;
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
    std::deque<void*> _data;
    /** @brief List of helpers data array storages for the Python objects
     * @see getPythonObject()
     * @see _data
     */
    std::deque<PyObject*> _objects;
};

// ---------------------------------------------------------------------------
// Variables manager
// ---------------------------------------------------------------------------

/** @class VariableManager VariableManager.h VariableManager.h
 * @brief Variables manager, which can interpret the types on the fly.
 */
class Variables
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
     * @return false if all gone right, true otherwise
     */
    bool registerVariable(const char* name,
                          const char* type,
                          const char* length,
                          const char* value);

    /** Get a variable.
     * @param index Index of the variable.
     * @return Variable, NULL if the variable cannot be found.
     */
    Variable* get(unsigned int index);

    /** Get a variable.
     * @param name Name of the variable.
     * @return Variable, NULL if the variable cannot be found.
     */
    Variable* get(const char* name);

    /** Get all the registered variables.
     * @return Variable, NULL if the variable cannot be found.
     */
    std::deque<Variable*> getAll(){return _vars;}

    /** Get the number of variables.
     * @return Number of registered variables.
     */
    unsigned int size() const {return _vars.size();}

    /** Get the allocated memory.
     * @return Allocated memory on device. Just the arrays can contribute to this value.
     */
    size_t allocatedMemory();

    /** Convert a type name to bytes.
     * @param type Type name.
     * @return Type size in bytes, 0 if the type is not recognized.
     */
    size_t typeToBytes(const char* type) const;

    /** Get the number of components of a type name.
     * @param type Type name.
     * @return Number of components (1 for not recognized types).
     */
    unsigned int typeToN(const char* type) const;

    /** Get if two types strings are the same one.
     * @param type_a First type name.
     * @param type_b Second type name.
     * @param ignore_asterisk true to ignore the asterisk of arrays.
     * @return true if it is the same type, false otherwise.
     */
    bool isSameType(const char* type_a,
                    const char* type_b,
                    bool ignore_asterisk=true);

    /** Solve a string, interpreting the variables.
     * @param type_name Type of the output desired value.
     * @param value Expression to evaluate.
     * @param data Allocated memory where the result should be stored.
     * @param name Variable name to register in the tokenizer.
     * @return false if all gone right, true otherwise.
     * @note typeToBytes(type) bytes should be allocated in data.
     */
    bool solve(const char *type_name,
               const char *value,
               void *data,
               const char* name="NULL");

    /** @brief Populate variables in order that the tokenizer may get the
     * updated value.
     * @param name Name of the variable to be populated, NULL if all the
     * variables should be populated.
     * @return false if all gone right, true otherwise
     */
    bool populate(const char* name=NULL);

    /** @brief Populate a variable in order that the tokenizer may know the
     * updated value.
     * @param var Variable to be populated.
     * @return false if all gone right, true otherwise
     */
    bool populate(Variable* var);
private:

    /** Register a scalar variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param value Variable value, NULL for arrays. It is optional for
     * scalar variables.
     * @return false if all gone right, true otherwise
     */
    bool registerScalar(const char* name,
                        const char* type,
                        const char* value);
    /** Register a cl_mem variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param length Array length, 1 for scalars, 0 for arrays that will
     * not be allocated at the start (for instance the heads of chains,
     * which requires the number of cells).
     * @return false if all gone right, true otherwise
     */
    bool registerClMem(const char* name,
                       const char* type,
                       const char* length);

    /** Read a set of components from a value array.
     * @param name Name of the variable. It is used to register variables in
     * the tokenizer and to report errors.
     * @param value Value string where the components are extracted from.
     * @param n Number of components to read.
     * @param v Allocated array where the components should be stored
     * @return false if all gone right, true otherwise.
     */
    bool readComponents(const char* name,
                        const char* value,
                        unsigned int n,
                        float* v);


    /// Set of available variables
    std::deque<Variable*> _vars;
    /// Tokenizer to evaluate variables
    Tokenizer tok;
};

}}  // namespace

#endif // VARIABLE_H_INCLUDED
