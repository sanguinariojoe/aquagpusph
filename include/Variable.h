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

#include <CL/cl.h>

#include <string>
#include <vector>
#include <sphPrerequisites.h>
#include <PyCast.h>
#include <Tokenizer/Tokenizer.h>

namespace Aqua{ namespace InputOutput{

/** @class Variable Variable.h Variable.h
 * @brief A generic variable
 * @note In general, this class assumes the behaviour of a scalar variable, in
 * such a way that you can call almost everysingle method without caring to cast
 * the variable to its actual type before. A remarkable exception is set()
 * method, which is purely virtual for safe purposes.
 */
class Variable
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    Variable(const std::string varname, const std::string vartype);

    /** @brief Destructor
     */
    virtual ~Variable() {};

    /** @brief Let efficiently know whether the variable is an array or not
     *
     * @return true if the variable is an array, false if the variable is a
     * scalar. true in this base class.
     */
    virtual const bool isArray() const {return true;}

    /** @brief Let efficiently know whether the variable is a scalar or not
     *
     * @return true if the variable is a scalar, false if the variable is an
     * array
     */
    inline const bool isScalar() const {return !isArray();}

    /** @brief Name of the variable
     * @return The name of the variable
     */
    const std::string& name() const {return _name;}

    /** @brief Type of the variable
     * @return The type of the variable
     */
    const std::string& type() const {return _typename;}

    /** @brief Get the variable type size.
     * @return Variable inner type size (in bytes). 0 for this base class.
     */
    virtual const size_t typesize() const {return 0;}

    /** @brief Get the variable size
     *
     * For scalar variables this matchs typesize(). Conversely, for array
     * variables this returns the length of the array multiplied by typesize()
     *
     * @return Variable outer size (in bytes). typesize() in this case
     */
    virtual const size_t size() const {return typesize();}

    /** @brief Get variable pointer basis pointer
     * @return Implementation pointer, NULL for this class.
     */
    virtual const void* get(){return NULL;}

    /** @brief Set variable from memory
     * @param ptr Memory to copy.
     */
    virtual void set(const void* ptr)=0;

    /** @brief Get a Python interpretation of the variable
     * @param i0 First component to be read, just for array variables.
     * @param n Number of component to be read, just for array variables.
     * @return Python object, NULL for this class.
     */
    virtual PyObject* getPythonObject(int i0=0, int n=0){return NULL;}

    /** @brief Set the variable from a Python object
     * @param obj Python object.
     * @param i0 First component to be set, just for array variables.
     * @param n Number of component to be set, just for array variables.
     * @return true, which means that it finished with errors.
     */
    virtual const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0)
    {
        return true;
    }

    /** @brief Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    virtual const std::string asString() {return "";}

    /** @brief Set the variable current event
     *
     * clRetainEvent() is called on top of the provided event, while any
     * previously stored event is conveniently released calling clReleaseEvent()
     *
     * Might exists other events linked to this variable, but since those are
     * considered predecessors, we can just forgive about them, tracking just
     * the last one
     *
     * @remarks Events are used even for non-OpenCL variables
     * @param event OpenCL event
     */
    void setEvent(cl_event event);    

    /** @brief Returns the last event associated to this variable
     * @return OpenCL event
     */
    inline const cl_event& getEvent() const {return _event;}

protected:
    /** @brief Set the variable type
     * @param tname New type
     */
    void type(const std::string& tname) {_typename = tname;}

    /** @brief Wait for variable underlying event to be complete
     * 
     * This function is tracking the syncing state, avoiding calling
     * clWaitForEvents() more than once, while the event is not updated.
     *
     * This is obviously a synced function, which is blocking program execution
     * until the underlying event is complete
     */
    void sync();

private:
    /// Name of the variable
    std::string _name;

    /// Type of the variable
    std::string _typename;

    /// List of events affecting this variable
    cl_event _event;

    /// Shortcut to avoid calling the expensive OpenCL API
    bool _synced;
};

/** @class ScalarVariable Variable.h Variable.h
 * @brief A generic Scalar variable
 */
template <class T>
class ScalarVariable : public Variable
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    ScalarVariable(const std::string varname);

    /** @brief Destructor
     */
    ~ScalarVariable() {};

    /** @brief Report that the varaible is not an array
     *
     * @return false
     */
    inline const bool isArray() const {return false;}

    /** @brief Get the variable type size.
     * @return Variable type size (in bytes)
     */
    inline const size_t typesize() const {return sizeof(T);}

    /** @brief Get variable pointer basis pointer
     *
     * This is a synced operation, that will block the program execution until
     * the underlying variable event is complete.
     *
     * @return Implementation pointer
     */
    inline const void* get(){sync(); return &_value;}

    /** @brief Set variable from memory
     *
     * This is a synced operation, that will block the program execution until
     * the underlying variable event is complete.
     *
     * @param ptr Memory to copy
     */
    inline void set(const void* ptr){sync(); memcpy(&_value, ptr, sizeof(T));}

    /** @brief Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    virtual const std::string asString() = 0;
protected:
    /// Variable value
    T _value;
};

/** @class ScalarVariable Variable.h Variable.h
 * @brief A generic Scalar variable of one single component
 */
template <class T>
class ScalarNumberVariable : public ScalarVariable<T>, public PyCast<T>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    ScalarNumberVariable(const std::string varname);

    /** @brief Destructor
     */
    ~ScalarNumberVariable() {};

    /** @brief Get a PyLongObject interpretation of the variable value
     * @param i0 ignored parameter.
     * @param n ignored parameter.
     * @return PyLongObject Python object.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyLongObject object.
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** @brief Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    virtual const std::string asString();
};

/** @brief An integer variable
 */
using IntVariable = ScalarNumberVariable<int>;

/** @brief An integer variable
 */
using UIntVariable = ScalarNumberVariable<unsigned int>;

/** @brief An integer variable
 */
using FloatVariable = ScalarNumberVariable<float>;

/** @class ScalarVecVariable Variable.h Variable.h
 * @brief A generic Scalar variable of 2 or more components
 *
 * Scalar varaibles are those that are stored in the host, while array variables
 * are conversely stored in the computational device
 */
template <class T>
class ScalarVecVariable : public ScalarVariable<T>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable
     * @param vartype Type of the variable
     * @param dims Number of components
     */
    ScalarVecVariable(const std::string varname,
                      const unsigned int dims);

    /** @brief Destructor
     */
    ~ScalarVecVariable() {};

    /** @brief Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    virtual const std::string asString();
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
    const bool checkPyhonObjectDims(PyObject* obj);
private:
    /// Number of components of the type
    unsigned int _dims;
};

/** @class Vec2Variable Variable.h Variable.h
 * @brief A vec2 variable
 */
class Vec2Variable : public ScalarVecVariable<vec2>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    Vec2Variable(const std::string varname);

    /** @brief Destructor
     */
    ~Vec2Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class Vec3Variable Variable.h Variable.h
 * @brief A vec3 variable
 */
class Vec3Variable : public ScalarVecVariable<vec3>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    Vec3Variable(const std::string varname);

    /** @brief Destructor
     */
    ~Vec3Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class Vec4Variable Variable.h Variable.h
 * @brief A vec4 variable
 */
class Vec4Variable : public ScalarVecVariable<vec4>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    Vec4Variable(const std::string varname);

    /** @brief Destructor
     */
    ~Vec4Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_FLOAT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_FLOAT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class IVec2Variable Variable.h Variable.h
 * @brief A ivec2 variable
 */
class IVec2Variable : public ScalarVecVariable<ivec2>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    IVec2Variable(const std::string varname);

    /** @brief Destructor
     */
    ~IVec2Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class IVec3Variable Variable.h Variable.h
 * @brief A ivec3 variable
 */
class IVec3Variable : public ScalarVecVariable<ivec3>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    IVec3Variable(const std::string varname);

    /** @brief Destructor
     */
    ~IVec3Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class IVec4Variable Variable.h Variable.h
 * @brief A ivec4 variable
 */
class IVec4Variable : public ScalarVecVariable<ivec4>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    IVec4Variable(const std::string varname);

    /** @brief Destructor
     */
    ~IVec4Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_INT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_INT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class UIVec2Variable Variable.h Variable.h
 * @brief A uivec2 variable
 */
class UIVec2Variable : public ScalarVecVariable<uivec2>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    UIVec2Variable(const std::string varname);

    /** @brief Destructor
     */
    ~UIVec2Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class IVec3Variable Variable.h Variable.h
 * @brief A uivec3 variable
 */
class UIVec3Variable : public ScalarVecVariable<uivec3>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    UIVec3Variable(const std::string varname);

    /** @brief Destructor
     */
    ~UIVec3Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

/** @class UIVec4Variable Variable.h Variable.h
 * @brief A uivec4 variable
 */
class UIVec4Variable : public ScalarVecVariable<uivec4>
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     */
    UIVec4Variable(const std::string varname);

    /** @brief Destructor
     */
    ~UIVec4Variable() {};

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return PyArrayObject Python object (PyArray_UINT subtype).
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object (PyArray_UINT subtype).
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);
};

#ifdef HAVE_3D
    /** @brief A real vector variable.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    using VecVariable = Vec4Variable;
    /** @brief A real vector variable.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    using IVecVariable = IVec4Variable;
    /** @brief A real vector variable.
     *
     * The number of components depends on weather the 2D version or 3D
     * version is compiled:
     *   - 2D = 2 components
     *   - 3D = 4 components
     */
    using UIVecVariable = UIVec4Variable;
#else
    using VecVariable = Vec2Variable;
    using IVecVariable = IVec2Variable;
    using UIVecVariable = UIVec2Variable;
#endif

/** @class FloatVariable Variable.h Variable.h
 * @brief A float variable.
 */
class ArrayVariable : public Variable
{
public:
    /** @brief Constructor
     * @param varname Name of the variable.
     * @param vartype Type of the variable.
     */
    ArrayVariable(const std::string varname, const std::string vartype);

    /** @brief Destructor
     */
    ~ArrayVariable();

    /** @brief Report that the varaible is an array
     *
     * @return true
     */
    inline const bool isArray() const {return true;}

    /** Get the cl_mem type size.
     * @return cl_mem type size (in bytes)
     * @note In order to know the typesize of the components into the array you
     * may use Variables::typeToBytes()
     */
    inline const size_t typesize() const {return sizeof(cl_mem);}

    /** Get the array size.
     * @return Array allocated memory (in bytes)
     * @note In order to get the length of the array the command
     * size() / Variables::typeToBytes(type()) can be used
     */
    const size_t size() const;

    /** Get the OpenCL buffer memory pointer
     * @return Implementation pointer
     */
    inline const void* get(){return &_value;}

    /** Set variable from memory
     * @param ptr Memory to copy to the OpenCL memory buffer
     */
    inline void set(const void* ptr){_value = *(cl_mem*)ptr;}

    /** @brief Get a PyArrayObject interpretation of the variable
     * @param i0 First component to be read.
     * @param n Number of component to be read, 0 to read all available memory,
     * i.e. All the array after i0.
     * @return PyArrayObject Python object. NULL if the memory cannot be read.
     */
    PyObject* getPythonObject(int i0=0, int n=0);

    /** @brief Set the variable value from a Python object
     * @param obj PyArrayObject object.
     * @param i0 ignored parameter.
     * @param n ignored parameter
     * @return false if all gone right, true otherwise.
     */
    const bool setFromPythonObject(PyObject* obj, int i0=0, int n=0);

    /** Get the variable text representation
     * @return The variable represented as a string, NULL in case of errors.
     */
    const std::string asString();

    /** Get a component text representation
     * @param i Index of the component to be extracted.
     * @return The component represented as a string, NULL in case of errors.
     */
    const std::string asString(size_t i);
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
class Variables
{
public:
    /** @brief Constructor
     */
    Variables();

    /** @brief Destructor
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
    void registerVariable(const std::string& name,
                          const std::string& type,
                          const std::string& length,
                          const std::string& value);

    /** Get a variable.
     * @param index Index of the variable.
     * @return Variable, NULL if the variable cannot be found.
     */
    Variable* get(const unsigned int& index) const;

    /** Get a variable.
     * @param name Name of the variable.
     * @return Variable, NULL if the variable cannot be found.
     */
    Variable* get(const std::string& name) const;

    /** Get all the registered variables.
     * @return Variable, NULL if the variable cannot be found.
     */
    inline const std::vector<Variable*>& getAll() const {return _vars;}

    /** Get the number of variables.
     * @return Number of registered variables.
     */
    inline const unsigned int size() const {return _vars.size();}

    /** Get the allocated memory.
     * @return Allocated memory on device. Just the arrays can contribute to this value.
     */
    const size_t allocatedMemory() const;

    /** Convert a type name to bytes.
     * @param type Type name.
     * @return Type size in bytes, 0 if the type is not recognized.
     */
    static const size_t typeToBytes(const std::string& type);

    /** Get the number of components of a type name.
     * @param type Type name.
     * @return Number of components (1 for not recognized types).
     */
    static const unsigned int typeToN(const std::string& type);

    /** Get if two types strings are the same one.
     * @param type_a First type name.
     * @param type_b Second type name.
     * @param ignore_asterisk true to ignore the asterisk of arrays.
     * @return true if it is the same type, false otherwise.
     */
    static const bool isSameType(const std::string& type_a,
                                 const std::string& type_b,
                                 const bool ignore_asterisk=true);

    /** Solve a string, interpreting the variables.
     * @param type_name Type of the output desired value.
     * @param value Expression to evaluate.
     * @param data Allocated memory where the result should be stored.
     * @param name Variable name to register in the tokenizer.
     * @note typeToBytes(type) bytes should be allocated in data.
     */
    void solve(const std::string& type_name,
               const std::string& value,
               void *data,
               const std::string name="NULL");

    /** @brief Populate all the variables in order that the tokenizer may get
     * the updated value
     */
    void populate();

    /** @brief Populate a variable in order that the tokenizer may know the
     * updated value
     * @param name Name of the variable to be populated, "" if all the
     * variables should be populated.
     */
    void populate(const std::string& name);

    /** @brief Populate a variable in order that the tokenizer may know the
     * updated value
     * @param var Variable to be populated.
     */
    void populate(Variable* var);

private:
    /** Register a scalar variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param value Variable value, NULL for arrays. It is optional for
     * scalar variables.
     */
    void registerScalar(const std::string& name,
                        const std::string& type,
                        const std::string& value);
    /** Register a cl_mem variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param length Array length, 1 for scalars, 0 for arrays that will
     * not be allocated at the start (for instance the heads of chains,
     * which requires the number of cells).
     */
    void registerClMem(const std::string& name,
                       const std::string& type,
                       const std::string& length);

    /** Read a set of components from a value array.
     * @param name Name of the variable. It is used to register variables in
     * the tokenizer and to report errors.
     * @param value Value string where the components are extracted from.
     * @param n Number of components to read.
     * @param v Allocated array where the components should be stored
     */
    void readComponents(const std::string& name,
                        const std::string& value,
                        const unsigned int& n,
                        float* v);


    /// Set of available variables
    std::vector<Variable*> _vars;
    /// Tokenizer to evaluate variables
    Tokenizer tok;
};

}}  // namespace

#endif // VARIABLE_H_INCLUDED
