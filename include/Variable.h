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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	Variable(const char *varname, const char *vartype, bool varsave=false);

	/** Destructor.
	 */
	~Variable();

	/** Name of the variable
     * @return The name of the variable
	 */
    const char* name() const {return (const char*)_name;}

    /** Get if the variable should be printed in output files.
     * @return true if the variable should be printed in output files, false
     * otherwise.
     */
    const bool save() const {return _save;}

    /** Set if the variable should be printed in output files.
     * @param varsave true if the variable should be printed in output files,
     * false otherwise.
     */
    void save(bool varsave) {_save = varsave;}

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
private:
	/// Name of the variable
	char* _name;

    /// Type of the variable
    char* _typename;

    /** true if the variable should be printed in output files, false
     * otherwise.
     */
    bool _save;
};

/** @class IntVariable Variable.h Variable.h
 * @brief An integer variable.
 */
class IntVariable : public Variable
{
public:
	/** Constructor.
	 * @param varname Name of the variable.
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	IntVariable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	UIntVariable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	FloatVariable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	Vec2Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	Vec3Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	Vec4Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	IVec2Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	IVec3Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	IVec4Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	UIVec2Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	UIVec3Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	UIVec4Variable(const char *varname, bool varsave=false);

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
	 * @param varsave true if the variable should be printed in output files,
	 * false otherwise.
	 */
	ArrayVariable(const char *varname, const char *vartype, bool varsave=false);

	/** Destructor.
	 */
	~ArrayVariable();

    /** Get the variable type size.
     * @return Variable type size (in bytes)
     */
    size_t typesize() const {return sizeof(cl_mem);}

    /** Get the variable type size.
     * @return Variable type size (in bytes)
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
private:
    /// Variable value
    cl_mem _value;
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
     * @param save true if the variable should be saved, false otherwise.
     * @return false if all gone right, true otherwise
     */
    bool registerVariable(const char* name,
                          const char* type,
                          const char* length,
                          const char* value,
                          const bool save);

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

    /** Get the number of variables.
     * @return Number of registered variables.
     */
    unsigned int size() const {return _vars.size();}

    /** Get the allocated memory.
     * @return Allocated memory on device. Just the arrays can contribute to this value.
     */
    size_t allocatedMemory();
private:
    /** Register a scalar variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param value Variable value, NULL for arrays. It is optional for
     * scalar variables.
     * @param save true if the variable should be saved, false otherwise.
     * @return false if all gone right, true otherwise
     */
    bool registerScalar(const char* name,
                        const char* type,
                        const char* value,
                        const bool save);
    /** Register a cl_mem variable
     * @param name Name of the variable.
     * @param type Type of the variable.
     * @param length Array length, 1 for scalars, 0 for arrays that will
     * not be allocated at the start (for instance the heads of chains,
     * which requires the number of cells).
     * @param save true if the variable should be saved, false otherwise.
     * @return false if all gone right, true otherwise
     */
    bool registerClMem(const char* name,
                       const char* type,
                       const char* length,
                       const bool save);

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
