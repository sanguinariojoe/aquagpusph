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

static std::string str_val;
static std::ostringstream pyerr;

Variable::Variable(const std::string varname, const std::string vartype)
    : _name(varname)
    , _typename(vartype)
    , _event(NULL)
    , _synced(true)
{
    cl_int err_code;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    // Create a dummy starting event for the queue
    _event = clCreateUserEvent(C->context(), &err_code);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure creating user event for \"" <<
               varname << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
    err_code = clSetUserEventStatus(_event, CL_COMPLETE);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure setting user event status for \"" <<
               varname << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void Variable::setEvent(cl_event event)
{
    cl_int err_code;
    // Forgive the former/predecessor event
    err_code = clReleaseEvent(_event);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure releasing the predecessor event for \"" <<
               name() << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
    // And get retained the current event
    err_code = clRetainEvent(event);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure reteaning the event for \"" <<
               name() << "\" variable." << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
    _event = event;
    _synced = false;
}

void Variable::sync()
{
    if(_synced)
        return;

    cl_int err_code;
    err_code = clWaitForEvents(1, &_event);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure syncing variable \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
    _synced = true;
}

template <class T>
ScalarVariable<T>::ScalarVariable(const std::string varname)
    : Variable(varname, "")
    , PyCast<T>()
{
    type(type_name<decltype(_value)>());
}

template <class T>
PyObject* ScalarVariable<T>::getPythonObject(int i0, int n)
{
    PyObject* obj = this->valToPython(*((T*)this->get()));
    if(!obj){
        pyerr.str("");
        pyerr << "Failure creating Python object for \"" << this->name()
            << "\" variable" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
    }
    return obj;
}

template <class T>
const bool ScalarVariable<T>::setFromPythonObject(PyObject* obj,
                                                        int i0,
                                                        int n)
{
    void *value = this->PythonToPtr(obj);
    if(!value){
        pyerr.str("");
        pyerr << "Failure casting Python object for \"" << this->name()
            << "\" variable" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }

    this->set(value);
    return false;
}

template <class T>
ScalarNumberVariable<T>::ScalarNumberVariable(const std::string varname)
    : ScalarVariable<T>(varname)
{
    this->_value = 0;
}

template <class T>
const std::string ScalarNumberVariable<T>::asString() {
    std::ostringstream msg;
    msg << *((T*)this->get());
    str_val = msg.str();
    return str_val;
}

template <class T, size_t N>
ScalarVecVariable<T, N>::ScalarVecVariable(const std::string varname)
    : ScalarVariable<T>(varname)
    , _dims(N)
{
    memset(&(this->_value), 0, sizeof(T));
}

template <class T, size_t N>
const std::string ScalarVecVariable<T, N>::asString() {
    T *val = (T*)(this->get());
    std::ostringstream msg;
    msg << "(";
    for (unsigned int i = 0; i < _dims; i++) {
        msg << val->s[i] << ",";
    }
    str_val = msg.str();
    str_val.back() = ')';
    return str_val;
}

ArrayVariable::ArrayVariable(const std::string varname, const std::string vartype)
    : Variable(varname, vartype)
{
    _value = NULL;
}

ArrayVariable::~ArrayVariable()
{
    for(auto object : _objects){
        if(object) Py_DECREF(object);
    }
    _objects.clear();
    for(auto data : _data){
        if(data) free(data);
    }
    _data.clear();
    if(_value) clReleaseMemObject(_value); _value=NULL;
}

const size_t ArrayVariable::size() const
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
        std::ostringstream msg;
        msg << "Failure getting allocated memory from variable \"" << name()
            << "\"." << std::endl,
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(status);
    }
    return memsize;
}

PyObject* ArrayVariable::getPythonObject(int i0, int n)
{
    if(i0 < 0){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" cannot handle \"offset\" lower than 0" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    if(n < 0){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" cannot handle \"n\" lower than 0" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
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
        pyerr.str("");
        pyerr << "Failure reading variable \"" << name()
            << "\" out of bounds" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    size_t len = memsize / typesize - offset;
    if(n != 0){
        len = n;
    }
    if(len == 0){
        pyerr.str("");
        pyerr << "0 bytes asked to be read from variable \"" << name()
            << "\"" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    if((offset + len) * typesize > memsize){
        pyerr.str("");
        pyerr << "Failure reading variable \"" << name()
            << "\" out of bounds" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    npy_intp dims[] = {static_cast<npy_intp>(len), components};
    // Get the appropiate type
    int pytype = PyArray_FLOAT;
    if(!type().compare("unsigned int") ||
       !type().compare("unsigned int*") ||
       !type().compare("uivec") ||
       !type().compare("uivec*")){
       pytype = PyArray_UINT;
    }
    else if(!type().compare("int") ||
            !type().compare("int*") ||
            !type().compare("ivec") ||
            !type().compare("ivec*")){
       pytype = PyArray_INT;
    }
    else if(!type().compare("float") ||
            !type().compare("float*") ||
            !type().compare("vec") ||
            !type().compare("vec*")){
       pytype = PyArray_FLOAT;
    }
    else{
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" is of type \"" << type()
            << "\", which can't be handled by Python" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    // Reallocate memory
    void *data = malloc(len * typesize);
    if(!data){
        pyerr.str("");
        pyerr << "Failure allocating " << len * typesize
            << " bytes for variable \"" << name()
            << "\"" << std::endl;
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
    if(err_code != CL_SUCCESS){
        pyerr.str("");
        pyerr << "Failure downloading variable \"" << name()
            << "\"" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    // Build and return the Python object
    PyObject *obj = PyArray_SimpleNewFromData(2, dims, pytype, data);
    if(!obj){
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

const bool ArrayVariable::setFromPythonObject(PyObject* obj, int i0, int n)
{
    if(i0 < 0){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" cannot handle \"offset\" lower than 0" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return NULL;
    }
    if(n < 0){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" cannot handle \"n\" lower than 0" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
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
        pyerr.str("");
        pyerr << "Failure writing variable \"" << name()
            << "\" out of bounds" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }
    size_t len = memsize / typesize - offset;
    if(n != 0){
        len = n;
    }
    if(len == 0){
        pyerr.str("");
        pyerr << "0 bytes asked to be written into variable \"" << name()
            << "\"" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }
    if((offset + len) * typesize > memsize){
        pyerr.str("");
        pyerr << "Failure writing variable \"" << name()
            << "\" out of bounds" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }

    if(!PyObject_TypeCheck(obj, &PyArray_Type)){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" expected a PyArrayObject" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }

    PyArrayObject* array_obj = (PyArrayObject*) obj;
    if(array_obj->nd != 2){
        pyerr.str("");
        pyerr << "Variable \"" << name()
            << "\" expected a 2D array" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }
    npy_intp *dims = array_obj->dimensions;
    if((size_t)dims[0] != len){
        pyerr.str("");
        pyerr << len << " elements have been asked to be written in variable \""
              << name() << "\" but " << (size_t)dims[0]
              << " have been provided" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }
    if((size_t)dims[1] != components){
        pyerr.str("");
        pyerr << components << " components per elements are expected by variable \""
              << name() << "\" but " << (size_t)dims[1]
              << " have been provided" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }

    void *data = array_obj->data;
    cl_event event, event_wait = getEvent();
    err_code =  clEnqueueWriteBuffer(C->command_queue(),
                                     _value,
                                     CL_FALSE,
                                     offset * typesize,
                                     len * typesize,
                                     data,
                                     1,
                                     &event_wait,
                                     &event);
    if(err_code != CL_SUCCESS){
        pyerr.str("");
        pyerr << "Failure uploading variable \""
              << name() << "\"" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }
    setEvent(event);
    err_code = clReleaseEvent(event);
    if(err_code != CL_SUCCESS){
        pyerr.str("");
        pyerr << "Failure releasing variable \""
              << name() << "\" event" << std::endl;
        PyErr_SetString(PyExc_ValueError, pyerr.str().c_str());
        return true;
    }

    return false;
}

const std::string ArrayVariable::asString()
{
    const cl_mem* val = (cl_mem*)get();
    std::ostringstream msg;
    msg << val;
    str_val = msg.str();
    return str_val;
}

const std::string ArrayVariable::asString(size_t i)
{
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    size_t type_size = Variables::typeToBytes(type());
    size_t length = size() / type_size;
    if(i > length){
        std::ostringstream msg;
        msg << "Out of bounds (length = " << length << ")" << std::endl;
        LOG(L_ERROR, msg.str());
        return NULL;
    }
    void *ptr = malloc(type_size);
    if(!ptr){
        std::ostringstream msg;
        msg << "Failure allocating memory to download the variable \""
            << name() << "\"" << std::endl;
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
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure downloading the variable \"" << name() << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        Logger::singleton()->printOpenCLError(err_code);
        free(ptr);
        return NULL;
    }

    std::ostringstream str_stream;

    if(!type().compare("unsigned int*")){
        str_stream << ((unsigned int*)ptr)[0];
    }
    else if(!type().compare("uivec2*")){
        str_stream << "(" << ((unsigned int*)ptr)[0] << ","
                          << ((unsigned int*)ptr)[1] << ")";
    }
    else if(!type().compare("uivec3*")){
        str_stream << "(" << ((unsigned int*)ptr)[0] << ","
                          << ((unsigned int*)ptr)[1] << ","
                          << ((unsigned int*)ptr)[2] << ")";
    }
    else if(!type().compare("uivec4*")){
        str_stream << "(" << ((unsigned int*)ptr)[0] << ","
                          << ((unsigned int*)ptr)[1] << ","
                          << ((unsigned int*)ptr)[2] << ","
                          << ((unsigned int*)ptr)[3] << ")";
    }
    else if(!type().compare("uivec*")){
        #ifdef HAVE_3D
            str_stream << "(" << ((unsigned int*)ptr)[0] << ","
                              << ((unsigned int*)ptr)[1] << ","
                              << ((unsigned int*)ptr)[2] << ","
                              << ((unsigned int*)ptr)[3] << ")";
        #else
            str_stream << "(" << ((unsigned int*)ptr)[0] << ","
                              << ((unsigned int*)ptr)[1] << ")";
        #endif
    }
    else if(!type().compare("int*")){
        str_stream << ((int*)ptr)[0];
    }
    else if(!type().compare("ivec2*")){
        str_stream << "(" << ((int*)ptr)[0] << ","
                          << ((int*)ptr)[1] << ")";
    }
    else if(!type().compare("ivec3*")){
        str_stream << "(" << ((int*)ptr)[0] << ","
                          << ((int*)ptr)[1] << ","
                          << ((int*)ptr)[2] << ")";
    }
    else if(!type().compare("ivec4*")){
        str_stream << "(" << ((int*)ptr)[0] << ","
                          << ((int*)ptr)[1] << ","
                          << ((int*)ptr)[2] << ","
                          << ((int*)ptr)[3] << ")";
    }
    else if(!type().compare("ivec*")){
        #ifdef HAVE_3D
            str_stream << "(" << ((int*)ptr)[0] << ","
                              << ((int*)ptr)[1] << ","
                              << ((int*)ptr)[2] << ","
                              << ((int*)ptr)[3] << ")";
        #else
            str_stream << "(" << ((int*)ptr)[0] << ","
                              << ((int*)ptr)[1] << ")";
        #endif
    }
    else if(!type().compare("float*")){
        str_stream << ((float*)ptr)[0];
    }
    else if(!type().compare("vec2*")){
        str_stream << "(" << ((float*)ptr)[0] << ","
                          << ((float*)ptr)[1] << ")";
    }
    else if(!type().compare("vec3*")){
        str_stream << "(" << ((float*)ptr)[0] << ","
                          << ((float*)ptr)[1] << ","
                          << ((float*)ptr)[2] << ")";
    }
    else if(!type().compare("vec4*")){
        str_stream << "(" << ((float*)ptr)[0] << ","
                          << ((float*)ptr)[1] << ","
                          << ((float*)ptr)[2] << ","
                          << ((float*)ptr)[3] << ")";
    }
    else if(!type().compare("vec*")){
        #ifdef HAVE_3D
            str_stream << "(" << ((float*)ptr)[0] << ","
                              << ((float*)ptr)[1] << ","
                              << ((float*)ptr)[2] << ","
                              << ((float*)ptr)[3] << ")";
        #else
            str_stream << "(" << ((float*)ptr)[0] << ","
                              << ((float*)ptr)[1] << ")";
        #endif
    }
    else{
        std::ostringstream msg;
        msg << "Variable \"" << name()
            << "\" has unknown type \"" << type()
            << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        free(ptr);
        return NULL;
    }
    free(ptr);
    
    str_val = str_stream.str();
    return str_val;
}

void ArrayVariable::cleanMem()
{
    for(int i = _objects.size() - 1; i >= 0; i--){
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
    for(auto var : _vars){
        delete var;
    }
    _vars.clear();
}

void Variables::registerVariable(const std::string& name,
                                 const std::string& type,
                                 const std::string& length,
                                 const std::string& value)
{
    // Look for an already existing variable with the same name
    for(unsigned int i = 0; i < _vars.size(); i++){
        if(!_vars.at(i)->name().compare(name)){
            delete _vars.at(i);
            _vars.erase(_vars.begin() + i);
        }
    }

    // Discriminate scalar vs. array
    if(type.find('*') != std::string::npos){
        registerClMem(name, type, length);
    }
    else{
        registerScalar(name, type, value);
    }
}

Variable* Variables::get(const unsigned int& index) const
{
    if(index >= _vars.size()){
        return NULL;
    }
    return _vars.at(index);
}

Variable* Variables::get(const std::string& name) const
{
    for(auto var : _vars){
        if(!name.compare(var->name())){
            return var;
        }
    }
    return NULL;
}

const size_t Variables::allocatedMemory() const{
    size_t allocated_mem = 0;
    for(auto var : _vars){
        if(var->type().find('*') == std::string::npos){
            continue;
        }
        allocated_mem += var->size();
    }
    return allocated_mem;
}

const size_t Variables::typeToBytes(const std::string& type)
{
    const unsigned int n = typeToN(type);
    size_t type_size = 0;

    if(type.find("unsigned int") != std::string::npos ||
       type.find("uivec") != std::string::npos) {
        type_size = sizeof(unsigned int);
    }
    else if(type.find("int") != std::string::npos ||
            type.find("ivec") != std::string::npos){
        type_size = sizeof(int);
    }
    else if(type.find("float") != std::string::npos ||
            type.find("vec") != std::string::npos ||
            type.find("matrix") != std::string::npos){
        type_size = sizeof(float);
    }
    else{
        std::ostringstream msg;
        msg << "Unvalid type \"" << type << "\"" << std::endl;
        LOG(L_ERROR, msg.str());
        return 0;
    }
    return n * type_size;
}

const unsigned int Variables::typeToN(const std::string& type)
{
    unsigned int n = 1;
    if(type.find("vec2") != std::string::npos) {
        n = 2;
    }
    else if(type.find("vec3") != std::string::npos) {
        n = 3;
    }
    else if(type.find("vec4") != std::string::npos) {
        n = 4;
    }
    else if(type.find("vec") != std::string::npos) {
        #ifdef HAVE_3D
            n = 4;
        #else
            n = 2;
        #endif // HAVE_3D
    }
    else if(type.find("matrix") != std::string::npos) {
        #ifdef HAVE_3D
            n = 16;
        #else
            n = 4;
        #endif // HAVE_3D
    }
    return n;
}

const bool Variables::isSameType(const std::string& type_a,
                                 const std::string& type_b,
                                 const bool ignore_asterisk)
{
    if(typeToN(type_a) != typeToN(type_b)){
        return false;
    }

    if(!ignore_asterisk){
        if((type_a.find('*') != std::string::npos) &&
           (type_b.find('*') == std::string::npos)) {
            return false;
        }
        else if((type_a.find('*') == std::string::npos) &&
                (type_b.find('*') != std::string::npos)) {
            return false;
        }
    }

    std::string ta = trimCopy(type_a);
    if((ta.back() == '*')){
        ta.pop_back();
    }
    if((ta.back() == '2') ||
       (ta.back() == '3') ||
       (ta.back() == '4')){
        ta.pop_back();
    }
    std::string tb = trimCopy(type_b);
    if((tb.back() == '*')){
        tb.pop_back();
    }
    if((tb.back() == '2') ||
       (tb.back() == '3') ||
       (tb.back() == '4')){
        tb.pop_back();
    }

    if(ta.compare(tb)){
        return false;
    }

    return true;
}

void Variables::solve(const std::string& type_name,
                      const std::string& value,
                      void *data,
                      const std::string name)
{
    size_t typesize = typeToBytes(type_name);
    if(!typesize){
        throw std::runtime_error("0 bytes size");
    }
    if(!value.compare("")){
        LOG(L_ERROR, "Empty value received\n");
        throw std::runtime_error("Empty value string");
    }

    // Ignore whether it is an array or a scalar
    std::string type = trimCopy(type_name);
    if(type.back() == '*')
        type.pop_back();

    if(!type.compare("int")){
        int val;
        float auxval;
        readComponents(name, value, 1, &auxval);
        val = round(auxval);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("unsigned int")){
        unsigned int val;
        float auxval;
        readComponents(name, value, 1, &auxval);
        val = (unsigned int)round(auxval);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("float")){
        float val;
        readComponents(name, value, 1, &val);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("vec")){
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
    }
    else if(!type.compare("vec2")){
        vec2 val;
        float auxval[2];
        readComponents(name, value, 2, auxval);
        val.x = auxval[0];
        val.y = auxval[1];
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("vec3")){
        vec3 val;
        float auxval[3];
        readComponents(name, value, 3, auxval);
        val.x = auxval[0];
        val.y = auxval[1];
        val.z = auxval[2];
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("vec4")){
        vec4 val;
        float auxval[4];
        readComponents(name, value, 4, auxval);
        val.x = auxval[0];
        val.y = auxval[1];
        val.z = auxval[2];
        val.w = auxval[3];
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("ivec")){
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
    }
    else if(!type.compare("ivec2")){
        ivec2 val;
        float auxval[2];
        readComponents(name, value, 2, auxval);
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("ivec3")){
        ivec3 val;
        float auxval[3];
        readComponents(name, value, 3, auxval);
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        val.z = round(auxval[2]);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("ivec4")){
        ivec4 val;
        float auxval[4];
        readComponents(name, value, 4, auxval);
        val.x = round(auxval[0]);
        val.y = round(auxval[1]);
        val.z = round(auxval[2]);
        val.w = round(auxval[3]);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("uivec")){
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
    }
    else if(!type.compare("uivec2")){
        uivec2 val;
        float auxval[2];
        readComponents(name, value, 2, auxval);
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("uivec3")){
        uivec3 val;
        float auxval[3];
        readComponents(name, value, 3, auxval);
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        val.z = (unsigned int)round(auxval[2]);
        memcpy(data, &val, typesize);
    }
    else if(!type.compare("uivec4")){
        uivec4 val;
        float auxval[4];
        readComponents(name, value, 4, auxval);
        val.x = (unsigned int)round(auxval[0]);
        val.y = (unsigned int)round(auxval[1]);
        val.z = (unsigned int)round(auxval[2]);
        val.w = (unsigned int)round(auxval[3]);
        memcpy(data, &val, typesize);
    }
    else{
        throw std::runtime_error("Invalid variable type");
    }
}

void Variables::populate()
{
    for(auto var : _vars){
        populate(var);
    }
}

void Variables::populate(const std::string& name)
{
    Variable *var = get(name);
    if(!var){
        std::ostringstream msg;
        msg << "Variable \"" << name << "\" cannot be found" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("No such variable");
    }
    populate(var);
}

void Variables::populate(Variable* var)
{
    std::ostringstream name;
    const std::string type = trimCopy(var->type());
    if(!type.compare("int")){
        int val = *(int*)var->get();
        tok.registerVariable(var->name(), (float)val);
    }
    else if(!type.compare("unsigned int")){
        unsigned int val = *(unsigned int*)var->get();
        tok.registerVariable(var->name(), (float)val);
    }
    else if(!type.compare("float")){
        float val = *(float*)var->get();
        tok.registerVariable(var->name(), (float)val);
    }
    else if(!type.compare("vec")){
        vec val = *(vec*)var->get();
        #ifdef HAVE_3D
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
            name.str(""); name << var->name() << "_z";
            tok.registerVariable(name.str(), (float)(val.z));
            name.str(""); name << var->name() << "_w";
            tok.registerVariable(name.str(), (float)(val.w));
        #else
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
        #endif // HAVE_3D
    }
    else if(!type.compare("vec2")){
        vec2 val = *(vec2*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
    }
    else if(!type.compare("vec3")){
        vec3 val = *(vec3*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
    }
    else if(!type.compare("vec4")){
        vec4 val = *(vec4*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
        name.str(""); name << var->name() << "_w";
        tok.registerVariable(name.str(), (float)(val.w));
    }
    else if(!type.compare("ivec")){
        ivec val = *(ivec*)var->get();
        #ifdef HAVE_3D
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
            name.str(""); name << var->name() << "_z";
            tok.registerVariable(name.str(), (float)(val.z));
            name.str(""); name << var->name() << "_w";
            tok.registerVariable(name.str(), (float)(val.w));
        #else
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
        #endif // HAVE_3D
    }
    else if(!type.compare("ivec2")){
        ivec2 val = *(ivec2*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
    }
    else if(!type.compare("ivec3")){
        ivec3 val = *(ivec3*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
    }
    else if(!type.compare("ivec4")){
        ivec4 val = *(ivec4*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
        name.str(""); name << var->name() << "_w";
        tok.registerVariable(name.str(), (float)(val.w));
    }
    else if(!type.compare("uivec")){
        uivec val = *(uivec*)var->get();
        #ifdef HAVE_3D
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
            name.str(""); name << var->name() << "_z";
            tok.registerVariable(name.str(), (float)(val.z));
            name.str(""); name << var->name() << "_w";
            tok.registerVariable(name.str(), (float)(val.w));
        #else
            name.str(""); name << var->name() << "_x";
            tok.registerVariable(name.str(), (float)(val.x));
            name.str(""); name << var->name() << "_y";
            tok.registerVariable(name.str(), (float)(val.y));
        #endif // HAVE_3D
    }
    else if(!type.compare("uivec2")){
        uivec2 val = *(uivec2*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
    }
    else if(!type.compare("uivec3")){
        uivec3 val = *(uivec3*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
    }
    else if(!type.compare("uivec4")){
        uivec4 val = *(uivec4*)var->get();
        name.str(""); name << var->name() << "_x";
        tok.registerVariable(name.str(), (float)(val.x));
        name.str(""); name << var->name() << "_y";
        tok.registerVariable(name.str(), (float)(val.y));
        name.str(""); name << var->name() << "_z";
        tok.registerVariable(name.str(), (float)(val.z));
        name.str(""); name << var->name() << "_w";
        tok.registerVariable(name.str(), (float)(val.w));
    }
    else{
        std::ostringstream msg;
        msg << "\"" << var->name() << "\" declared as \""
            << type << "\", which is not a scalar type" << std::endl;
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

void Variables::registerScalar(const std::string& name,
                               const std::string& type_name,
                               const std::string& value)
{
    const std::string type = trimCopy(type_name);
    if(!type.compare("int")){
        IntVariable *var = new IntVariable(name);
        if(value.compare("")){
            int val = round(tok.solve(value));
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("unsigned int")){
        UIntVariable *var = new UIntVariable(name);
        if(value.compare("")){
            unsigned int val = (unsigned int)round(tok.solve(value));
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("float")){
        FloatVariable *var = new FloatVariable(name);
        if(value.compare("")){
            float val = tok.solve(value);
            tok.registerVariable(name, val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("vec")){
        VecVariable *var = new VecVariable(name);
        if(value.compare("")){
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
    }
    else if(!type.compare("vec2")){
        Vec2Variable *var = new Vec2Variable(name);
        if(value.compare("")){
            vec2 val;
            float auxval[2];
            readComponents(name, value, 2, auxval);
            val.x = auxval[0];
            val.y = auxval[1];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("vec3")){
        Vec3Variable *var = new Vec3Variable(name);
        if(value.compare("")){
            vec3 val;
            float auxval[3];
            readComponents(name, value, 3, auxval);
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("vec4")){
        Vec4Variable *var = new Vec4Variable(name);
        if(value.compare("")){
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
    }
    else if(!type.compare("ivec")){
        IVecVariable *var = new IVecVariable(name);
        if(value.compare("")){
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
    }
    else if(!type.compare("ivec2")){
        IVec2Variable *var = new IVec2Variable(name);
        if(value.compare("")){
            ivec2 val;
            float auxval[2];
            readComponents(name, value, 2, auxval);
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("ivec3")){
        IVec3Variable *var = new IVec3Variable(name);
        if(value.compare("")){
            ivec3 val;
            float auxval[3];
            readComponents(name, value, 3, auxval);
            val.x = round(auxval[0]);
            val.y = round(auxval[1]);
            val.z = round(auxval[2]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("ivec4")){
        IVec4Variable *var = new IVec4Variable(name);
        if(value.compare("")){
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
    }
    else if(!type.compare("uivec")){
        UIVecVariable *var = new UIVecVariable(name);
        if(value.compare("")){
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
    }
    else if(!type.compare("uivec2")){
        UIVec2Variable *var = new UIVec2Variable(name);
        if(value.compare("")){
            uivec2 val;
            float auxval[2];
            readComponents(name, value, 2, auxval);
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("uivec3")){
        UIVec3Variable *var = new UIVec3Variable(name);
        if(value.compare("")){
            uivec3 val;
            float auxval[3];
            readComponents(name, value, 3, auxval);
            val.x = (unsigned int)round(auxval[0]);
            val.y = (unsigned int)round(auxval[1]);
            val.z = (unsigned int)round(auxval[2]);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!type.compare("uivec4")){
        UIVec4Variable *var = new UIVec4Variable(name);
        if(value.compare("")){
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
    }
    else{
        std::ostringstream msg;
        msg << "\"" << name << "\" declared as \""
            << type << "\", which is not a valid scalar type" << std::endl;
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

void Variables::registerClMem(const std::string& name,
                              const std::string& type_name,
                              const std::string& length)
{
    unsigned int n;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    // Get the type size
    std::string type = trimCopy(type_name);
    type.pop_back();  // Remove the asterisk
    const size_t typesize = typeToBytes(type);
    if(!typesize){
        std::ostringstream msg;
        msg << "\"" << name << "\" declared as \""
            << type << "*\", which is not a valid array type" << std::endl;
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
    if(!length.compare("")){
        std::ostringstream msg;
        msg << "No length specified for variable \"" << name
            << "\", declared as array (" << type << "*)" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid array length");        
    }

    try{
        n = (unsigned int)round(tok.solve(length));
    }catch(...){
        std::ostringstream msg;
        msg << "Failure evaluating variable \"" << name
            << "\" length" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid array length");                
    }

    if(!n){
        std::ostringstream msg;
        msg << "Variable \"" << name << "\" has null length" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid array length");                
    }

    // Generate the variable
    ArrayVariable *var = new ArrayVariable(name, trimCopy(type_name));

    // Allocate memory on device
    cl_int status;
    cl_mem mem;
    mem = clCreateBuffer(C->context(),
                            CL_MEM_READ_WRITE,
                            n * typesize,
                            NULL,
                            &status);
    if(status != CL_SUCCESS) {
        LOG(L_ERROR, "Allocation failure.\n");
        Aqua::InputOutput::Logger::singleton()->printOpenCLError(status);
        throw std::bad_alloc();
    }
    var->set(&mem);

    _vars.push_back(var);
}

void Variables::readComponents(const std::string& name,
                               const std::string& value,
                               const unsigned int& n,
                               float* v)
{
    float val;
    unsigned int i, j;
    if(n == 0){
        std::ostringstream msg;
        msg << n << " components required for the variable \""
            << name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("0 components variable registration");
    }
    if(n > 4){
        LOG(L_ERROR, "No more than 4 components can be required\n");
        std::ostringstream msg;
        msg << n << " components required for the variable \""
            << name << "\"." << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("5+ components variable registration");
    }

    // Replace all the commas outside functions by semicolons, to be taken into
    // account as separators
    std::string edited_val = value;
    int parenthesis_counter = 0;
    for (auto it = edited_val.begin(); it != edited_val.end(); ++it) {
        // We does not care about unbalanced parenthesis, muparser will do it
        if(*it == '(')
            parenthesis_counter++;
        else if(*it == ')')
            parenthesis_counter--;
        else if((*it == ',') && (parenthesis_counter == 0)){
            *it = ';';
        }        
    }

    // Split the expression by semicolons, and parse each one
    const char* extensions[4] = {"_x", "_y", "_z", "_w"};
    i = 0;
    std::istringstream f(edited_val);
    std::string s;
    while (getline(f, s, ';')) {
        try {
            val = tok.solve(s);
        }
        catch(...){
            std::ostringstream msg;
            msg << "parsing variable \""
                << name << "\"" << std::endl;
            LOG0(L_DEBUG, msg.str());
            throw std::runtime_error("Variable evaluation error");
        }
        std::ostringstream tok_name;
        if (n != 1) {
            tok_name << name << extensions[i];
        }
        else {
            tok_name << name;
        }
        tok.registerVariable(tok_name.str(), val);
        v[i] = val;
        if(++i == n) {
            break;
        }
    }
    if (i != n) {
        std::ostringstream msg;
        msg << "Failure parsing the variable \""
            << name << "\" value" << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << n << " fields expected, "
            << i << " received" << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid number of fields");
    }
}

}}  // namespace
