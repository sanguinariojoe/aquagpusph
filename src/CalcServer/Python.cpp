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
 * @brief Python script execution tool.
 * (See Aqua::CalcServer::Python for details)
 */

#include <string.h>
#include <AuxiliarMethods.h>
#include <ScreenManager.h>
#include <CalcServer/Python.h>

/** @def PY_ARRAY_UNIQUE_SYMBOL
 * @brief Define the extension module which this Python stuff should be linked
 * to.
 *
 * In AQUAgpusph all the Python stuff is linked in the same group AQUA_ARRAY_API
 * @see http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define PY_ARRAY_UNIQUE_SYMBOL AQUA_ARRAY_API
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>

/** @brief stdout Python redirector.
 * @see logMsg
 */
const char* _stdout_redirect = "         \n\
class stdoutWriter(object):              \n\
    def write(self, data):               \n\
        aquagpusph.log(0, data.rstrip()) \n\
    def flush(self):                     \n\
        pass                             \n\
\n";

/** @brief stderr Python redirector.
 * @see logMsg
 */
const char* _stderr_redirect = "         \n\
class stderrWriter(object):              \n\
    def write(self, data):               \n\
        aquagpusph.log(3, data.rstrip()) \n\
    def flush(self):                     \n\
        pass                             \n\
\n";

/** @brief Get a variable by its name.
 * @param self Module.
 * @param args Positional arguments.
 * @param keywds Keyword arguments.
 * @return Computed value, NULL if errors have been detected.
 */
static PyObject* get(PyObject *self, PyObject *args, PyObject *keywds)
{
    Aqua::CalcServer::CalcServer *C = Aqua::CalcServer::CalcServer::singleton();
    Aqua::InputOutput::Variables *V = C->variables();
    const char* varname;

    int i0 = 0;
    int n = 0;

    static char *kwlist[] = {"varname", "offset", "n", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|ii", kwlist,
                                     &varname, &i0, &n)){
        return NULL;
    }

    Aqua::InputOutput::Variable *var = V->get(varname);
    if(!var){
        char errstr[64 + strlen(varname)];
        sprintf(errstr, "Variable \"%s\" has not been declared", varname);
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }

    PyObject *result = var->getPythonObject(i0, n);
    return result;
}

/** @brief Set a variable by its name.
 * @param self Module.
 * @param args Positional arguments.
 * @param keywds Keyword arguments.
 * @return Computed value, NULL if errors have been detected.
 */
static PyObject* set(PyObject *self, PyObject *args, PyObject *keywds)
{
    Aqua::CalcServer::CalcServer *C = Aqua::CalcServer::CalcServer::singleton();
    Aqua::InputOutput::Variables *V = C->variables();
    const char* varname;
    PyObject *value;

    int i0 = 0;
    int n = 0;

    static char *kwlist[] = {"varname", "value", "offset", "n", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "sO|ii", kwlist,
                                     &varname, &value, &i0, &n)){
        return NULL;
    }

    Aqua::InputOutput::Variable *var = V->get(varname);
    if(!var){
        char errstr[64 + strlen(varname)];
        sprintf(errstr, "Variable \"%s\" has not been declared", varname);
        PyErr_SetString(PyExc_ValueError, errstr);
        return NULL;
    }

    if(var->setFromPythonObject(value, i0, n)){
        return NULL;
    }

    // Populate the variable if it is a scalar one
    if(var->type().find('*') == std::string::npos){
        try {
            V->populate(var);
        } catch(...) {
            return NULL;
        }
    }

    Py_RETURN_NONE;
}

/** @brief Log a message from the Python.
 *
 * In AQUAgpusph the Python stdout and stderr are redirected to this function,
 * such that:
 *     - stdout messages will be logged with level 0
 *     - stderr messages will be logged with level 3
 * @param self Module.
 * @param args Positional arguments.
 * @param keywds Keyword arguments.
 * @return Computed value, NULL if errors have been detected.
 */
static PyObject* logMsg(PyObject *self, PyObject *args, PyObject *keywds)
{
    Aqua::InputOutput::ScreenManager *S =
        Aqua::InputOutput::ScreenManager::singleton();
    int level;
    const char* msg;

    static char *kwlist[] = {"log_level", "message", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "is", kwlist,
                                     &level, &msg)){
        return NULL;
    }

    switch(level) {
        case 0:
            S->addMessage(Aqua::L_DEBUG, msg); break;
        case 1:
            S->addMessage(Aqua::L_INFO, msg); break;
        case 2:
            S->addMessage(Aqua::L_WARNING, msg); break;
        case 3:
            S->addMessage(Aqua::L_ERROR, msg); break;            
    }

    Py_RETURN_NONE;
}

/// List of methods declared in the module
static PyMethodDef methods[] = {
    {"get", (PyCFunction)get, METH_VARARGS | METH_KEYWORDS, "Get a variable"},
    {"set", (PyCFunction)set, METH_VARARGS | METH_KEYWORDS, "Set a variable"},
    {"log", (PyCFunction)logMsg, METH_VARARGS | METH_KEYWORDS, "Log a message"},
    {NULL, NULL, 0, NULL}
};

// Python 3.0 or above
#if PY_VERSION_HEX >= 0x03000000

/// Module definition
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "aquagpusph",
    NULL,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

/** @brief Module initialization.
 */
PyMODINIT_FUNC PyInit_aquagpusph(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    return m;
}

// Python 2.7 or below
#else

/** @brief Module initialization.
 */
PyMODINIT_FUNC PyInit_aquagpusph(void)
{
    PyObject *m;

    m = Py_InitModule("aquagpusph", methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();
}
#endif

namespace Aqua{ namespace CalcServer{

Python::Python(const char *tool_name, const char *script)
    : Tool(tool_name)
    , _script(NULL)
    , _module(NULL)
    , _func(NULL)
{
    _script = new char[strlen(script) + 1];
    strcpy(_script, script);
    // Look for a .py extension to remove it
    char *dot = strrchr(_script, '.');
    if(dot && !strcmp(dot, ".py")){
        strcpy(dot, "");
    }
}

Python::~Python()
{
    if(_script) delete[] _script; _script=NULL;
    if(_module) Py_DECREF(_module); _module=0;
    if(_func) Py_DECREF(_func); _func=0;
}

bool Python::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    if(initPython())
        return true;

    if(load())
        return true;

    return false;
}

bool Python::_execute()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    PyObject *result;

    result = PyObject_CallObject(_func, NULL);
    if(!result) {
        S->addMessageF(L_ERROR, "main() function execution failed.\n");
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }

    if(!PyObject_TypeCheck(result, &PyBool_Type)){
        S->addMessageF(L_ERROR,
                       "main() function returned non boolean variable.\n");
        return true;
    }

    if(result == Py_False){
        S->addMessageF(L_ERROR, "main() function returned False.\n");
        return true;
    }

    return false;
}

bool Python::initPython()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if(Py_IsInitialized()){
        return false;
    }

    PyImport_AppendInittab("aquagpusph", PyInit_aquagpusph);

    Py_Initialize();
    if(!Py_IsInitialized()){
        S->addMessageF(L_ERROR, "Failure calling Py_Initialize().\n");
        return true;
    }

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("curdir = os.getcwd()");
    PyRun_SimpleString("sys.path.append(curdir)");

    PyImport_ImportModule("aquagpusph");

    PyRun_SimpleString("import aquagpusph");
    PyRun_SimpleString(_stdout_redirect);
    PyRun_SimpleString("logger = stdoutWriter()");
    PyRun_SimpleString("sys.stdout = logger");
    PyRun_SimpleString(_stderr_redirect);
    PyRun_SimpleString("logger = stderrWriter()");
    PyRun_SimpleString("sys.stderr = logger");

    return false;
}

bool Python::load()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char msg[1024];
    PyObject *modName;

    char comm[512];
    strcpy(comm, "");
    PyRun_SimpleString("curdir = os.getcwd()");
    const char *path = getFolderFromFilePath((const char*)_script);
    const char *filename = getFileNameFromFilePath((const char*)_script);
    if(path[0]=='.')   // "./" prefix has been set
        sprintf(comm, "modulePath = curdir + \"%s\"", &path[1]);
    else
        sprintf(comm, "modulePath = \"%s\"", path);
    PyRun_SimpleString(comm);
    PyRun_SimpleString("sys.path.append(modulePath)");

    modName = PyUnicode_FromString(filename);
    _module = PyImport_Import(modName);
    Py_DECREF(modName); modName=0;
    if(!_module){
        sprintf(msg,
                "Python module \"%s\" cannot be imported.\n",
                filename);
        S->addMessageF(L_ERROR, msg);
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }

    _func = PyObject_GetAttrString(_module, "main");
    if(!_func || !PyCallable_Check(_func)) {
        S->addMessageF(L_ERROR, "main() function cannot be found.\n");
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }
    return false;
}

}}  // namespaces
