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

/** @def PY_ARRAY_UNIQUE_SYMBOL
 * @brief Define the extension module which this Python stuff should be linked
 * to.
 *
 * In AQUAgpusph all the Python stuff is linked in the same group AQUA_ARRAY_API
 * @see
 * http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
 */
#define PY_ARRAY_UNIQUE_SYMBOL AQUA_ARRAY_API
#include <numpy/npy_no_deprecated_api.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>
#include <numpy/npy_3kcompat.h>
#include <sstream>

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Python.hpp"
#include "SetScalar.hpp"

/** @brief stdout Python redirector.
 * @see logMsg
 */
const char* _stdout_redirect = "         \n\
class stdoutWriter(object):              \n\
    def write(self, data):               \n\
        aquagpusph.log(0, data)          \n\
    def flush(self):                     \n\
        pass                             \n\
\n";

/** @brief stderr Python redirector.
 * @see logMsg
 */
const char* _stderr_redirect = "         \n\
class stderrWriter(object):              \n\
    def write(self, data):               \n\
        aquagpusph.log(0, data)          \n\
    def flush(self):                     \n\
        pass                             \n\
\n";

/** @brief Get a variable by its name.
 * @param self Module.
 * @param args Positional arguments.
 * @param keywds Keyword arguments.
 * @return Computed value, NULL if errors have been detected.
 */
static PyObject*
get(PyObject* self, PyObject* args, PyObject* keywds)
{
	auto C = Aqua::CalcServer::CalcServer::singleton();
	auto vars = C->variables();
	const char* varname;

	int i0 = 0;
	int n = 0;

	static char* kwlist[] = { "varname", "offset", "n", NULL };

	if (!PyArg_ParseTupleAndKeywords(
	        args, keywds, "s|ii", kwlist, &varname, &i0, &n)) {
		return NULL;
	}

	auto var = vars->get(varname);
	if (!var) {
		std::ostringstream errstr;
		errstr << "Variable \"" << varname << "\" has not been declared";
		PyErr_SetString(PyExc_ValueError, errstr.str().c_str());
		return NULL;
	}

	PyObject* result = var->getPythonObject(i0, n);
	return result;
}

/** @brief Set a variable by its name.
 * @param self Module.
 * @param args Positional arguments.
 * @param keywds Keyword arguments.
 * @return Computed value, NULL if errors have been detected.
 */
static PyObject*
set(PyObject* self, PyObject* args, PyObject* keywds)
{
	auto C = Aqua::CalcServer::CalcServer::singleton();
	auto vars = C->variables();
	const char* varname;
	PyObject* value;

	int i0 = 0;
	int n = 0;

	static char* kwlist[] = { "varname", "value", "offset", "n", NULL };

	if (!PyArg_ParseTupleAndKeywords(
	        args, keywds, "sO|ii", kwlist, &varname, &value, &i0, &n)) {
		return NULL;
	}

	auto var = vars->get(varname);
	if (!var) {
		std::ostringstream errstr;
		errstr << "Variable \"" << varname << "\" has not been declared";
		PyErr_SetString(PyExc_ValueError, errstr.str().c_str());
		return NULL;
	}

	if (var->setFromPythonObject(value, i0, n)) {
		return NULL;
	}

	// Populate the variable if it is a scalar one
	if (var->type().find('*') == std::string::npos) {
		try {
			vars->populate(var);
		} catch (...) {
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
static PyObject*
logMsg(PyObject* self, PyObject* args, PyObject* keywds)
{
	int level;
	const char* msg;

	static char* kwlist[] = { "log_level", "message", NULL };

	if (!PyArg_ParseTupleAndKeywords(
	        args, keywds, "is", kwlist, &level, &msg)) {
		return NULL;
	}

	switch (level) {
		case 0:
			LOG0(Aqua::L_DEBUG, msg);
			break;
		case 1:
			LOG0(Aqua::L_INFO, msg);
			break;
		case 2:
			LOG0(Aqua::L_WARNING, msg);
			break;
		case 3:
			LOG0(Aqua::L_ERROR, msg);
			break;
	}

	Py_RETURN_NONE;
}

/// List of methods declared in the module
static PyMethodDef methods[] = {
	{ "get", (PyCFunction)get, METH_VARARGS | METH_KEYWORDS, "Get a variable" },
	{ "set", (PyCFunction)set, METH_VARARGS | METH_KEYWORDS, "Set a variable" },
	{ "log",
	  (PyCFunction)logMsg,
	  METH_VARARGS | METH_KEYWORDS,
	  "Log a message" },
	{ NULL, NULL, 0, NULL }
};

// Python 3.0 or above
#if PY_VERSION_HEX >= 0x03000000

/// Module definition
static struct PyModuleDef moduledef = { PyModuleDef_HEAD_INIT,
	                                    "aquagpusph",
	                                    NULL,
	                                    -1,
	                                    methods,
	                                    NULL,
	                                    NULL,
	                                    NULL,
	                                    NULL };

/** @brief Module initialization.
 */
PyMODINIT_FUNC
PyInit_aquagpusph(void)
{
	PyObject* m;
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
PyMODINIT_FUNC
PyInit_aquagpusph(void)
{
	PyObject* m;

	m = Py_InitModule("aquagpusph", methods);
	if (m == NULL) {
		return;
	}

	import_array();
	import_umath();
}
#endif

namespace Aqua {
namespace CalcServer {

Python::Python(const std::string tool_name, const std::string script, bool once)
  : Tool(tool_name, once)
  , _script(script)
  , _module(NULL)
  , _func(NULL)
{
	// Look for a .py extension to remove it
	std::size_t last_sep = _script.find_last_of(".");
	if (last_sep != std::string::npos &&
	    !_script.substr(last_sep).compare(".py")) {
		_script = _script.substr(0, last_sep);
	}

	Profiler::substages({ new ScalarProfile("script", this) });
}

Python::~Python()
{
	if (_module)
		Py_DECREF(_module);
	_module = 0;
	if (_func)
		Py_DECREF(_func);
	_func = 0;
}

void
Python::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	initPython();
	load();
}

cl_event
Python::_execute(const std::vector<cl_event> events)
{
	PyObject* result;

	dynamic_cast<ScalarProfile*>(Profiler::substages().front())->start();

	result = PyObject_CallObject(_func, NULL);
	if (!result) {
		LOG(L_ERROR, "main() function execution failed.\n");
		LOG0(L_DEBUG, "\n--- Python report --------------------------\n\n");
		PyErr_Print();
		LOG0(L_DEBUG, "\n-------------------------- Python report ---\n\n");
		throw std::runtime_error("Python execution error");
	}

	if (!PyObject_TypeCheck(result, &PyBool_Type)) {
		LOG(L_ERROR, "main() function returned non boolean variable.\n");
		throw std::runtime_error("Python execution error");
	}

	if (result == Py_False) {
		LOG(L_ERROR, "main() function returned False.\n");
		throw std::runtime_error("Python invoked simulation stop");
	}

	dynamic_cast<ScalarProfile*>(Profiler::substages().front())->end();

	// This function is not pruducing events by itself. This work is relayed
	// to the setters
	return NULL;
}

void
Python::initPython()
{
	if (Py_IsInitialized()) {
		return;
	}

	PyImport_AppendInittab("aquagpusph", PyInit_aquagpusph);

	Py_Initialize();
	if (!Py_IsInitialized()) {
		LOG(L_ERROR, "Failure calling Py_Initialize().\n");
		throw std::runtime_error("Python error");
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
}

void
Python::load()
{
	PyObject* modName;

	std::ostringstream comm;
	PyRun_SimpleString("curdir = os.getcwd()");
	std::string path = trimCopy(getFolderFromFilePath(_script));
	std::string filename = trimCopy(getFileNameFromFilePath(_script));
	if (path.at(0) == '.') // "./" prefix has been set
		comm << "modulePath = curdir + \"" << path.substr(1) << "\"";
	else
		comm << "modulePath = \"" << path << "\"";
	PyRun_SimpleString(comm.str().c_str());
	PyRun_SimpleString("sys.path.append(modulePath)");

	modName = PyUnicode_FromString(filename.c_str());
	_module = PyImport_Import(modName);
	Py_DECREF(modName);
	modName = 0;
	if (!_module) {
		std::ostringstream msg;
		msg << "Python module \"" << filename << "\" cannot be imported."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		LOG0(L_DEBUG, "\n--- Python report --------------------------\n\n");
		PyErr_Print();
		LOG0(L_DEBUG, "\n-------------------------- Python report ---\n\n");
		throw std::runtime_error("Python execution error");
	}

	_func = PyObject_GetAttrString(_module, "main");
	if (!_func || !PyCallable_Check(_func)) {
		LOG(L_ERROR, "main() function cannot be found.\n");
		LOG0(L_DEBUG, "\n--- Python report --------------------------\n\n");
		PyErr_Print();
		LOG0(L_DEBUG, "\n-------------------------- Python report ---\n\n");
		throw std::runtime_error("Python execution error");
	}
}

}
} // namespaces
