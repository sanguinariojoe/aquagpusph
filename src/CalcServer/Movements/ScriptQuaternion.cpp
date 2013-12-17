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

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/ScriptQuaternion.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

// ----------------------------------------------------------------------------
// Include the Time manager
// ----------------------------------------------------------------------------
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the Screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include auxiliar methods
// ----------------------------------------------------------------------------
#include <AuxiliarMethods.h>

#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )
using namespace xercesc;

namespace Aqua{ namespace CalcServer{ namespace Movement{

ScriptQuaternion::ScriptQuaternion()
	: Quaternion()
	, mScript(0)
	, mTorque(0)
	, mModule(0)
	, mFunc(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	mScript = new char[256];
	if(!mScript){
	    S->addMessage(3, "(ScriptQuaternion::ScriptQuaternion) Can't allocate memory at host.\n");
	    exit(EXIT_FAILURE);
	}
	mTorque = new Torque();
}

ScriptQuaternion::~ScriptQuaternion()
{
	if(mScript) delete[] mScript;   mScript=0;
	if(mTorque) delete[] mTorque;   mTorque=0;
	if(mModule) Py_DECREF(mModule); mModule=0;
	if(mFunc)   Py_DECREF(mFunc);   mFunc=0;
}

bool ScriptQuaternion::execute()
{
	unsigned int i,j,nI,nJ;
	PyObject *arg, *quat;
	vec torque, force, cor;
	mat axis;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	CalcServer *C = CalcServer::singleton();
	//! Calculate Torque
	mTorque->cor(mCOR);
	mTorque->execute();
	torque = mTorque->torque();
	force = mTorque->force();
	//! Call to script
	arg = args(torque, force);
	quat = PyObject_CallObject(mFunc, arg);
	if(PyErr_Occurred() || !quat){
	    S->addMessage(3, "(ScriptQuaternion::execute): perform() execution fail.\n");
	    printf("\n--- Python report --------------------------\n\n");
	    PyErr_Print();
	    printf("\n-------------------------- Python report ---\n");
	    return true;
	}
	//! Retrieve result
	if(!isValidOutput(quat)){
	    S->addMessage(3, "(ScriptQuaternion::execute): perform() returned quaternion is not valid.\n");
	    return true;
	}
	#ifdef HAVE_3D
	    // [ [COR.x,COR.y,COR.z], [X.x,X.y,X.z], [Y.x,Y.y,Y.z], [Z.x,Z.y,Z.z] ]
	    nI = 4;
	    nJ = 3;
	    float array[4][3];
	#else
	    // [ [COR.x,COR.y], [X.x,X.y], [Y.x,Y.y] ]
	    nI = 3;
	    nJ = 2;
	    float array[3][2];
	#endif
	for(i=0;i<nI;i++){
	    PyObject *row = PyList_GetItem(quat, i);
	    for(j=0;j<nJ;j++){
	        PyObject *item = PyList_GetItem(row, j);
	        if(PyObject_TypeCheck(item,&PyLong_Type)){
	            array[i][j] = (float)PyLong_AsLong(item);
	        }
	        if(PyObject_TypeCheck(item,&PyFloat_Type)){
	            array[i][j] = (float)PyFloat_AsDouble(item);
	        }
	    }
	}
	Py_DECREF(quat); quat=0;
	Py_DECREF(arg); arg=0;
	#ifdef HAVE_3D
	    // [ [COR.x,COR.y,COR.z], [X.x,X.y,X.z], [Y.x,Y.y,Y.z], [Z.x,Z.y,Z.z] ]
	    cor.x = array[0][0]; cor.y = array[0][1]; cor.z = array[0][2]; cor.w = 0.f;
	    axis[0].x = array[1][0]; axis[0].y = array[1][1]; axis[0].z = array[1][2]; axis[0].w = 0.f;
	    axis[1].x = array[2][0]; axis[1].y = array[2][1]; axis[1].z = array[2][2]; axis[1].w = 0.f;
	    axis[2].x = array[3][0]; axis[2].y = array[3][1]; axis[2].z = array[3][2]; axis[2].w = 0.f;
	#else
	    // [ [COR.x,COR.y], [X.x,X.y], [Y.x,Y.y] ]
	    cor.x = array[0][0]; cor.y = array[0][1];
	    axis[0].x = array[1][0]; axis[0].y = array[1][1];
	    axis[1].x = array[2][0]; axis[1].y = array[2][1];
	#endif
	//! Move quaternion
	set(cor, axis);
	Quaternion::execute();
	return false;
}

bool ScriptQuaternion::_parse(xercesc::DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[1024];
	//! Read Script section
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("PyScript"));
	if(!nodes->getLength()){
	    S->addMessage(3, "(ScriptQuaternion::_parse): Python script has not been specified.\n");
	    S->addMessage(0, "\tPyScript tag is mandatory.\n");
	    return true;
	}
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
		//! Get script path
		strcpy( mScript, XMLString::transcode( elem->getAttribute(XMLString::transcode("file")) ) );
	    sprintf(msg, "(ScriptQuaternion::_parse): Using \"%s\" script file.\n", mScript);
	    S->addMessage(1, msg);
	}
	//! Get initial quaternion
	if(initQuaternion()){
	    return true;
	}
	//! Prepare python function
	if(initPython()){
	    return true;
	}
	return false;
}

bool ScriptQuaternion::initQuaternion()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[1024];
	const char* function = "init";
	PyObject *modName, *func, *quat;
	vec cor;
	mat axis;
	//! Initialize Python
	if(!Py_IsInitialized()){
	    Py_Initialize();
	}
	if(!Py_IsInitialized()){
	    S->addMessage(3, "(ScriptQuaternion::initQuaternion): Python can't be initializaed.\n");
	    return true;
	}
	//! Add script path to Python modules path
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("import os");
	PyRun_SimpleString("curdir = os.getcwd()");
	char comm[512];
	strcpy(comm, "");
	const char *path = getFolderFromFilePath(mScript);
	if(path[0]=='.')   // "./" has been append
	    sprintf(comm, "modulePath = curdir + \"%s\"", &path[1]);
	else
	    sprintf(comm, "modulePath = \"%s\"", path);
	PyRun_SimpleString(comm);
	PyRun_SimpleString("sys.path.append(modulePath)");
	//! Import module
	unsigned int start = strlen(path);
	if( (path[0]=='.') && (mScript[0]!='.') )   // "./" has been append
	    start -= 2;
	modName = PyUnicode_FromString(&mScript[start]); // Module name starts at dir ends
	mModule = PyImport_Import(modName);
	Py_DECREF(modName); modName=0;
	if(!mModule){
	    sprintf(msg, "(ScriptQuaternion::initQuaternion): Python module \"%s\" can't be imported.\n", &mScript[start]);
	    S->addMessage(3, msg);
	    printf("\n--- Python report --------------------------\n\n");
	    PyErr_Print();
	    printf("\n-------------------------- Python report ---\n");
	    return true;
	}
	//! Load function
	func = PyObject_GetAttrString(mModule, "init");
	if(!func || !PyCallable_Check(func)) {
	    S->addMessage(3, "(ScriptQuaternion::initQuaternion): Can't find init() function.\n");
	    printf("\n--- Python report --------------------------\n\n");
	    PyErr_Print();
	    printf("\n-------------------------- Python report ---\n");
	    return true;
	}
	//! Call function
	quat = PyObject_CallObject(func, NULL);
	if(!quat) {
	    S->addMessage(3, "(ScriptQuaternion::initQuaternion): init() execution fail.\n");
	    printf("\n--- Python report --------------------------\n\n");
	    PyErr_Print();
	    printf("\n-------------------------- Python report ---\n");
	    return true;
	}
	//! Retrieve result
	if(!isValidOutput(quat)){
	    S->addMessage(3, "(ScriptQuaternion::initQuaternion): init() returned quaternion is not valid.\n");
	    return true;
	}
	unsigned int i,j,nI,nJ;
	#ifdef HAVE_3D
	    // [ [COR.x,COR.y,COR.z], [X.x,X.y,X.z], [Y.x,Y.y,Y.z], [Z.x,Z.y,Z.z] ]
	    nI = 4;
	    nJ = 3;
	    float array[4][3];
	#else
	    // [ [COR.x,COR.y], [X.x,X.y], [Y.x,Y.y] ]
	    nI = 3;
	    nJ = 2;
	    float array[3][2];
	#endif
	for(i=0;i<nI;i++){
	    PyObject *row = PyList_GetItem(quat, i);
	    for(j=0;j<nJ;j++){
	        PyObject *item = PyList_GetItem(row, j);
	        if(PyObject_TypeCheck(item,&PyLong_Type)){
	            array[i][j] = (float)PyLong_AsLong(item);
	        }
	        if(PyObject_TypeCheck(item,&PyFloat_Type)){
	            array[i][j] = (float)PyFloat_AsDouble(item);
	        }
	    }
	}
	#ifdef HAVE_3D
	    // [ [COR.x,COR.y,COR.z], [X.x,X.y,X.z], [Y.x,Y.y,Y.z], [Z.x,Z.y,Z.z] ]
	    cor.x = array[0][0]; cor.y = array[0][1]; cor.z = array[0][2]; cor.w = 0.f;
	    axis[0].x = array[1][0]; axis[0].y = array[1][1]; axis[0].z = array[1][2]; axis[0].w = 0.f;
	    axis[1].x = array[2][0]; axis[1].y = array[2][1]; axis[1].z = array[2][2]; axis[1].w = 0.f;
	    axis[2].x = array[3][0]; axis[2].y = array[3][1]; axis[2].z = array[3][2]; axis[2].w = 0.f;
	#else
	    // [ [COR.x,COR.y], [X.x,X.y], [Y.x,Y.y] ]
	    cor.x = array[0][0]; cor.y = array[0][1];
	    axis[0].x = array[1][0]; axis[0].y = array[1][1];
	    axis[1].x = array[2][0]; axis[1].y = array[2][1];
	#endif
	set(cor, axis, true);
	//! Destroy all Python C object references.
	Py_DECREF(quat); quat=0;
	Py_DECREF(func); func=0;
	return false;
}

bool ScriptQuaternion::initPython()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	//! Load function
	mFunc = PyObject_GetAttrString(mModule, "perform");
	if(!mFunc || !PyCallable_Check(mFunc)) {
	    S->addMessage(3, "(ScriptQuaternion::initPython): Can't find perform() function.\n");
	    printf("\n--- Python report --------------------------\n\n");
	    PyErr_Print();
	    printf("\n-------------------------- Python report ---\n\n");
	    return true;
	}
	return false;
}

bool ScriptQuaternion::isValidOutput(PyObject *quat)
{
	unsigned int i,j,nI,nJ;
	// Check that is a tuple
	if(!PyList_Check(quat)){
	    return false;
	}
	#ifdef HAVE_3D
	    // [ [COR.x,COR.y,COR.z], [X.x,X.y,X.z], [Y.x,Y.y,Y.z], [Z.x,Z.y,Z.z] ]
	    nI = 4;
	    nJ = 3;
	#else
	    // [ [COR.x,COR.y], [X.x,X.y], [Y.x,Y.y] ]
	    nI = 3;
	    nJ = 2;
	#endif
	if(PyList_Size(quat) != nI){
	    return true;
	}
	for(i=0;i<nI;i++){
	    PyObject *row = PyList_GetItem(quat, i);
	    if(!PyList_Check(row)){
	        return false;
	    }
	    if(PyList_Size(row) != nJ){
	        return false;
	    }
	    for(j=0;j<nJ;j++){
	        PyObject *item = PyList_GetItem(row, j);
	        if((!PyObject_TypeCheck(item,&PyLong_Type)) &&
	           (!PyObject_TypeCheck(item,&PyFloat_Type)) ){
	            return true;
	        }
	    }
	}
	return true;
}

PyObject* ScriptQuaternion::args(vec torque, vec force)
{
	PyObject *arg, *tuple, *val;
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	#ifdef HAVE_3D
	    // [[COR], [X], [Y], [Z], [oldCOR], [oldX], [oldY], [oldZ], torque, time, dt]
	    arg = PyTuple_New(8);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)mCOR.x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mCOR.y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)mCOR.z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 0, tuple);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)mAxis[0].x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mAxis[0].y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)mAxis[0].z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 1, tuple);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)mAxis[1].x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mAxis[1].y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)mAxis[1].z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 2, tuple);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)mAxis[2].x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mAxis[2].y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)mAxis[2].z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 3, tuple);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)torque.x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)torque.y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)torque.z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 4, tuple);
	    tuple = PyTuple_New(3);
	    val = PyFloat_FromDouble((double)force.x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)force.y); PyTuple_SetItem(tuple, 1, val);
	    val = PyFloat_FromDouble((double)force.z); PyTuple_SetItem(tuple, 2, val);
	    PyTuple_SetItem(arg, 5, tuple);
	    unsigned int start = 6;
	#else
	    // [[COR], [X], [Y], [oldCOR], [oldX], [oldY], torque, time, dt]
	    arg = PyTuple_New(7);
	    tuple = PyTuple_New(2);
	    val = PyFloat_FromDouble((double)mCOR.x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mCOR.y); PyTuple_SetItem(tuple, 1, val);
	    PyTuple_SetItem(arg, 0, tuple);
	    tuple = PyTuple_New(2);
	    val = PyFloat_FromDouble((double)mAxis[0].x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mAxis[0].y); PyTuple_SetItem(tuple, 1, val);
	    PyTuple_SetItem(arg, 1, tuple);
	    tuple = PyTuple_New(2);
	    val = PyFloat_FromDouble((double)mAxis[1].x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)mAxis[1].y); PyTuple_SetItem(tuple, 1, val);
	    PyTuple_SetItem(arg, 2, tuple);
	    // In 2D torque only have one component.
	    val = PyFloat_FromDouble((double)torque.x); PyTuple_SetItem(arg, 3, val);
	    tuple = PyTuple_New(2);
	    val = PyFloat_FromDouble((double)force.x); PyTuple_SetItem(tuple, 0, val);
	    val = PyFloat_FromDouble((double)force.y); PyTuple_SetItem(tuple, 1, val);
	    PyTuple_SetItem(arg, 4, tuple);
	    unsigned int start = 5;
	#endif
	val = PyFloat_FromDouble((double)T->time()); PyTuple_SetItem(arg, start, val);
	val = PyFloat_FromDouble((double)T->dt()); PyTuple_SetItem(arg, start+1, val);
	return arg;
}


}}} // namespaces
