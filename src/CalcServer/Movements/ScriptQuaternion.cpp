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

#include <CalcServer/Movements/ScriptQuaternion.h>
#include <CalcServer.h>
#include <TimeManager.h>
#include <ScreenManager.h>
#include <AuxiliarMethods.h>

#include <vector>
#include <deque>
static std::deque<char*> cpp_str;
static std::deque<XMLCh*> xml_str;

static char *xmlTranscode(const XMLCh *txt)
{
    char *str = xercesc::XMLString::transcode(txt);
    cpp_str.push_back(str);
    return str;
}

static XMLCh *xmlTranscode(const char *txt)
{
    XMLCh *str = xercesc::XMLString::transcode(txt);
    xml_str.push_back(str);
    return str;
}

static void xmlClear()
{
    unsigned int i;
    for(i = 0; i < cpp_str.size(); i++){
        xercesc::XMLString::release(&cpp_str.at(i));
    }
    cpp_str.clear();
    for(i = 0; i < xml_str.size(); i++){
        xercesc::XMLString::release(&xml_str.at(i));
    }
    xml_str.clear();
}

#ifdef xmlS
    #undef xmlS
#endif // xmlS
#define xmlS(txt) xmlTranscode(txt)

#ifdef xmlAttribute
    #undef xmlAttribute
#endif
#define xmlAttribute(elem, att) xmlS( elem->getAttribute(xmlS(att)) )

#ifdef xmlHasAttribute
    #undef xmlHasAttribute
#endif
#define xmlHasAttribute(elem, att) elem->hasAttribute(xmlS(att))

using namespace xercesc;

namespace Aqua{ namespace CalcServer{ namespace Movement{

ScriptQuaternion::ScriptQuaternion()
    : Quaternion()
    , _script(0)
    , _torque(0)
    , _module(0)
    , _func(0)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    _script = new char[256];
    if(!_script){
        S->addMessageF(3, "Failure while allocating memory at host.\n");
        exit(EXIT_FAILURE);
    }
    _torque = new Torque();
}

ScriptQuaternion::~ScriptQuaternion()
{
    if(_script) delete[] _script; _script=0;
    if(_torque) delete _torque; _torque=0;
    if(_module) Py_DECREF(_module); _module=0;
    if(_func) Py_DECREF(_func); _func=0;
    Py_Finalize();
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

    _torque->cor(_cor);
    _torque->execute();
    torque = _torque->torque();
    force = _torque->force();

    arg = args(torque, force);
    quat = PyObject_CallObject(_func, arg);
    if(PyErr_Occurred() || !quat){
        S->addMessageF(3, "perform() execution fail.\n");
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }

    if(!isValidOutput(quat)){
        S->addMessageF(3, "perform() returned an invalid quaternion.\n");
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

    set(cor, axis);
    Quaternion::execute();
    return false;
}

bool ScriptQuaternion::_parse(xercesc::DOMElement *root)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char msg[1024];

    DOMNodeList* nodes = root->getElementsByTagName(xmlS("PyScript"));
    if(!nodes->getLength()){
        S->addMessageF(3, "None Python script has not been specified.\n");
        S->addMessage(0, "\tPyScript tag is mandatory.\n");
        xmlClear();
        return true;
    }
    for( XMLSize_t i=0; i<nodes->getLength();i++ ){
        DOMNode* node = nodes->item(i);
        if( node->getNodeType() != DOMNode::ELEMENT_NODE )
            continue;
        DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );

        strcpy( _script, xmlS( elem->getAttribute(xmlS("file")) ) );
        sprintf(msg, "Using \"%s\" script file.\n", _script);
        S->addMessageF(1, msg);
    }

    if(initQuaternion()){
        xmlClear();
        return true;
    }

    if(initPython()){
        xmlClear();
        return true;
    }

    xmlClear();
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

    if(!Py_IsInitialized()){
        Py_Initialize();
    }
    if(!Py_IsInitialized()){
        S->addMessageF(3, "Python cannot be initializaed.\n");
        return true;
    }

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("curdir = os.getcwd()");
    char comm[512];
    strcpy(comm, "");
    const char *path = getFolderFromFilePath(_script);
    if(path[0]=='.')   // "./" has been append
        sprintf(comm, "modulePath = curdir + \"%s\"", &path[1]);
    else
        sprintf(comm, "modulePath = \"%s\"", path);
    PyRun_SimpleString(comm);
    PyRun_SimpleString("sys.path.append(modulePath)");

    unsigned int start = strlen(path);
    if( (path[0]=='.') && (_script[0]!='.') )   // "./" has been append
        start -= 2;
    modName = PyUnicode_FromString(&_script[start]); // Module name starts at dir ends
    _module = PyImport_Import(modName);
    Py_DECREF(modName); modName=0;
    if(!_module){
        sprintf(msg, "Python module \"%s\" cannot be imported.\n", &_script[start]);
        S->addMessageF(3, msg);
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }

    func = PyObject_GetAttrString(_module, "init");
    if(!func || !PyCallable_Check(func)) {
        S->addMessageF(3, "init() function cannot be found.\n");
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }
    //! Call function
    quat = PyObject_CallObject(func, NULL);
    if(!quat) {
        S->addMessageF(3, "init() function execution fail.\n");
        printf("\n--- Python report --------------------------\n\n");
        PyErr_Print();
        printf("\n-------------------------- Python report ---\n");
        return true;
    }

    if(!isValidOutput(quat)){
        S->addMessage(3, "init() returned an invalid quaternion.\n");
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

    Py_DECREF(quat); quat=0;
    Py_DECREF(func); func=0;
    return false;
}

bool ScriptQuaternion::initPython()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    _func = PyObject_GetAttrString(_module, "perform");
    if(!_func || !PyCallable_Check(_func)) {
        S->addMessageF(3, "perform() function cannot be found.\n");
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
        val = PyFloat_FromDouble((double)_cor.x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_cor.y); PyTuple_SetItem(tuple, 1, val);
        val = PyFloat_FromDouble((double)_cor.z); PyTuple_SetItem(tuple, 2, val);
        PyTuple_SetItem(arg, 0, tuple);
        tuple = PyTuple_New(3);
        val = PyFloat_FromDouble((double)_axis[0].x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_axis[0].y); PyTuple_SetItem(tuple, 1, val);
        val = PyFloat_FromDouble((double)_axis[0].z); PyTuple_SetItem(tuple, 2, val);
        PyTuple_SetItem(arg, 1, tuple);
        tuple = PyTuple_New(3);
        val = PyFloat_FromDouble((double)_axis[1].x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_axis[1].y); PyTuple_SetItem(tuple, 1, val);
        val = PyFloat_FromDouble((double)_axis[1].z); PyTuple_SetItem(tuple, 2, val);
        PyTuple_SetItem(arg, 2, tuple);
        tuple = PyTuple_New(3);
        val = PyFloat_FromDouble((double)_axis[2].x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_axis[2].y); PyTuple_SetItem(tuple, 1, val);
        val = PyFloat_FromDouble((double)_axis[2].z); PyTuple_SetItem(tuple, 2, val);
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
        val = PyFloat_FromDouble((double)_cor.x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_cor.y); PyTuple_SetItem(tuple, 1, val);
        PyTuple_SetItem(arg, 0, tuple);
        tuple = PyTuple_New(2);
        val = PyFloat_FromDouble((double)_axis[0].x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_axis[0].y); PyTuple_SetItem(tuple, 1, val);
        PyTuple_SetItem(arg, 1, tuple);
        tuple = PyTuple_New(2);
        val = PyFloat_FromDouble((double)_axis[1].x); PyTuple_SetItem(tuple, 0, val);
        val = PyFloat_FromDouble((double)_axis[1].y); PyTuple_SetItem(tuple, 1, val);
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
