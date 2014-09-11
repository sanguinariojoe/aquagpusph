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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <Variable.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

Variable::Variable(const char *varname, bool varsave)
    : _name(NULL)
    , _save(varsave)
{
    unsigned int len;

    len = strlen(varname) + 1;
    _name = new char[len];
    strcpy(_name, varname);
}

Variable::~Variable()
{
	delete[] _name; _name = NULL;
}

IntVariable::IntVariable(const char *varname, bool varsave)
    : Variable(varname, varsave)
    , _value(0)
{
}

IntVariable::~IntVariable()
{
}

UIntVariable::UIntVariable(const char *varname, bool varsave)
    : Variable(varname, varsave)
    , _value(0)
{
}

UIntVariable::~UIntVariable()
{
}

FloatVariable::FloatVariable(const char *varname, bool varsave)
    : Variable(varname, varsave)
    , _value(0.f)
{
}

FloatVariable::~FloatVariable()
{
}

Vec2Variable::Vec2Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0.f;
    _value.y = 0.f;
}

Vec2Variable::~Vec2Variable()
{
}

Vec3Variable::Vec3Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0.f;
    _value.y = 0.f;
    _value.z = 0.f;
}

Vec3Variable::~Vec3Variable()
{
}

Vec4Variable::Vec4Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0.f;
    _value.y = 0.f;
    _value.z = 0.f;
    _value.w = 0.f;
}

Vec4Variable::~Vec4Variable()
{
}

IVec2Variable::IVec2Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
}

IVec2Variable::~IVec2Variable()
{
}

IVec3Variable::IVec3Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
}

IVec3Variable::~IVec3Variable()
{
}

IVec4Variable::IVec4Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
    _value.w = 0;
}

IVec4Variable::~IVec4Variable()
{
}

UIVec2Variable::UIVec2Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
}

UIVec2Variable::~UIVec2Variable()
{
}

UIVec3Variable::UIVec3Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
}

UIVec3Variable::~UIVec3Variable()
{
}

UIVec4Variable::UIVec4Variable(const char *varname, bool varsave)
    : Variable(varname, varsave)
{
    _value.x = 0;
    _value.y = 0;
    _value.z = 0;
    _value.w = 0;
}

UIVec4Variable::~UIVec4Variable()
{
}

ArrayVariable::ArrayVariable(const char *varname, bool varsave)
    : Variable(varname, varsave)
    , _value(NULL)
{
}

ArrayVariable::~ArrayVariable()
{
    if(_value != NULL) clReleaseMemObject(_value); _value=NULL;
}

size_t ArrayVariable::size() const
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
        char msg[256];
        ScreenManager *S = ScreenManager::singleton();
        sprintf(msg,
                "Failure getting allocated memory from variable \"%s\"\n",
                name());
        S->addMessageF(3, msg);
	    S->printOpenCLError(status);
    }
    return memsize;
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
    for(i = 0; i < _vars.size(); i++){
        delete _vars.at(i);
    }
    _vars.clear();
}

bool Variables::registerVariable(const char* name,
                                 const char* type,
                                 const char* length,
                                 const char* value,
                                 const bool save)
{
    // Look for an already existing variable with the same name
    unsigned int i;
    for(i = 0; i < _vars.size(); i++){
        if(!strcmp(_vars.at(i)->name(), name)){
            delete _vars.at(i);
            _vars.erase(_vars.begin() + i);
        }
    }

    // Discriminate scalar vs. array
    if(strstr(type, "*")){
        return registerClMem(name, type, length, save);
    }
    else{
        return registerScalar(name, type, value, save);
    }
    return false;
}

Variable* Variables::get(unsigned int index)
{
    if(index >= _vars.size()){
        return NULL;
    }
    return _vars.at(index);
}

Variable* Variables::get(const char* name)
{
    unsigned int i;
    for(i = 0; i < _vars.size(); i++){
        if(!strcmp(name, _vars.at(i)->name())){
            return _vars.at(i);
        }
    }
    return NULL;
}

bool Variables::registerScalar(const char* name,
                               const char* type,
                               const char* value,
                               const bool save)
{
    if(!strcmp(type, "int")){
        IntVariable *var = new IntVariable(name, save);
        if(strcmp(value, "")){
            int val = (int)tok.solve(value);
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "unsigned int")){
        UIntVariable *var = new UIntVariable(name, save);
        if(strcmp(value, "")){
            unsigned int val = (unsigned int)tok.solve(value);
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "float")){
        FloatVariable *var = new FloatVariable(name, save);
        if(strcmp(value, "")){
            unsigned int val = (unsigned int)tok.solve(value);
            tok.registerVariable(name, (float)val);
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec")){
        VecVariable *var = new VecVariable(name, save);
        if(strcmp(value, "")){
            vec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = auxval[0];
                val.y = auxval[1];
                val.z = auxval[2];
                val.w = auxval[3];
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = auxval[0];
                val.y = auxval[1];
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec2")){
        Vec2Variable *var = new Vec2Variable(name, save);
        if(strcmp(value, "")){
            vec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec3")){
        Vec3Variable *var = new Vec3Variable(name, save);
        if(strcmp(value, "")){
            vec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "vec4")){
        Vec4Variable *var = new Vec4Variable(name, save);
        if(strcmp(value, "")){
            vec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = auxval[0];
            val.y = auxval[1];
            val.z = auxval[2];
            val.w = auxval[3];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec")){
        IVecVariable *var = new IVecVariable(name, save);
        if(strcmp(value, "")){
            ivec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = (int)auxval[0];
                val.y = (int)auxval[1];
                val.z = (int)auxval[2];
                val.w = (int)auxval[3];
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = (int)auxval[0];
                val.y = (int)auxval[1];
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec2")){
        IVec2Variable *var = new IVec2Variable(name, save);
        if(strcmp(value, "")){
            ivec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = (int)auxval[0];
            val.y = (int)auxval[1];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec3")){
        IVec3Variable *var = new IVec3Variable(name, save);
        if(strcmp(value, "")){
            ivec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = (int)auxval[0];
            val.y = (int)auxval[1];
            val.z = (int)auxval[2];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "ivec4")){
        IVec4Variable *var = new IVec4Variable(name, save);
        if(strcmp(value, "")){
            ivec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = (int)auxval[0];
            val.y = (int)auxval[1];
            val.z = (int)auxval[2];
            val.w = (int)auxval[3];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec")){
        UIVecVariable *var = new UIVecVariable(name, save);
        if(strcmp(value, "")){
            uivec val;
            #ifdef HAVE_3D
                float auxval[4];
                if(readComponents(name, value, 4, auxval))
                    return true;
                val.x = (unsigned int)auxval[0];
                val.y = (unsigned int)auxval[1];
                val.z = (unsigned int)auxval[2];
                val.w = (unsigned int)auxval[3];
            #else
                float auxval[2];
                if(readComponents(name, value, 2, auxval))
                    return true;
                val.x = (unsigned int)auxval[0];
                val.y = (unsigned int)auxval[1];
            #endif // HAVE_3D
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec2")){
        UIVec2Variable *var = new UIVec2Variable(name, save);
        if(strcmp(value, "")){
            uivec2 val;
            float auxval[2];
            if(readComponents(name, value, 2, auxval))
                return true;
            val.x = (unsigned int)auxval[0];
            val.y = (unsigned int)auxval[1];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec3")){
        UIVec3Variable *var = new UIVec3Variable(name, save);
        if(strcmp(value, "")){
            uivec3 val;
            float auxval[3];
            if(readComponents(name, value, 3, auxval))
                return true;
            val.x = (unsigned int)auxval[0];
            val.y = (unsigned int)auxval[1];
            val.z = (unsigned int)auxval[2];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else if(!strcmp(type, "uivec4")){
        UIVec4Variable *var = new UIVec4Variable(name, save);
        if(strcmp(value, "")){
            uivec4 val;
            float auxval[4];
            if(readComponents(name, value, 4, auxval))
                return true;
            val.x = (unsigned int)auxval[0];
            val.y = (unsigned int)auxval[1];
            val.z = (unsigned int)auxval[2];
            val.w = (unsigned int)auxval[3];
            var->set(&val);
        }
        _vars.push_back(var);
    }
    else{
        char msg[256];
        ScreenManager *S = ScreenManager::singleton();
        sprintf(msg,
                "\"%s\" declared as \"%s\", wich is not a valid scalar type.\n",
                name,
                type);
        S->addMessageF(3, msg);
        S->addMessageF(0, "Valid types are:\n");
        S->addMessageF(0, "\tint\n");
        S->addMessageF(0, "\tunsigned int\n");
        S->addMessageF(0, "\tfloat\n");
        S->addMessageF(0, "\tvec\n");
        S->addMessageF(0, "\tvec2\n");
        S->addMessageF(0, "\tvec3\n");
        S->addMessageF(0, "\tvec4\n");
        S->addMessageF(0, "\tivec\n");
        S->addMessageF(0, "\tivec2\n");
        S->addMessageF(0, "\tivec3\n");
        S->addMessageF(0, "\tivec4\n");
        S->addMessageF(0, "\tuivec\n");
        S->addMessageF(0, "\tuivec2\n");
        S->addMessageF(0, "\tuivec3\n");
        S->addMessageF(0, "\tuivec4\n");
        return true;
    }
    return false;
}

bool Variables::registerClMem(const char* name,
                              const char* type,
                              const char* length,
                              const bool save)
{
    return false;
}

bool Variables::readComponents(const char* name,
                               const char* value,
                               unsigned int n,
                               float* v)
{
    float val;
    unsigned int i;
    char msg[256];
    ScreenManager *S = ScreenManager::singleton();
    if(n == 0){
        sprintf(msg,
                "%u components required for the variable \"%s\".\n",
                n,
                name);
        S->addMessageF(3, msg);
    }
    if(n > 4){
        S->addMessageF(3, "No more than 4 components can be required\n");
        sprintf(msg,
                "%u components required for the variable \"%s\".\n",
                n,
                name);
        S->addMessageF(0, msg);
    }
    char remain[strlen(value) + 1];
    char aux[strlen(value) + 1];
    char nameaux[strlen(name) + 3];
    const char* extensions[4] = {".x", ".y", ".z", ".w"};
    strcpy(remain, value);
    for(i = 0; i < n - 1; i++){
        strcpy(aux, remain);
        if(!strchr(aux, ',')){
            sprintf(msg, "Failure reading the variable \"%s\" value", name);
            S->addMessageF(3, msg);
            sprintf(msg, "%u fields expected, %u received.\n", n, i);
            S->addMessageF(0, msg);
            return true;
        }
        strcpy(strchr(aux, ','), "");
        strcpy(remain, strchr(remain, ',') + 1);
        val = tok.solve(aux);
        strcpy(nameaux, name);
        strcat(nameaux, extensions[i]);
        tok.registerVariable(nameaux, val);
        v[i] = val;
    }

    strcpy(aux, remain);
    if(strchr(aux, ','))
        strcpy(strchr(aux, ','), "");
    val = tok.solve(aux);
    strcpy(nameaux, name);
    if(n != 1)
        strcat(nameaux, extensions[n - 1]);
    tok.registerVariable(nameaux, val);
    v[n - 1] = val;
    return false;
}

}}  // namespace

