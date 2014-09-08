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
 * @brief Continuous C1 quaternion based motion.
 * (See Aqua::CalcServer::Movement::C1Quaternion for details)
 */

#include <CalcServer/Movements/C1Quaternion.h>
#include <CalcServer.h>
#include <TimeManager.h>
#include <ScreenManager.h>

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

C1Quaternion::C1Quaternion()
    : Quaternion()
    , _data(0)
{
    _data = new C1Interpolation();
}

C1Quaternion::~C1Quaternion()
{
    if(_data)delete _data; _data=0;
}

bool C1Quaternion::execute()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
    CalcServer *C = CalcServer::singleton();
    char msg[1025];
    vec cor;
    mat axis;

    float t = T->time();
    std::deque<float> data = _data->update(t);
    #ifdef HAVE_3D
        if(data.size() != 10){
            S->addMessageF(3, "Invalid data file line.\n");
            sprintf(msg,
                    "\tLine with %lu varaibles reached, 10 variables are mandatory at 3D cases.\n",
                    data.size());
            S->addMessage(0, msg);
            return true;
        }
        cor.x = data.at(1);
        cor.y = data.at(2);
        cor.z = data.at(3);
        cor.w = 0.f;
        axis[0].x = data.at(4);
        axis[0].y = data.at(5);
        axis[0].z = data.at(6);
        axis[0].w = 0.f;
        axis[1].x = data.at(7);
        axis[1].y = data.at(8);
        axis[1].z = data.at(9);
        axis[1].w = 0.f;
    #else
        if(data.size() != 5){
            S->addMessageF(3, "Invalid data file line.\n");
            sprintf(msg,
                    "\tLine with %lu varaibles reached, 5 variables are mandatory at 2D cases.\n",
                    data.size());
            S->addMessage(0, msg);
            return true;
        }
        cor.x = data.at(1);
        cor.y = data.at(2);
        axis[0].x = data.at(3);
        axis[0].y = data.at(4);
    #endif

    #ifdef HAVE_3D
        axis[2].x = axis[0].y*axis[1].z - axis[0].z*axis[1].y;
        axis[2].y = axis[0].z*axis[1].x - axis[0].x*axis[1].z;
        axis[2].z = axis[0].x*axis[1].y - axis[0].y*axis[1].x;
        axis[2].w = 0.f;
    #else
        axis[1].x = -axis[0].y;
        axis[1].y =  axis[0].x;
    #endif
    set(cor, axis);
    Quaternion::execute();
    return false;
}

bool C1Quaternion::_parse(xercesc::DOMElement *root)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char msg[1024];

    DOMNodeList* nodes = root->getElementsByTagName(xmlS("DataFile"));
    if(!nodes->getLength()){
        S->addMessageF(3, "The data file has not been specified.\n");
        S->addMessage(0, "\tDataFile tag is mandatory.\n");
        xmlClear();
        return true;
    }
    for( XMLSize_t i=0; i<nodes->getLength();i++ ){
        DOMNode* node = nodes->item(i);
        if( node->getNodeType() != DOMNode::ELEMENT_NODE )
            continue;
        DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
        // Open it
        sprintf(msg,
                "Using \"%s\" data file.\n",
                xmlS( elem->getAttribute(xmlS("file")) ));
        S->addMessageF(1, msg);
        if( !_data->open(xmlS( elem->getAttribute(xmlS("file")) )) ){
            xmlClear();
            return true;
        }
    }

    xmlClear();
    setInitialPositions();
    return false;
}

bool C1Quaternion::setInitialPositions()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();

    // Perform the process of the motion guided by the time step
    float t = 0.f;
    float dt = P->time_opts.dt0;
    while(t < T->time() - dt){
        fflush(stdout);
        _data->update(t);
        if(!dt){
            dt = P->SPH_opts.h / P->SPH_opts.cs * P->time_opts.courant;
        }
        t += dt;
    }
     std::deque<float> data = _data->update(T->time() - dt);

    vec cor;
    mat axis;
    #ifdef HAVE_3D
        cor.x = data[1]; axis[0].x = data[4]; axis[1].x = data[7];
        cor.y = data[2]; axis[0].y = data[5]; axis[1].y = data[8];
        cor.z = data[3]; axis[0].z = data[6]; axis[1].z = data[9];
        axis[2].x = axis[0].y*axis[1].z - axis[0].z*axis[1].y;
        axis[2].y = axis[0].z*axis[1].x - axis[0].x*axis[1].z;
        axis[2].z = axis[0].x*axis[1].y - axis[0].y*axis[1].x;
        cor.w = 0.f; axis[0].w = 0.f; axis[1].w = 0.f; axis[2].w = 0.f;
    #else
        cor.x = data[1]; axis[0].x = data[3];
        cor.y = data[2]; axis[0].y = data[4];
        axis[1].x = -axis[0].y;
        axis[1].y =  axis[0].x;
    #endif
    set(cor, axis, true);
    return false;
}

}}} // namespaces
