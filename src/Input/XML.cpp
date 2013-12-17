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
#include <Input/XML.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// xerces XML parser stuff
// ----------------------------------------------------------------------------
//! @def xmlAttribute Short definition of a method that returns an attribute of an element as char*
#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )
using namespace xercesc;

namespace Aqua{ namespace InputOutput{ namespace Input{

int loadXML(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	//! 1st.- Open xml file
    sprintf(msg, "(loadXML): Loading fluid from XML file \"%s\"\n", path);
    S->addMessage(1, msg);
    sprintf(msg, "\tParticles from %u to %u\n", i0, i0+n-1);
    S->addMessage(0, msg);
    sprintf(msg, "\tCan take some time, please wait...\n");
    S->addMessage(0, msg);
	// Try to open as ascii file in order to know if file exist
	FILE *dummy=0; dummy = fopen(path, "r");
	if(!dummy){
        sprintf(msg, "(loadXML): File don't exist.\n");
        S->addMessage(3, msg);
	    return 1;
	}
	xercesc::XercesDOMParser *mParser = new XercesDOMParser;
	mParser->setValidationScheme( XercesDOMParser::Val_Never );
	mParser->setDoNamespaces( false );
	mParser->setDoSchema( false );
	mParser->setLoadExternalDTD( false );
	mParser->parse( path );
	//! 2nd.- Handling
	DOMDocument* doc = mParser->getDocument();
	DOMElement* root = doc->getDocumentElement();
	if( !root ){
        sprintf(msg, "(loadXML): Empty XML file\n");
        S->addMessage(3, msg);
	    return 2;
	}
	//! 3rd.- Look for particles
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Particle"));
	unsigned int count=0;
	unsigned int Percentage=-1;
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    count++;
	    if(count > n){
            sprintf(msg, "(loadXML): File contains more particles than specified for this fluid.\n");
            S->addMessage(2, msg);
            sprintf(msg, "\t%u particles will be discarded.\n", n-count);
            S->addMessage(0, msg);
	        break;
	    }
	    if(Percentage != count*100/n){
	        Percentage = count*100/n;
            sprintf(msg, "\t\t%u%%\n", Percentage);
            S->addMessage(0, msg);
	    }
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Position")); // Mandatory field
	    bool haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->pos[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->pos[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->pos[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->pos[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!haveMe){
            sprintf(msg, "(loadXML): Particle %u position don't specified.\n", count-1);
            S->addMessage(3, msg);
	        return 4;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Normal")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->normal[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->normal[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->normal[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->normal[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!haveMe){
            sprintf(msg, "(loadXML): Particle %u normal don't specified.\n", count-1);
            S->addMessage(3, msg);
	        return 4;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Velocity")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->v[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->v[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->v[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->v[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!haveMe){
            sprintf(msg, "(loadXML): Particle %u position don't specified.\n", count-1);
            S->addMessage(3, msg);
	        return 5;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Force"));
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->f[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->f[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->f[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->f[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!haveMe){
	        F->f[count+i0-1].x = 0.f;
	        F->f[count+i0-1].y = 0.f;
	        #if HAVE_3D
	            F->f[count+i0-1].z = 0.f;
	            F->f[count+i0-1].w = 0.f;
	        #endif
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Density"));
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->dens[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
	        F->dens[count+i0-1] = refd;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("DensityRate"));
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->drdt[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
	        F->drdt[count+i0-1] = 0.f;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Mass")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->mass[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
            sprintf(msg, "(loadXML): Particle %u mass don't specified.\n", count-1);
            S->addMessage(3, msg);
	        return 5;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("KernelH")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->hp[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
	        F->hp[count+i0-1] = h;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Imove")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->imove[count+i0-1] = atoi(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
	        F->imove[count+i0-1] = 1;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Ifluid")); // Mandatory field
	    haveMe = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveMe = true;
	        F->ifluid[count+i0-1] = atoi(xmlAttribute(s_elem, "value"));
	    }
	    if(!haveMe){
	        F->ifluid[count+i0-1] = ifluid;
	    }
	}
	if(count < n){
        sprintf(msg, "(loadXML): Not enought particles in the file\n");
        S->addMessage(3, msg);
        sprintf(msg, "\t%u of %u particles readed\n", count, n);
        S->addMessage(0, msg);
	    return 6;
	}
    sprintf(msg, "(loadXML): Fluid set!\n");
    S->addMessage(1, msg);
	return 0;
}

}}} // namespace
