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

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/XMLUni.hpp>
#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )

using namespace xercesc;

#include <Input/XML.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

bool loadXML(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];

    sprintf(msg, "Loading fluid from XML file \"%s\"\n", path);
    S->addMessageF(1, msg);
    sprintf(msg, "\tParticles from %u to %u\n", i0, i0+n-1);
    S->addMessage(0, msg);
    sprintf(msg, "\tCan take some time, please wait...\n");
    S->addMessage(0, msg);
	// Try to open as an ascii file in order to know if the file already exist
	FILE *dummy=0; dummy = fopen(path, "r");
	if(!dummy){
        sprintf(msg, "File don't exist.\n");
        S->addMessageF(3, msg);
	    return true;
	}
	xercesc::XercesDOMParser *mParser = new XercesDOMParser;
	mParser->setValidationScheme( XercesDOMParser::Val_Never );
	mParser->setDoNamespaces( false );
	mParser->setDoSchema( false );
	mParser->setLoadExternalDTD( false );
	mParser->parse( path );

	DOMDocument* doc = mParser->getDocument();
	DOMElement* root = doc->getDocumentElement();
	if( !root ){
        sprintf(msg, "Empty XML file\n");
        S->addMessageF(3, msg);
	    return true;
	}

	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Particle"));
	unsigned int count=0;
	unsigned int progress=-1;
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    count++;
	    if(count > n){
            sprintf(msg, "File contains more particles than specified for this fluid.\n");
            S->addMessageF(2, msg);
            sprintf(msg, "\t%u particles will be discarded.\n", n-count);
            S->addMessage(0, msg);
	        break;
	    }
	    if(progress != count*100/n){
	        progress = count*100/n;
            sprintf(msg, "\t\t%u%%\n", progress);
            S->addMessage(0, msg);
	    }
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Position"));
	    bool data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->pos[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->pos[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->pos[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->pos[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!data_exist){
            sprintf(msg, "Particle %u position don't specified.\n", count-1);
            S->addMessageF(3, msg);
	        return true;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Normal"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->normal[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->normal[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->normal[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->normal[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!data_exist){
            sprintf(msg, "Particle %u normal don't specified.\n", count-1);
            S->addMessageF(3, msg);
	        return true;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Velocity"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->v[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->v[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->v[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->v[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!data_exist){
            sprintf(msg, "Particle %u position don't specified.\n", count-1);
            S->addMessageF(3, msg);
	        return true;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Force"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->f[count+i0-1].x = atof(xmlAttribute(s_elem, "x"));
	        F->f[count+i0-1].y = atof(xmlAttribute(s_elem, "y"));
	        #if HAVE_3D
	            F->f[count+i0-1].z = atof(xmlAttribute(s_elem, "z"));
	            F->f[count+i0-1].w = 0.f;
	        #endif
	    }
	    if(!data_exist){
	        F->f[count+i0-1].x = 0.f;
	        F->f[count+i0-1].y = 0.f;
	        #if HAVE_3D
	            F->f[count+i0-1].z = 0.f;
	            F->f[count+i0-1].w = 0.f;
	        #endif
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Density"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->dens[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
	        F->dens[count+i0-1] = refd;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("DensityRate"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->drdt[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
	        F->drdt[count+i0-1] = 0.f;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Mass"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->mass[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
            sprintf(msg, "Particle %u mass don't specified.\n", count-1);
            S->addMessageF(3, msg);
	        return true;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("KernelH"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->hp[count+i0-1] = atof(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
	        F->hp[count+i0-1] = h;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Imove"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->imove[count+i0-1] = atoi(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
	        F->imove[count+i0-1] = 1;
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Ifluid"));
	    data_exist = false;
	    for( XMLSize_t j=0; j<nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        data_exist = true;
	        F->ifluid[count+i0-1] = atoi(xmlAttribute(s_elem, "value"));
	    }
	    if(!data_exist){
	        F->ifluid[count+i0-1] = ifluid;
	    }
	}
	if(count < n){
        sprintf(msg, "Not enought particles in the file\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\t%u of %u particles readed\n", count, n);
        S->addMessage(0, msg);
	    return true;
	}
    sprintf(msg, "Fluid set!\n");
    S->addMessageF(1, msg);
	return false;
}

}}} // namespace
