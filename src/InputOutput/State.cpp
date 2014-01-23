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
#include <string.h>
#include <map>
#include <limits>

#include <InputOutput/State.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <FileManager.h>
#include <TimeManager.h>
#include <CalcServer.h>

#ifdef xmlS
    #undef xmlS
#endif // xmlS
#define xmlS(txt) XMLString::transcode(txt)

#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) xmlS( elem->getAttribute(xmlS(att)) )

#ifdef xmlHasAttribute
	#undef xmlHasAttribute
#endif
#define xmlHasAttribute(elem, att) elem->hasAttribute(xmlS(att))

using namespace xercesc;
using namespace std;

namespace Aqua{ namespace InputOutput{

State::State()
    : _output_file(NULL)
{
    unsigned int i;
	ScreenManager *S = ScreenManager::singleton();
	// Start the XML parser
	try {
	    XMLPlatformUtils::Initialize();
	}
	catch( XMLException& e ){
	    char* message = xmlS(e.getMessage());
        S->addMessageF(3, "XML toolkit initialization error.\n");
	    char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
	    XMLString::release( &message );
	    exit(255);
	}
	catch( ... ){
        S->addMessageF(3, "XML toolkit initialization error.\n");
        S->addMessage(0, "\tUnhandled exception\n");
	    exit(255);
	}

	// Look ofr the first available file place
    i = 0;
    _output_file = new char[256];
    while(true){
        sprintf(_output_file, "AQUAgpusph.save.%u.xml", i);
        FILE *f = fopen(_output_file, "r");
        if(f){
            // The file already exist, look for another one
            i++;
            fclose(f);
            continue;
        }
        break;
    }
}

State::~State()
{
	unsigned int i;
	ScreenManager *S = ScreenManager::singleton();
	// Terminate Xerces
	try{
	    XMLPlatformUtils::Terminate();
	}
	catch( xercesc::XMLException& e ){
	    char* message = xercesc::xmlS( e.getMessage() );
        S->addMessageF(3, "XML toolkit exit error.\n");
	    char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
	    XMLString::release( &message );
	    exit(255);
	}
	catch( ... ){
        S->addMessageF(3, "XML toolkit exit error.\n");
        S->addMessage(0, "\tUnhandled exception\n");
	    exit(255);
	}

	delete[] _output_file;
}

bool State::save()
{
    return write(_output_file);
}

bool State::load()
{
	FileManager *F = FileManager::singleton();
    return parse(F->inputFile());
}

bool State::parse(const char* filepath)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024];
	strcpy(msg, "");

	sprintf(msg, "Parsing the XML file \"%s\"\n", filepath);
	S->addMessageF(1, msg);
	// Try to open as ascii file, just to know if the file already exist
	FILE *dummy=0; dummy = fopen(filepath, "r");
	if(!dummy){
	    S->addMessageF(3, "File inaccessible!\n");
	    return true;
	}
	fclose(dummy);
	XercesDOMParser *parser = new XercesDOMParser;
	parser->setValidationScheme(XercesDOMParser::Val_Never);
	parser->setDoNamespaces(false);
	parser->setDoSchema(false);
	parser->setLoadExternalDTD(false);
	parser->parse(filepath);
 	DOMDocument* doc = parser->getDocument();
	DOMElement* root = doc->getDocumentElement();
	if( !root ){
	    S->addMessageF(3, "Empty XML file\n");
	    return true;
	}
	// Look for XML files included to parse them before
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Include"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
	    const char* included_file = xmlS(elem->getAttribute(xmlS("file")));
	    if(parse(included_file)){
            return true;
        }
	}

	if(parseSettings(root))
	    return true;
	if(parseOpenCL(root))
	    return true;
	if(parseTiming(root))
        return true;
	if(parseTiming(root))
	    return true;
	if(parseSPH(root))
	    return true;
	if(parseFluid(root))
	    return true;
	if(parseSensors(root))
	    return true;
	if(parseMotions(root))
	    return true;
	if(parsePortals(root))
	    return true;
	if(parseGhostParticles(root))
	    return true;
	return false;
}

bool State::parseSettings(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Settings"));
	for(XMLSize_t i=0; i<nodes->getLength();i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Verbose"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        P->settings.verbose_level = atoi(xmlAttribute(s_elem, "level"));
	    }

	    s_nodes = elem->getElementsByTagName(xmlS("Device"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        P->settings.platform_id = atoi(xmlAttribute(s_elem, "platform"));
	        P->settings.device_id   = atoi(xmlAttribute(s_elem, "device"));
	        if(!strcmp("ALL", xmlAttribute(s_elem, "type")))
	            P->settings.device_type = CL_DEVICE_TYPE_ALL;
	        else if(!strcmp("CPU", xmlAttribute(s_elem, "type")))
	            P->settings.device_type = CL_DEVICE_TYPE_CPU;
	        else if(!strcmp("GPU", xmlAttribute(s_elem, "type")))
	            P->settings.device_type = CL_DEVICE_TYPE_GPU;
	        else if(!strcmp("ACCELERATOR", xmlAttribute(s_elem, "type")))
	            P->settings.device_type = CL_DEVICE_TYPE_ACCELERATOR;
	        else if(!strcmp("DEFAULT", xmlAttribute(s_elem, "type")))
	            P->settings.device_type = CL_DEVICE_TYPE_DEFAULT;
	        else{
	            sprintf(msg,
                        "Unknow \"%s\" type of device\n",
                        xmlAttribute(s_elem, "type"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tThe valid options are:\n");
	            S->addMessage(0, "\t\tALL\n");
	            S->addMessage(0, "\t\tCPU\n");
	            S->addMessage(0, "\t\tGPU\n");
	            S->addMessage(0, "\t\tACCELERATOR\n");
	            S->addMessage(0, "\t\tDEFAULT\n");
	            return true;
	        }
	    }
	}
	return false;
}

bool State::parseOpenCL(DOMElement *root)
{
	ProblemSetup *P = ProblemSetup::singleton();
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("OpenCL"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

        std::map<char*, char*> tags;
        std::map<char*, char*>::iterator it;
        tags["Predictor"] = P->OpenCL_kernels.predictor;
        tags["LinkList"] = P->OpenCL_kernels.link_list;
        tags["Rates"] = P->OpenCL_kernels.rates;
        tags["Corrector"] = P->OpenCL_kernels.corrector;
        tags["TimeStep"] = P->OpenCL_kernels.time_step;
        tags["Reduction"] = P->OpenCL_kernels.reduction;
        tags["RadixSort"] = P->OpenCL_kernels.radix_sort;
        tags["ElasticBounce"] = P->OpenCL_kernels.elastic_bounce;
        tags["DeLeffe"] = P->OpenCL_kernels.de_Leffe;
        tags["GhostParticles"] = P->OpenCL_kernels.ghost;
        tags["Shepard"] = P->OpenCL_kernels.shepard;
        tags["DensInterpolation"] = P->OpenCL_kernels.dens_int;
        tags["Torque"] = P->OpenCL_kernels.torque;
        tags["Energy"] = P->OpenCL_kernels.energy;
        tags["Bounds"] = P->OpenCL_kernels.bounds;
        tags["Domain"] = P->OpenCL_kernels.domain;
        tags["Portal"] = P->OpenCL_kernels.portal;

        for(it=tags.begin(); it!=tags.end(); it++){
            DOMNodeList* s_nodes = elem->getElementsByTagName(
                xmlS(it->first));
            for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
                DOMNode* s_node = s_nodes->item(j);
                if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                    continue;
                DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
                strcpy(it->second, xmlAttribute(s_elem, "file"));
            }
        }
	}
	return false;
}

bool State::parseTiming(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Timing"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
	    // Get options
	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Option"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        const char *name = xmlAttribute(s_elem, "name");
			if(!strcmp(name, "Start") || strcmp(name, "SimulationStart")){
	            P->time_opts.t0 = atof(xmlAttribute(s_elem, "value"));
			}

			else if(!strcmp(name, "End") || strcmp(name, "SimulationStop")){
			    const char *type = xmlAttribute(s_elem, "type");
				if(!strcmp(type, "Time") || !strcmp(type, "T")){
					P->time_opts.sim_end_mode =
                        P->time_opts.sim_end_mode | __TIME_MODE__;
					P->time_opts.sim_end_time =
                        atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "Steps") || !strcmp(type, "S")){
					P->time_opts.sim_end_mode =
                        P->time_opts.sim_end_mode | __ITER_MODE__;
					P->time_opts.sim_end_step =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "Frames") || !strcmp(type, "F")){
					P->time_opts.sim_end_mode =
                        P->time_opts.sim_end_mode | __FRAME_MODE__;
					P->time_opts.sim_end_frame =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg,
                            "Unknow simulation stop criteria \"%s\"\n",
                            type);
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tTime\n");
	                S->addMessage(0, "\t\tSteps\n");
	                S->addMessage(0, "\t\tFrames\n");
	                return true;
				}
			}

			else if(!strcmp(name, "TimeStep")){
				if(   !strcmp(xmlAttribute(s_elem, "value"), "Courant")
	               || !strcmp(xmlAttribute(s_elem, "value"), "Variable") ){
					P->time_opts.dt_mode = __DT_VARIABLE__;
				}
				else if(   !strcmp(xmlAttribute(s_elem, "value"), "Fix")
	                    || !strcmp(xmlAttribute(s_elem, "value"), "Fixed") ){
					P->time_opts.dt_mode = __DT_FIXCALCULATED__;
				}
				else {
					P->time_opts.dt_mode = __DT_FIX__;
					P->time_opts.dt = atof(xmlAttribute(s_elem, "value"));
					if(P->time_opts.dt <= 0.f){
	                    sprintf(msg,
                                "%g s is not a valid time step\n",
                                P->time_opts.dt);
	                    S->addMessageF(2, msg);
	                    S->addMessage(
                            0,
                            "\tAutomatic calculated time step will be used\n");
	                    P->time_opts.dt_mode = __DT_FIXCALCULATED__;
					}
				}
			}
			else if(!strcmp(name, "MinTimeStep")){
	            P->time_opts.dt_min = atof(xmlAttribute(s_elem, "value"));
	            if(P->time_opts.dt_min < 0.f){
	                sprintf(msg,
                            "The minimum time step is lower than 0 s\n");
	                S->addMessageF(2, msg);
	            }
			}
			else if(!strcmp(name, "Courant")){
	            P->time_opts.courant = atof(xmlAttribute(s_elem, "value"));
	            if(P->time_opts.courant <= 0.f){
	                sprintf(msg,
                            "The courant factor cannot be <= 0\n");
	                S->addMessageF(3, msg);
	                return true;
	            }
			}
			else if(!strcmp(name, "ClampVel")){
	            if(!strcmp(xmlAttribute(s_elem, "value"), "True") ||
	               !strcmp(xmlAttribute(s_elem, "value"), "true"))
	            {
	                P->time_opts.velocity_clamp = true;
	                if(P->time_opts.dt_min <= 0.f){
	                    sprintf(msg,
                                "Velocity clamping has been activated, but the minimum time step is %g s\n",
                                P->time_opts.dt_min);
	                    S->addMessageF(2, msg);
	                    S->addMessage(0,
                                      "\tVelocity clamping will not take effect\n");
	                }
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "False") ||
	                    !strcmp(xmlAttribute(s_elem, "value"), "false"))
	            {
	                P->time_opts.velocity_clamp = false;
	            }
	            else{
	                sprintf(msg,
                            "%s is not a valid velocity clamping option\n",
                            xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\ttrue\n");
	                S->addMessage(0, "\t\tfalse\n");
	                return true;
	            }
			}

			else if(!strcmp(name, "LogFile")){
			    const char *type = xmlAttribute(s_elem, "type");
				if(!strcmp(type, "No")){
					P->time_opts.log_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(type, "FPS")){
					P->time_opts.log_mode =
                        P->time_opts.log_mode | __FPS_MODE__;
					P->time_opts.log_fps =
                        atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "IPF")){
					P->time_opts.log_mode =
                        P->time_opts.log_mode | __IPF_MODE__;
					P->time_opts.log_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg,
                            "Unknow log file print criteria \"%s\"\n",
                            xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(name, "EnFile")){
			    const char *type = xmlAttribute(s_elem, "type");
				if(!strcmp(type, "No")){
					P->time_opts.energy_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(type, "FPS")){
					P->time_opts.energy_mode =
                        P->time_opts.energy_mode | __FPS_MODE__;
					P->time_opts.energy_fps =
                        atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "IPF")){
					P->time_opts.energy_mode =
                        P->time_opts.energy_mode | __IPF_MODE__;
					P->time_opts.energy_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg,
                            "Unknow energy report print criteria \"%s\"\n",
                            xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(name, "BoundsFile")){
			    const char *type = xmlAttribute(s_elem, "type");
				if(!strcmp(type, "No")){
					P->time_opts.bounds_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(type, "FPS")){
					P->time_opts.bounds_mode =
                        P->time_opts.bounds_mode | __FPS_MODE__;
					P->time_opts.bounds_fps =
                        atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "IPF")){
					P->time_opts.bounds_mode =
                        P->time_opts.bounds_mode | __IPF_MODE__;
					P->time_opts.bounds_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg,
                            "Unknow bounds file print criteria \"%s\"\n",
                            type);
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(name, "Output")){
			    const char *type = xmlAttribute(s_elem, "type");
				if(!strcmp(type, "No")){
					P->time_opts.output_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(type, "FPS")){
					P->time_opts.output_mode =
                        P->time_opts.output_mode | __FPS_MODE__;
					P->time_opts.output_fps =
                        atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(type, "IPF")){
					P->time_opts.output_mode =
                        P->time_opts.output_mode | __IPF_MODE__;
					P->time_opts.output_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg,
                            "Unknow output file print criteria \"%s\"\n",
                            type);
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}

			else {
	            sprintf(msg, "Unknow timing option \"%s\"\n", xmlAttribute(s_elem, "name"));
	            S->addMessageF(3, msg);
	            return true;
			}
	    }
	}
	return false;
}

bool State::parseSPH(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("SPH"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
	    // Get options
	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Option"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        const char *name = xmlAttribute(s_elem, "name");

			if(!strcmp(name, "g")){
				P->SPH_opts.g.x = atof(xmlAttribute(s_elem, "x"));
				P->SPH_opts.g.y = atof(xmlAttribute(s_elem, "y"));
	            #ifdef HAVE_3D
	                P->SPH_opts.g.z = atof(xmlAttribute(s_elem, "z"));
	            #endif
			}

			else if(!strcmp(name, "hfac")){
				P->SPH_opts.hfac = atof(xmlAttribute(s_elem, "value"));
			}

			else if(!strcmp(name, "deltar")){
				P->SPH_opts.deltar.x = atof(xmlAttribute(s_elem, "x"));
	            P->SPH_opts.deltar.y = atof(xmlAttribute(s_elem, "y"));
	            #ifdef HAVE_3D
	                P->SPH_opts.deltar.z = atof(xmlAttribute(s_elem, "z"));
	            #endif
			}

			else if(!strcmp(name, "cs")){
				P->SPH_opts.cs = atof(xmlAttribute(s_elem, "value"));
			}

			else if(!strcmp(name, "LLSteps")){
				P->SPH_opts.link_list_steps =
                    atoi(xmlAttribute(s_elem, "value"));
				if(P->SPH_opts.link_list_steps < 1){
	                S->addMessageF(2, "LLSteps minor than 1 (it will be ignored).\n");
	                P->SPH_opts.link_list_steps = 1;
				}
			}

			else if(!strcmp(name, "DensSteps")){
				P->SPH_opts.dens_int_steps =
                    atoi(xmlAttribute(s_elem, "value"));
			}

			else if(!strcmp(name, "DensBounds")){
			    if(xmlHasAttribute(s_elem,"min"))
                    P->SPH_opts.rho_min = atof(xmlAttribute(s_elem, "min"));
			    if(xmlHasAttribute(s_elem,"max"))
                    P->SPH_opts.rho_max = atof(xmlAttribute(s_elem, "max"));
			}

			else if(!strcmp(name, "Boundary")){
	            if(!strcmp(xmlAttribute(s_elem, "value"), "ElasticBounce")){
	                P->SPH_opts.boundary_type = 0;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "FixedParticles")){
	                P->SPH_opts.boundary_type = 1;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "DeLeffe")){
	                P->SPH_opts.boundary_type = 2;
	            }
	            else{
	                sprintf(msg,
                            "Unknow boundary condition \"%s\"\n",
                            xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tElasticBounce\n");
	                S->addMessage(0, "\t\tFixedParticles\n");
	                S->addMessage(0, "\t\tDeLeffe\n");
	                return true;
	            }
			}

			else if(!strcmp(name, "SlipCondition")){
	            if(!strcmp(xmlAttribute(s_elem, "value"), "NoSlip")){
	                P->SPH_opts.slip_condition = 0;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "FreeSlip")){
	                P->SPH_opts.slip_condition = 1;
	            }
	            else{
	                sprintf(msg,
                            "Unknow slip condition \"%s\".\n",
                            xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tNoSlip\n");
	                S->addMessage(0, "\t\tFreeSlip\n");
	                return true;
	            }
			}

			else if(!strcmp(name, "BoundElasticFactor")){
			    P->SPH_opts.elastic_factor =
                    atof(xmlAttribute(s_elem, "value"));
			}

			else if(!strcmp(name, "BoundDist")){
			    P->SPH_opts.elastic_dist =
                    fabs(atof(xmlAttribute(s_elem, "value")));
			}

			else if(!strcmp(name, "Shepard")){
			    if(!strcmp(xmlAttribute(s_elem, "value"), "None")){
	                P->SPH_opts.has_shepard = 0;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "Force")){
	                P->SPH_opts.has_shepard = 1;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "Dens")){
	                P->SPH_opts.has_shepard = 2;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "ForceDens")){
	                P->SPH_opts.has_shepard = 3;
			    }
			    else {
	                sprintf(msg,
                            "Unknow shepard application \"%s\".\n",
                            xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tNone\n");
	                S->addMessage(0, "\t\tForce\n");
	                S->addMessage(0, "\t\tDens\n");
	                S->addMessage(0, "\t\tForceDens\n");
	                return true;
			    }
			}

			else if(!strcmp(name, "Domain")){
			    P->SPH_opts.has_domain = true;
			    P->SPH_opts.domain_min.x = atof(xmlAttribute(s_elem, "x"));
			    P->SPH_opts.domain_min.y = atof(xmlAttribute(s_elem, "y"));
			    P->SPH_opts.domain_max.x = atof(xmlAttribute(s_elem, "l")) +
                    P->SPH_opts.domain_min.x;
	            #ifdef HAVE_3D
	                P->SPH_opts.domain_min.z =
                        atof(xmlAttribute(s_elem, "z"));
	                P->SPH_opts.domain_max.y =
                        atof(xmlAttribute(s_elem, "d")) +
                        P->SPH_opts.domain_min.y;
	                P->SPH_opts.domain_max.z =
                        atof(xmlAttribute(s_elem, "h")) +
                        P->SPH_opts.domain_min.z;
	            #else
	                P->SPH_opts.domain_max.y =
                        atof(xmlAttribute(s_elem, "h")) +
                        P->SPH_opts.domain_min.y;
	            #endif
	            if(!strcmp(xmlAttribute(s_elem, "move"), "true") ||
                   !strcmp(xmlAttribute(s_elem, "move"), "True")){
                    P->SPH_opts.domain_motion = true;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "move"), "false") ||
                        !strcmp(xmlAttribute(s_elem, "move"), "False")){
                    P->SPH_opts.domain_motion = false;
	            }
	            // Fix negative lengths
	            if(P->SPH_opts.domain_min.x > P->SPH_opts.domain_max.x){
	                float aux = P->SPH_opts.domain_min.x;
	                P->SPH_opts.domain_min.x = P->SPH_opts.domain_max.x;
	                P->SPH_opts.domain_max.x = aux;
	            }
	            if(P->SPH_opts.domain_min.y > P->SPH_opts.domain_max.y){
	                float aux = P->SPH_opts.domain_min.y;
	                P->SPH_opts.domain_min.y = P->SPH_opts.domain_max.y;
	                P->SPH_opts.domain_max.y = aux;
	            }
	            #ifdef HAVE_3D
	                if(P->SPH_opts.domain_min.z > P->SPH_opts.domain_max.z){
	                    float aux = P->SPH_opts.domain_min.z;
	                    P->SPH_opts.domain_min.z = P->SPH_opts.domain_max.z;
	                    P->SPH_opts.domain_max.z = aux;
	                }
	            #endif
	            // Look for errors
	            if(    (P->SPH_opts.domain_min.x == P->SPH_opts.domain_max.x)
	                || (P->SPH_opts.domain_min.y == P->SPH_opts.domain_max.y)
	                #ifdef HAVE_3D
	                    || (P->SPH_opts.domain_min.z == P->SPH_opts.domain_max.z)
	                #endif
	              )
	            {
	                sprintf(msg, "Null domain volume.\n");
	                S->addMessageF(3, msg);
	                return true;
	            }
			}
			else {
	            sprintf(msg, "Unknow SPH option\n\t\t%s\n", name);
	            S->addMessageF(3, msg);
	            return true;
			}
	    }
	}
	return false;
}

bool State::parseFluid(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Fluid"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        if(!xmlHasAttribute(elem, "n")){
            S->addMessageF(3, "Found a fluid without the \"n\" attribute.\n");
            return true;
        }

	    P->addFluid();
	    P->fluids[P->n_fluids-1].n = atoi(xmlAttribute(elem, "n"));

	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Option"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);

	        const char *name = xmlAttribute(s_elem, "name");
	        const char *value = xmlAttribute(s_elem, "value");

			if(!strcmp(name, "gamma"))
				P->fluids[P->n_fluids-1].gamma = atof(value);
			else if(!strcmp(name, "refd"))
				P->fluids[P->n_fluids-1].refd = atof(value);
			else if(!strcmp(name, "Viscdyn"))
				P->fluids[P->n_fluids-1].visc_dyn = atof(value);
			else if(!strcmp(name, "alpha"))
				P->fluids[P->n_fluids-1].alpha = atof(value);
			else if(!strcmp(name, "delta"))
				P->fluids[P->n_fluids-1].delta = atof(value);
			else{
	            sprintf(msg, "Unknow Fluid option \"%s\"\n", name);
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tgamma\n");
	            S->addMessage(0, "\t\trefd\n");
	            S->addMessage(0, "\t\tViscdyn\n");
	            S->addMessage(0, "\t\talpha\n");
	            S->addMessage(0, "\t\tdelta\n");
	            return true;
			}
	    }

	    s_nodes = elem->getElementsByTagName(xmlS("Load"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			strcpy(P->fluids[P->n_fluids-1].in_path,
                   xmlAttribute(s_elem, "file"));
            if(xmlHasAttribute(s_elem, "format")){
                strcpy(P->fluids[P->n_fluids-1].in_format,
                       xmlAttribute(s_elem, "format"));
            }
	    }

	    s_nodes = elem->getElementsByTagName(xmlS("Save"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			strcpy(P->fluids[P->n_fluids-1].out_path,
                   xmlAttribute(s_elem, "file"));
            if(xmlHasAttribute(s_elem, "format")){
                strcpy(P->fluids[P->n_fluids-1].out_format,
                       xmlAttribute(s_elem, "format"));
            }
	    }
	}
	return false;
}

bool State::parseSensors(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Sensors"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("FPS"));
	    bool haveFPS = false;
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        haveFPS = true;
	        P->SensorsParameters.fps = atof(xmlAttribute(s_elem, "value"));
	    }
		if(!haveFPS){
	        sprintf(msg, "Sensors section without FPS output frecuency.\n");
	        S->addMessageF(2, msg);
		}

	    s_nodes = elem->getElementsByTagName(xmlS("Script"));
	    bool haveScript = false;
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        haveScript = true;
	        strcpy(P->SensorsParameters.script, xmlAttribute(s_elem, "file"));
	    }
		if(!haveScript){
	        sprintf(msg, "Sensors section without script.\n");
	        S->addMessageF(2, msg);
		}

	    s_nodes = elem->getElementsByTagName(xmlS("Sensor"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
		    vec pos;
	        pos.x = atof(xmlAttribute(s_elem, "x"));
	        pos.y = atof(xmlAttribute(s_elem, "y"));
	        #ifdef HAVE_3D
	            pos.z = atof(xmlAttribute(s_elem, "z"));
	            pos.w = 0.f;
	        #endif
	        // Sensor type
	        cl_ushort mode=0;
	        if(!strcmp(xmlAttribute(s_elem, "type"), "Interpolated"))
	            mode = 0;
	        else if(!strcmp(xmlAttribute(s_elem, "type"), "Maximum"))
	            mode = 1;
	        else{
	            sprintf(msg,
                        "%s is not a valid type of sensor.\n",
                        xmlAttribute(s_elem, "type"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tInterpolated\n");
	            S->addMessage(0, "\t\tMaximum\n");
	            return true;
	        }
	        if(P->SensorsParameters.add(pos, mode))
	            return true;
	    }
	}
	return false;
}

bool State::parseMotions(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Movements"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
	    // Look for movements
	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Movement"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        // Build new movement
	        ProblemSetup::sphMoveParameters *Move =
                new ProblemSetup::sphMoveParameters();
	        // Get movement type
	        if(!strcmp(xmlAttribute(s_elem, "type"),"Quaternion")){
	            Move->type = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"LIQuaternion")){
	            Move->type = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"C1Quaternion")){
	            Move->type = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"ScriptQuaternion")){
	            Move->type = 3;
	        }
	        else{
	            sprintf(msg,
                        "Unknow type of movement \"%s\"\n",
                        xmlAttribute(s_elem, "type"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tLIQuaternion\n");
	            S->addMessage(0, "\t\tScriptQuaternion\n");
	            return true;
	        }

	        strcpy(Move->path, xmlAttribute(s_elem, "file"));

	        P->motions.push_back(Move);
	    }
	}
	return false;
}

bool State::parsePortals(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	bool hasNormal;
	ProblemSetup::sphPortal::Portal in, out;
	#ifdef HAVE_3D
	    in.corner.w = 0.f;
	    in.up.w     = 0.f;
	    in.side.w   = 0.f;
	    in.normal.w = 0.f;
	#endif
	DOMNode *ss_node;
	DOMNodeList *s_nodes;
	DOMElement *ss_elem, *s_elem;
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Portal"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
	    // Look for in/out planes
	    DOMNodeList* ss_nodes = elem->getElementsByTagName(xmlS("Plane"));
	    if(!ss_nodes->getLength()){
	        continue;
	    }
	    if(ss_nodes->getLength() != 2){
	        sprintf(msg, "Invalid planes data.\n");
	        S->addMessageF(3, msg);
	        sprintf(
                msg,
                "\t2 planes by portal must be stablished, but %u found.\n",
                ss_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }

	    // -----------------------------------------------------
	    // In portal
	    // -----------------------------------------------------
	    hasNormal = false;
	    ss_node = ss_nodes->item(0);
	    ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Corner"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg,
                    "\t1 corner must be provided, but %u found.\n",
                    s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    in.corner.x = atof(xmlAttribute(s_elem, "x"));
	    in.corner.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        in.corner.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Up"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg,
                    "\t1 up vector must be provided, but %u found.\n",
                    s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    in.up.x = atof(xmlAttribute(s_elem, "x"));
	    in.up.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        in.up.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    #ifdef HAVE_3D
	        s_nodes = ss_elem->getElementsByTagName(xmlS("Side"));
	        if(s_nodes->getLength() != 1){
	            sprintf(msg, "Invalid first plane data.\n");
	            S->addMessageF(3, msg);
	            sprintf(msg,
                        "\t1 side vector must be provided, but %u found.\n",
                        s_nodes->getLength());
	            S->addMessage(0, msg);
	            continue;
	        }
	        s_elem = dynamic_cast<xercesc::DOMElement*>(s_nodes->item(0));
	        in.side.x = atof(xmlAttribute(s_elem, "x"));
	        in.side.y = atof(xmlAttribute(s_elem, "y"));
	        in.side.z = atof(xmlAttribute(s_elem, "z"));
	    #else
	        in.side.x = 0.f;
	        in.side.y = 0.f;
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Normal"));
	    if(s_nodes->getLength() > 1){
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg,
                    "\tJust 1 normal vector can be provided, but %u found.\n",
                    s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    else if (s_nodes->getLength() == 1){
	        s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	        in.normal.x = atof(xmlAttribute(s_elem, "x"));
	        in.normal.y = atof(xmlAttribute(s_elem, "y"));
	        #ifdef HAVE_3D
	            in.normal.z = atof(xmlAttribute(s_elem, "z"));
	        #endif
	        hasNormal = true;
	    }
	    if(!hasNormal){
	        #ifdef HAVE_3D
	            in.normal = normalize(cross(in.up, in.side));
	        #else
	            in.normal.x = -in.side.y;
	            in.normal.y =  in.side.x;
	            in.normal   =  normalize(in.normal);
	        #endif
	    }
	    // -----------------------------------------------------
	    // Out portal
	    // -----------------------------------------------------
	    hasNormal = false;
	    ss_node = ss_nodes->item(1);
	    ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Corner"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid second plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg,
                    "\t1 corner must be provided, but %u found.\n",
                    s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    out.corner.x = atof(xmlAttribute(s_elem, "x"));
	    out.corner.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        out.corner.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Up"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid second plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg,
                    "\t1 up vector must be provided, but %u found.\n",
                    s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    out.up.x = atof(xmlAttribute(s_elem, "x"));
	    out.up.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        out.up.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    #ifdef HAVE_3D
	        s_nodes = ss_elem->getElementsByTagName(xmlS("Side"));
	        if(s_nodes->getLength() != 1){
                sprintf(msg, "Invalid second plane data.\n");
                S->addMessageF(3, msg);
                sprintf(msg,
                        "\t1 side vector must be provided, but %u found.\n",
                        s_nodes->getLength());
                S->addMessage(0, msg);
	            continue;
	        }
	        s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	        out.side.x = atof(xmlAttribute(s_elem, "x"));
	        out.side.y = atof(xmlAttribute(s_elem, "y"));
	        out.side.z = atof(xmlAttribute(s_elem, "z"));
	    #else
	        out.side.x = 0.f;
	        out.side.y = 0.f;
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(xmlS("Normal"));
	    if(s_nodes->getLength() > 1){
            sprintf(msg,
                    "Invalid second plane data.\n");
            S->addMessageF(3, msg);
            sprintf(msg,
                    "\tJust 1 normal vector can be provided, but %u found.\n",
                    s_nodes->getLength());
            S->addMessage(0, msg);
	        continue;
	    }
	    else if (s_nodes->getLength() == 1){
	        s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	        out.normal.x = atof(xmlAttribute(s_elem, "x"));
	        out.normal.y = atof(xmlAttribute(s_elem, "y"));
	        #ifdef HAVE_3D
	            out.normal.z = atof(xmlAttribute(s_elem, "z"));
	        #endif
	        hasNormal = true;
	    }
	    if(!hasNormal){
	        #ifdef HAVE_3D
	            out.normal = normalize(cross(out.up, out.side));
	        #else
	            out.normal.x = -out.side.y;
	            out.normal.y =  out.side.x;
	            out.normal   =  normalize(out.normal);
	        #endif
	    }

	    ProblemSetup::sphPortal *Portal = new ProblemSetup::sphPortal();
	    Portal->in  = in;
	    Portal->out = out;
	    P->portals.push_back(Portal);
	}
	return false;
}

bool State::parseGhostParticles(DOMElement *root)
{
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("GhostParticles"));
	for(XMLSize_t i=0; i<nodes->getLength(); i++){
	    DOMNode* node = nodes->item(i);
	    if(node->getNodeType() != DOMNode::ELEMENT_NODE)
	        continue;
	    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

	    DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("PressModel"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->ghost_particles.p_extension = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->ghost_particles.p_extension = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->ghost_particles.p_extension = 2;
	        }
	        else{
	            sprintf(msg,
                        "Invalid %s pressure extension model.\n",
                        xmlAttribute(s_elem, "value"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tASM\n");
	            S->addMessage(0, "\t\tSSM\n");
	            S->addMessage(0, "\t\tTakeda\n");
	            return true;
	        }
	    }
	    // Get normal velocity model
	    s_nodes = elem->getElementsByTagName(xmlS("NormalUModel"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->ghost_particles.vn_extension = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->ghost_particles.vn_extension = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->ghost_particles.vn_extension = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "U0M")){
	            P->ghost_particles.vn_extension = 3;
	        }
	        else{
	            sprintf(msg,
                        "Invalid %s normal velocity extension model.\n",
                        xmlAttribute(s_elem, "value"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tASM\n");
	            S->addMessage(0, "\t\tSSM\n");
	            S->addMessage(0, "\t\tTakeda\n");
	            S->addMessage(0, "\t\tU0M\n");
	            return true;
	        }
	    }
	    // Get normal velocity model
	    s_nodes = elem->getElementsByTagName(xmlS("TangentUModel"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->ghost_particles.vt_extension = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->ghost_particles.vt_extension = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->ghost_particles.vt_extension = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "U0M")){
	            P->ghost_particles.vt_extension = 3;
	        }
	        else{
	            sprintf(msg,
                        "Invalid %s tangent velocity extension model.\n",
                        xmlAttribute(s_elem, "value"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tASM\n");
	            S->addMessage(0, "\t\tSSM\n");
	            S->addMessage(0, "\t\tTakeda\n");
	            S->addMessage(0, "\t\tU0M\n");
	            return true;
	        }
	    }

	    s_nodes = elem->getElementsByTagName(xmlS("Wall"));
	    for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
	        DOMNode* s_node = s_nodes->item(j);
	        if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
	            continue;
	        DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
	        DOMNodeList* ss_nodes = s_elem->getElementsByTagName(
                xmlS("Vertex"));
	        #ifdef HAVE_3D
	            if((ss_nodes->getLength() != 3) &&
                   (ss_nodes->getLength() != 4)){
	                sprintf(msg,
                            "Found a wall with %u vertexes.\n",
                            ss_nodes->getLength());
	                S->addMessageF(3, msg);
	                S->addMessage(
                        0,
                        "\t3D simulations accepts walls of 3 or 4 vertexes\n"
                    );
	                return true;
	            }
	            if( ss_nodes->getLength() == 3 ){
	                vec p[3], v[3];
	                for(XMLSize_t k=0; k<ss_nodes->getLength(); k++){
	                    DOMNode* ss_node = ss_nodes->item(k);
	                    if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                        continue;
	                    DOMElement* ss_elem = dynamic_cast<xercesc::DOMElement*>(ss_node);
	                    p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                    p[k].y = atof(xmlAttribute(ss_elem, "y"));
	                    p[k].z = atof(xmlAttribute(ss_elem, "z"));
	                    p[k].w = 0.f;
	                    v[k].x = atof(xmlAttribute(ss_elem, "vx"));
	                    v[k].y = atof(xmlAttribute(ss_elem, "vy"));
	                    v[k].z = atof(xmlAttribute(ss_elem, "vz"));
	                    v[k].w = 0.f;
	                }
	                P->ghost_particles.add(p[0], p[1], p[2]
                                           v[0], v[1], v[2]);
	            }
	            if( ss_nodes->getLength() == 4 ){
	                vec p[4], v[4];
	                for(XMLSize_t k=0; k<ss_nodes->getLength(); k++){
	                    DOMNode* ss_node = ss_nodes->item(k);
	                    if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                        continue;
	                    DOMElement* ss_elem = dynamic_cast<xercesc::DOMElement*>(ss_node);
	                    p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                    p[k].y = atof(xmlAttribute(ss_elem, "y"));
	                    p[k].z = atof(xmlAttribute(ss_elem, "z"));
	                    p[k].w = 0.f;
	                    v[k].x = atof(xmlAttribute(ss_elem, "vx"));
	                    v[k].y = atof(xmlAttribute(ss_elem, "vy"));
	                    v[k].z = atof(xmlAttribute(ss_elem, "vz"));
	                    v[k].w = 0.f;
	                }
	                P->ghost_particles.add(p[0], p[1], p[2], p[3]
                                           v[0], v[1], v[2], v[3]);
	            }
	        #else
	            if(ss_nodes->getLength() != 2){
	                sprintf(msg,
                            "Found wall with %u vertexes.\n",
                            ss_nodes->getLength());
	                S->addMessageF(3, msg);
	                S->addMessage(
                        0,
                        "\t2D simulations accepts walls of 2 vertexes\n");
	                return true;
	            }
	            vec p[2], v[2];
	            for(XMLSize_t k=0; k<ss_nodes->getLength(); k++){
	                DOMNode* ss_node = ss_nodes->item(k);
	                if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                    continue;
	                DOMElement* ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	                p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                p[k].y = atof(xmlAttribute(ss_elem, "y"));
                    v[k].x = atof(xmlAttribute(ss_elem, "vx"));
                    v[k].y = atof(xmlAttribute(ss_elem, "vy"));
	            }
	            P->ghost_particles.add(p[0], p[1], v[0], v[1]);
	        #endif
	    }
	}
	return false;
}

bool State::write(const char* filepath)
{
    DOMImplementation* impl = DOMImplementationRegistry::getDOMImplementation(
        xmlS("Range"));
    DOMDocument* doc = impl->createDocument(
        NULL,
        xmlS("sphInput"),
        NULL);
    DOMElement* root = doc->getDocumentElement();

	if(writeSettings(doc, root))
	    return true;
	if(writeOpenCL(doc, root))
	    return true;
	if(writeTiming(doc, root))
        return true;
	if(writeTiming(doc, root))
	    return true;
	if(writeSPH(doc, root))
	    return true;
	if(writeFluid(doc, root))
	    return true;
	if(writeSensors(doc, root))
	    return true;
	if(writeMotions(doc, root))
	    return true;
	if(writePortals(doc, root))
	    return true;
	if(writeGhostParticles(doc, root))
	    return true;
    return false;
}

bool State::writeSettings(xercesc::DOMDocument* doc,
                          xercesc::DOMElement *root)
{
    DOMElement *elem, *s_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();

    elem = doc->createElement(xmlS("Settings"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("Verbose"));
    sprintf(att, "%d", P->settings.verbose_level);
    s_elem->setAttribute(xmlS("level"), xmlS(att));
    elem->appendChild(s_elem);
    s_elem = doc->createElement(xmlS("Device"));
    sprintf(att, "%u", P->settings.platform_id);
    s_elem->setAttribute(xmlS("platform"), xmlS(att));
    sprintf(att, "%u", P->settings.device_id);
    s_elem->setAttribute(xmlS("device"), xmlS(att));
    switch(P->settings.device_type){
    case CL_DEVICE_TYPE_ALL:
        strcpy(att, "ALL");
        break;
    case CL_DEVICE_TYPE_CPU:
        strcpy(att, "CPU");
        break;
    case CL_DEVICE_TYPE_GPU:
        strcpy(att, "GPU");
        break;
    case CL_DEVICE_TYPE_ACCELERATOR:
        strcpy(att, "ACCELERATOR");
        break;
    case CL_DEVICE_TYPE_DEFAULT:
        strcpy(att, "DEFAULT");
        break;
    }
    s_elem->setAttribute(xmlS("device"), xmlS(att));
    elem->appendChild(s_elem);
    return false;
}

bool State::writeOpenCL(xercesc::DOMDocument* doc,
                        xercesc::DOMElement *root)
{
    DOMElement *elem, *s_elem;
	ProblemSetup *P = ProblemSetup::singleton();

    elem = doc->createElement(xmlS("OpenCL"));
    root->appendChild(elem);

    std::map<char*, char*> tags;
    std::map<char*, char*>::iterator it;
    tags["Predictor"] = P->OpenCL_kernels.predictor;
    tags["LinkList"] = P->OpenCL_kernels.link_list;
    tags["Rates"] = P->OpenCL_kernels.rates;
    tags["Corrector"] = P->OpenCL_kernels.corrector;
    tags["TimeStep"] = P->OpenCL_kernels.time_step;
    tags["Reduction"] = P->OpenCL_kernels.reduction;
    tags["RadixSort"] = P->OpenCL_kernels.radix_sort;
    tags["ElasticBounce"] = P->OpenCL_kernels.elastic_bounce;
    tags["DeLeffe"] = P->OpenCL_kernels.de_Leffe;
    tags["GhostParticles"] = P->OpenCL_kernels.ghost;
    tags["Shepard"] = P->OpenCL_kernels.shepard;
    tags["DensInterpolation"] = P->OpenCL_kernels.dens_int;
    tags["Torque"] = P->OpenCL_kernels.torque;
    tags["Energy"] = P->OpenCL_kernels.energy;
    tags["Bounds"] = P->OpenCL_kernels.bounds;
    tags["Domain"] = P->OpenCL_kernels.domain;
    tags["Portal"] = P->OpenCL_kernels.portal;

    for(it=tags.begin(); it!=tags.end(); it++){
        s_elem = doc->createElement(xmlS(it->first));
        s_elem->setAttribute(xmlS("file"), xmlS(it->second));
        elem->appendChild(s_elem);
    }
    return false;
}

bool State::writeTiming(xercesc::DOMDocument* doc,
                        xercesc::DOMElement *root)
{
    DOMElement *elem, *s_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();
	TimeManager *T = TimeManager::singleton();

    elem = doc->createElement(xmlS("Timing"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("Start"));
    sprintf(att, "%g", T->time());
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    if(P->time_opts.sim_end_mode & __TIME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Time"));
        sprintf(att, "%g", P->time_opts.sim_end_time);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.sim_end_mode & __ITER_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Steps"));
        sprintf(att, "%d", P->time_opts.sim_end_step);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.sim_end_mode & __FRAME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Frames"));
        sprintf(att, "%d", P->time_opts.sim_end_frame);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("TimeStep"));
    switch(P->time_opts.dt_mode){
    case __DT_VARIABLE__:
        strcpy(att, "Courant");
        break;
    case __DT_FIXCALCULATED__:
        strcpy(att, "Fix");
        break;
    case __DT_FIX__:
        sprintf(att, "%g", P->time_opts.dt);
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("MinTimeStep"));
    sprintf(att, "%g", P->time_opts.dt_min);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("Courant"));
    sprintf(att, "%g", P->time_opts.courant);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("ClampVel"));
    strcpy(att, "false");
    if(P->time_opts.velocity_clamp)
        strcpy(att, "true");
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    if(P->time_opts.log_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("LogFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", P->time_opts.log_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.log_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("LogFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", P->time_opts.log_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    if(P->time_opts.energy_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("EnFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", P->time_opts.energy_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.energy_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("EnFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", P->time_opts.energy_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    if(P->time_opts.bounds_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("BoundsFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", P->time_opts.bounds_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.bounds_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("BoundsFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", P->time_opts.bounds_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    if(P->time_opts.output_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", P->time_opts.output_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(P->time_opts.output_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", P->time_opts.output_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    return false;
}

bool State::writeSPH(xercesc::DOMDocument* doc,
                     xercesc::DOMElement *root)
{
    DOMElement *elem, *s_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();

    elem = doc->createElement(xmlS("SPH"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("g"));
    sprintf(att, "%g", P->SPH_opts.g.x);
    s_elem->setAttribute(xmlS("x"), xmlS(att));
    sprintf(att, "%g", P->SPH_opts.g.y);
    s_elem->setAttribute(xmlS("y"), xmlS(att));
    #ifdef HAVE_3D
        sprintf(att, "%g", P->SPH_opts.g.z);
        s_elem->setAttribute(xmlS("z"), xmlS(att));
        elem->appendChild(s_elem);
    #endif

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("hfac"));
    sprintf(att, "%g", P->SPH_opts.hfac);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("deltar"));
    sprintf(att, "%g", P->SPH_opts.deltar.x);
    s_elem->setAttribute(xmlS("x"), xmlS(att));
    sprintf(att, "%g", P->SPH_opts.deltar.y);
    s_elem->setAttribute(xmlS("y"), xmlS(att));
    #ifdef HAVE_3D
        sprintf(att, "%g", P->SPH_opts.deltar.z);
        s_elem->setAttribute(xmlS("z"), xmlS(att));
        elem->appendChild(s_elem);
    #endif // HAVE_3D

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("cs"));
    sprintf(att, "%g", P->SPH_opts.cs);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("LLSteps"));
    sprintf(att, "%u", P->SPH_opts.link_list_steps);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("DensSteps"));
    sprintf(att, "%u", P->SPH_opts.dens_int_steps);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("DensBounds"));
    sprintf(att, "%g", P->SPH_opts.rho_min);
    s_elem->setAttribute(xmlS("min"), xmlS(att));
    if(P->SPH_opts.rho_max < std::numeric_limits<float>::max()){
        sprintf(att, "%g", P->SPH_opts.rho_max);
        s_elem->setAttribute(xmlS("max"), xmlS(att));
    }
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("Boundary"));
    switch(P->SPH_opts.boundary_type){
    case 0:
        strcpy(att, "ElasticBounce");
        break;
    case 1:
        strcpy(att, "FixedParticles");
        break;
    case 2:
        strcpy(att, "DeLeffe");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("SlipCondition"));
    switch(P->SPH_opts.boundary_type){
    case 0:
        strcpy(att, "NoSlip");
        break;
    case 1:
        strcpy(att, "FreeSlip");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("BoundElasticFactor"));
    sprintf(att, "%g", P->SPH_opts.elastic_factor);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("BoundDist"));
    sprintf(att, "%g", P->SPH_opts.elastic_factor);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Option"));
    s_elem->setAttribute(xmlS("name"), xmlS("Shepard"));
    switch(P->SPH_opts.has_shepard){
    case 0:
        strcpy(att, "None");
        break;
    case 1:
        strcpy(att, "Force");
        break;
    case 2:
        strcpy(att, "Dens");
        break;
    case 3:
        strcpy(att, "ForceDens");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    if(P->SPH_opts.has_domain){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Domain"));
        strcpy(att, "false");
        if(P->SPH_opts.domain_motion)
            strcpy(att, "true");
        s_elem->setAttribute(xmlS("move"), xmlS(att));

        sprintf(att, "%g", P->SPH_opts.domain_min.x);
        s_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->SPH_opts.domain_min.y);
        s_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->SPH_opts.domain_min.z);
            s_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D
        sprintf(att,
                "%g",
                P->SPH_opts.domain_max.x - P->SPH_opts.domain_min.x);
        s_elem->setAttribute(xmlS("l"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att,
                    "%g",
                    P->SPH_opts.domain_max.y - P->SPH_opts.domain_min.y);
            s_elem->setAttribute(xmlS("d"), xmlS(att));
            sprintf(att,
                    "%g",
                    P->SPH_opts.domain_max.z - P->SPH_opts.domain_min.z);
            s_elem->setAttribute(xmlS("h"), xmlS(att));
        #else
            sprintf(att,
                    "%g",
                    P->SPH_opts.domain_max.y - P->SPH_opts.domain_min.y);
            s_elem->setAttribute(xmlS("h"), xmlS(att));
        #endif // HAVE_3D
        elem->appendChild(s_elem);
    }

    return false;
}

bool State::writeFluid(xercesc::DOMDocument* doc,
                       xercesc::DOMElement *root)
{
    unsigned int i;
    DOMElement *elem, *s_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();

    for(i=0; i<P->n_fluids; i++){
        elem = doc->createElement(xmlS("Fluid"));
        sprintf(att, "%u", P->fluids[i].n);
        elem->setAttribute(xmlS("n"), xmlS(att));
        root->appendChild(elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("gamma"));
        sprintf(att, "%g", P->fluids[i].gamma);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("refd"));
        sprintf(att, "%g", P->fluids[i].refd);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Viscdyn"));
        sprintf(att, "%g", P->fluids[i].visc_dyn);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("alpha"));
        sprintf(att, "%g", P->fluids[i].alpha);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("delta"));
        sprintf(att, "%g", P->fluids[i].delta);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("delta"));
        sprintf(att, "%g", P->fluids[i].delta);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Load"));
        /// @todo Get the loading file from the fluid saver
        s_elem->setAttribute(xmlS("file"), xmlS(""));
        s_elem->setAttribute(xmlS("format"), xmlS(""));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Save"));
        s_elem->setAttribute(xmlS("file"), xmlS(P->fluids[i].out_path));
        s_elem->setAttribute(xmlS("format"), xmlS(P->fluids[i].out_format));
        elem->appendChild(s_elem);
    }

    return false;
}

bool State::writeSensors(xercesc::DOMDocument* doc,
                         xercesc::DOMElement *root)
{
    unsigned int i;
    DOMElement *elem, *s_elem, *ss_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    elem = doc->createElement(xmlS("Sensors"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("FPS"));
    sprintf(att, "%g", P->SensorsParameters.fps);
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("Script"));
    s_elem->setAttribute(xmlS("file"), xmlS(P->SensorsParameters.script));
    elem->appendChild(s_elem);

    vec *pos = C->sensors->positions();
    for(i=0; i<P->SensorsParameters.pos.size(); i++){
        ss_elem = doc->createElement(xmlS("Sensor"));
        sprintf(att, "%g", pos[i].x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", pos[i].y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", pos[i].z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        switch(P->SensorsParameters.mod.at(i)){
        case 0:
            strcpy(att, "Interpolated");
            break;
        case 1:
            strcpy(att, "Maximum");
            break;
        }
        s_elem->setAttribute(xmlS("type"), xmlS(att));

        s_elem->appendChild(ss_elem);
    }

    return false;
}

bool State::writeMotions(xercesc::DOMDocument* doc,
                         xercesc::DOMElement *root)
{
    unsigned int i;
    DOMElement *elem, *s_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();

    elem = doc->createElement(xmlS("Movements"));
    root->appendChild(elem);

    for(i=0; i<P->motions.size(); i++){
        s_elem = doc->createElement(xmlS("Movement"));
        switch(P->motions.at(i)->type){
        case 0:
            strcpy(att, "Quaternion");
            break;
        case 1:
            strcpy(att, "LIQuaternion");
            break;
        case 2:
            strcpy(att, "C1Quaternion");
            break;
        case 3:
            strcpy(att, "ScriptQuaternion");
            break;
        }
        s_elem->setAttribute(xmlS("type"), xmlS(att));
        s_elem->setAttribute(xmlS("file"), xmlS(P->motions.at(i)->path));

        elem->appendChild(s_elem);
    }

    return false;
}

bool State::writePortals(xercesc::DOMDocument* doc,
                         xercesc::DOMElement *root)
{
    unsigned int i;
    DOMElement *elem, *s_elem, *ss_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();

    elem = doc->createElement(xmlS("Portal"));
    root->appendChild(elem);

    for(i=0; i<P->portals.size(); i++){
        s_elem = doc->createElement(xmlS("Plane"));
        elem->appendChild(s_elem);

        ss_elem = doc->createElement(xmlS("Corner"));
        sprintf(att, "%g", P->portals.at(i)->in.corner.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->in.corner.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->in.corner.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);

        #ifdef HAVE_3D
            ss_elem = doc->createElement(xmlS("Side"));
            sprintf(att, "%g", P->portals.at(i)->in.side.x);
            ss_elem->setAttribute(xmlS("x"), xmlS(att));
            sprintf(att, "%g", P->portals.at(i)->in.side.y);
            ss_elem->setAttribute(xmlS("y"), xmlS(att));
            sprintf(att, "%g", P->portals.at(i)->in.side.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));

            s_elem->appendChild(ss_elem);
        #endif // HAVE_3D

        ss_elem = doc->createElement(xmlS("Up"));
        sprintf(att, "%g", P->portals.at(i)->in.up.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->in.up.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->in.up.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);

        ss_elem = doc->createElement(xmlS("Normal"));
        sprintf(att, "%g", P->portals.at(i)->in.normal.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->in.normal.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->in.normal.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);

        s_elem = doc->createElement(xmlS("Plane"));
        elem->appendChild(s_elem);

        ss_elem = doc->createElement(xmlS("Corner"));
        sprintf(att, "%g", P->portals.at(i)->out.corner.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->out.corner.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->out.corner.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);

        #ifdef HAVE_3D
            ss_elem = doc->createElement(xmlS("Side"));
            sprintf(att, "%g", P->portals.at(i)->out.side.x);
            ss_elem->setAttribute(xmlS("x"), xmlS(att));
            sprintf(att, "%g", P->portals.at(i)->out.side.y);
            ss_elem->setAttribute(xmlS("y"), xmlS(att));
            sprintf(att, "%g", P->portals.at(i)->out.side.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));

            s_elem->appendChild(ss_elem);
        #endif // HAVE_3D

        ss_elem = doc->createElement(xmlS("Up"));
        sprintf(att, "%g", P->portals.at(i)->out.up.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->out.up.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->out.up.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);

        ss_elem = doc->createElement(xmlS("Normal"));
        sprintf(att, "%g", P->portals.at(i)->out.normal.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->portals.at(i)->out.normal.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->portals.at(i)->out.normal.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D

        s_elem->appendChild(ss_elem);
    }

    return false;
}

bool State::writeGhostParticles(xercesc::DOMDocument* doc,
                                xercesc::DOMElement *root)
{
    unsigned int i;
    DOMElement *elem, *s_elem, *ss_elem;
    char att[1024];
	ProblemSetup *P = ProblemSetup::singleton();
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    elem = doc->createElement(xmlS("GhostParticles"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("PressModel"));
    switch(P->ghost_particles.p_extension){
    case 0:
        strcpy(att, "ASM");
        break;
    case 1:
        strcpy(att, "SSM");
        break;
    case 2:
        strcpy(att, "Takeda");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("NormalUModel"));
    switch(P->ghost_particles.vn_extension){
    case 0:
        strcpy(att, "ASM");
        break;
    case 1:
        strcpy(att, "SSM");
        break;
    case 2:
        strcpy(att, "Takeda");
        break;
    case 3:
        strcpy(att, "U0M");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    s_elem = doc->createElement(xmlS("TangentUModel"));
    switch(P->ghost_particles.vt_extension){
    case 0:
        strcpy(att, "ASM");
        break;
    case 1:
        strcpy(att, "SSM");
        break;
    case 2:
        strcpy(att, "Takeda");
        break;
    case 3:
        strcpy(att, "U0M");
        break;
    }
    s_elem->setAttribute(xmlS("value"), xmlS(att));
    elem->appendChild(s_elem);

    for(i=0; i<P->ghost_particles.walls.size(); i++){
        ss_elem = doc->createElement(xmlS("Vertex"));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->p1.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->p1.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p1.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->v1.x);
        ss_elem->setAttribute(xmlS("vx"), xmlS(att));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->v1.y);
        ss_elem->setAttribute(xmlS("vy"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v1.z);
            ss_elem->setAttribute(xmlS("vz"), xmlS(att));
        #endif // HAVE_3D
        s_elem->appendChild(ss_elem);

        ss_elem = doc->createElement(xmlS("Vertex"));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->p2.x);
        ss_elem->setAttribute(xmlS("x"), xmlS(att));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->p2.y);
        ss_elem->setAttribute(xmlS("y"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p2.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
        #endif // HAVE_3D
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->v2.x);
        ss_elem->setAttribute(xmlS("vx"), xmlS(att));
        sprintf(att, "%g", P->ghost_particles.walls.at(i)->v2.y);
        ss_elem->setAttribute(xmlS("vy"), xmlS(att));
        #ifdef HAVE_3D
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v2.z);
            ss_elem->setAttribute(xmlS("vz"), xmlS(att));
        #endif // HAVE_3D
        s_elem->appendChild(ss_elem);

        #ifdef HAVE_3D
            ss_elem = doc->createElement(xmlS("Vertex"));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p3.x);
            ss_elem->setAttribute(xmlS("x"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p3.y);
            ss_elem->setAttribute(xmlS("y"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p3.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v3.x);
            ss_elem->setAttribute(xmlS("vx"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v3.y);
            ss_elem->setAttribute(xmlS("vy"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v3.z);
            ss_elem->setAttribute(xmlS("vz"), xmlS(att));
            s_elem->appendChild(ss_elem);

            if((p3.x == p4.x) && (p3.y == p4.y) && (p3.z == p4.z))
                continue;

            ss_elem = doc->createElement(xmlS("Vertex"));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p4.x);
            ss_elem->setAttribute(xmlS("x"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p4.y);
            ss_elem->setAttribute(xmlS("y"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->p4.z);
            ss_elem->setAttribute(xmlS("z"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v4.x);
            ss_elem->setAttribute(xmlS("vx"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v4.y);
            ss_elem->setAttribute(xmlS("vy"), xmlS(att));
            sprintf(att, "%g", P->ghost_particles.walls.at(i)->v4.z);
            ss_elem->setAttribute(xmlS("vz"), xmlS(att));
            s_elem->appendChild(ss_elem);
        #endif // HAVE_3D
    }

    return false;
}

}}  // namespace
