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

#include <FileManager.h>
#include <ArgumentsManager.h>
#include <Fluid.h>
#include <CalcServer.h>
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// xerces XML parser stuff
// ----------------------------------------------------------------------------
#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )

#ifdef xmlHasAttribute
	#undef xmlHasAttribute
#endif
#define xmlHasAttribute(elem, att) elem->hasAttribute(XMLString::transcode(att))
using namespace xercesc;
using namespace std;

namespace Aqua{ namespace InputOutput{

FileManager::FileManager()
	: _in_file(0)
	, _out_prefix(0)
	, _n_out(0)
	, _log_file(0)
	, _en_file(0)
	, _bounds_file(0)
	, _log_file_id(0)
	, _en_file_id(0)
	, _bounds_file_id(0)
	#ifdef HAVE_H5PART
	    , _H5Part_file_id(0)
	#endif
	, _must_print_H5(false)
	, _must_print_VTK(false)
	, _must_print_Tecplot(false)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	_in_file     = new char[256];
	_log_file    = new char[256];
	_en_file     = new char[256];
	_bounds_file = new char[256];
	_out_prefix  = new char[64];
	// Default names for files
	strcpy(_in_file, "Input.xml");
	strcpy(_log_file, "Log.html");
	strcpy(_en_file, "Energy.dat");
	strcpy(_bounds_file, "Bounds.dat");
	strcpy(_out_prefix, "");
	// Start XML parser
	try {
	    XMLPlatformUtils::Initialize();
	}
	catch( XMLException& e ) {
	    char* message = XMLString::transcode( e.getMessage() );
        sprintf(msg, "XML toolkit initialization error.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
	    XMLString::release( &message );
	    exit(1);
	}
	catch( ... ) {
        sprintf(msg, "XML toolkit initialization error.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\tUnhandled exception\n");
        S->addMessage(0, msg);
	    exit(2);
	}
}

FileManager::~FileManager()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	delete[] _in_file; _in_file=0;
	delete[] _log_file; _log_file=0;
	delete[] _en_file; _en_file=0;
	delete[] _bounds_file; _bounds_file=0;
	delete[] _out_prefix; _out_prefix=0;
	for(i=0;i<_excluded_fields.size();i++){
	    delete[] _excluded_fields.at(i);
	}
	_excluded_fields.clear();
	// Terminate Xerces
	try {
	    XMLPlatformUtils::Terminate();
	}
	catch( xercesc::XMLException& e ) {
	    char* message = xercesc::XMLString::transcode( e.getMessage() );
        sprintf(msg, "XML toolkit exit error.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
	    XMLString::release( &message );
	}
	catch( ... ) {
        sprintf(msg, "XML toolkit exit error.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\tUnhandled exception\n");
        S->addMessage(0, msg);
	}
}

bool FileManager::parse()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024];
	strcpy(msg, "");
	// Open the xml file
	sprintf(msg, "Parsing the XML file \"%s\"\n", _in_file);
	S->addMessageF(1, msg);
	// Try to open as ascii file, just to know if the file already exist
	FILE *dummy=0; dummy = fopen(_in_file, "r");
	if(!dummy){
	    sprintf(msg, "The file doesn't exist.\n\t\t\"%s\"\n", _in_file);
	    S->addMessageF(3, msg);
	    return true;
	}
	fclose(dummy);
	XercesDOMParser *parser = new XercesDOMParser;
	parser->setValidationScheme( XercesDOMParser::Val_Never );
	parser->setDoNamespaces( false );
	parser->setDoSchema( false );
	parser->setLoadExternalDTD( false );
	parser->parse( _in_file );
 	DOMDocument* doc = parser->getDocument();
	DOMElement* root = doc->getDocumentElement();
	if( !root ){
	    sprintf(msg, "Empty XML file\n");
	    S->addMessageF(3, msg);
	    return true;
	}
	// Look for XML files included to parse them before
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Include"));
	char *orig_in_file = new char[256];
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
		// Backup the original input file
		strcpy(orig_in_file, _in_file);
		// Gets the included file
		strcpy(_in_file, XMLString::transcode( elem->getAttribute(XMLString::transcode("file")) ));
		// Calls recursively to parse it
		if(parse())
			return true;
		// Restore the backuped file
		strcpy(_in_file, orig_in_file);
	}
	delete[] orig_in_file; orig_in_file=0;
	// Parse
	parseTiming(root);
	int format = P->time_opts.output_format;
	if(format & __H5Part__) {
		_must_print_H5 = true;
	}
	if(format & __VTK__) {
		_must_print_VTK = true;
	}
	if(format & __TECPLOT__) {
		_must_print_Tecplot = true;
	}
	if(parseSettings(root))
	    return true;
	if(parseOpenCL(root))
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

bool FileManager::startLog()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	sprintf(_log_file, "%s%s", _out_prefix, _log_file);
	_log_file_id = fopen( _log_file, "w" );
	if(!_log_file_id)
	{
	    sprintf(msg, "The log file \"%s\" could not be writen\n",_log_file);
	    S->addMessageF(3, msg);
	    return true;
	}
	// Write the initial data (web header)
	initLogFile();
	return false;
}

bool FileManager::openFiles()
{
	ProblemSetup *P = P->singleton();
	ArgumentsManager *A = A->singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	sprintf(_log_file, "%s%s", _out_prefix, _log_file);
	sprintf(_en_file, "%s%s", _out_prefix, _en_file);
	sprintf(_bounds_file, "%s%s", _out_prefix, _bounds_file);
	// Open report files
	if(P->time_opts.energy_mode > __NO_OUTPUT_MODE__)
	{
	    _en_file_id = fopen( _en_file, "w" );
		if(!_en_file_id)
		{
            sprintf(msg, "The energy file \"%s\" could not be writen\n",_en_file);
            S->addMessageF(3, msg);
		    return true;
		}
	}
	if(P->time_opts.bounds_mode > __NO_OUTPUT_MODE__)
	{
	    _bounds_file_id = fopen( _bounds_file, "w" );
		if(!_bounds_file_id)
		{
            sprintf(msg, "The bounds report file \"%s\" could not be writen\n",_bounds_file);
            S->addMessageF(3, msg);
		    return true;
		}
	}
	// Open the main output files
	#ifdef HAVE_H5PART
	    #ifdef HAVE_MPI
	        // Init MPI for H5Part
	        int argc = A->argc();
	        char** argv = A->argv();
	        MPI_Init(&argc, &argv);
	    #endif
	    // Look for already existing files
	    char file_name[256];
	    if(P->settings.start_mode == 1){
	        // Look for the first empty file.
	        int index = 0;
	        char file_test[256];
	        sprintf(file_test, "%sParticles%d.h5part", _out_prefix, index);
	        while(isFile(file_test)){
	            index++;
	            sprintf(file_test, "%sParticles%d.h5part", _out_prefix, index);
	        }
	        if(!index){
                sprintf(msg, "No valid files found.\n");
                S->addMessageF(2, msg);
                sprintf(msg, "\tThe simulation will start at \"t = 0 s\" therefore.\n");
                S->addMessage(0, msg);
	            P->settings.start_mode = 0;
	        }
	        sprintf(file_name, "%sParticles%d.h5part", _out_prefix, index);
	        sprintf(P->settings.fluid_file, "%sParticles%d.h5part", _out_prefix, index-1);
	    }
	    else{
	        sprintf(file_name, "%sParticles0.h5part", _out_prefix);
	    }
	    // Open the output files
	    if(P->time_opts.output_mode > __NO_OUTPUT_MODE__)
	    {
	        if(P->time_opts.output_format & __H5Part__)
	        {
                sprintf(msg, "Saving into \"%s\"\n", file_name);
                S->addMessageF(1, msg);
	            #ifdef HAVE_MPI
	                _H5Part_file_id = H5PartOpenFileParallel(file_name,H5PART_WRITE,MPI_COMM_WORLD);
	            #else
	                _H5Part_file_id = H5PartOpenFile(file_name,H5PART_WRITE);
	            #endif
	            H5PartWriteFileAttribString(_H5Part_file_id, "rmsUnit", "m");
	        }
	    }
	#endif
	// Write initial data into the files (headers)
	initEnFile();
	initBoundsFile();
	return false;
}

void FileManager::excludeField(const char* field)
{
	char* str = NULL;
	str       = (char*)malloc((strlen(field)+1)*sizeof(char));
	strcpy(str, field);
	_excluded_fields.push_back(str);
}

bool FileManager::isField(const char* field)
{
	unsigned int i;
	for(i=0;i<_excluded_fields.size();i++){
	    if(!strcmp(_excluded_fields.at(i),field)){
	        return false;
	    }
	}
	return true;
}

void FileManager::initLogFile()
{
	if(!_log_file_id)
	    return;
	fprintf(_log_file_id, "<html>\n");
	fprintf(_log_file_id, "<head><title>AQUAgpusph log file.</title></head>\n");
	fprintf(_log_file_id, "<body bgcolor=\"#f0ffff\">\n");
	fprintf(_log_file_id, "<h1 align=\"center\">AQUAgpusph log file.</h1>\n");
	// Starting data
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	fprintf(_log_file_id, "<p align=\"left\">%s</p>\n", ctime(&seconds));
	fprintf(_log_file_id, "<hr><br>\n");
	fflush(_log_file_id);
}

void FileManager::endLogFile()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	// Compute fluid mass
	float fluidMass = 0.f;
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	InputOutput::Fluid *F = InputOutput::Fluid::singleton();
	int *imove = new int[C->N];
	float *mass = new float[C->N];
	int clFlag=0;
	clFlag |= C->getData(imove, C->imove, C->N * sizeof( int ));
	clFlag |= C->getData(mass, C->mass, C->N * sizeof( float ));
	for(i=0;i<C->N;i++) {
	    if(imove[i] > 0)
	        fluidMass += mass[i];
	}
	delete[] imove; imove=0;
	delete[] mass; mass=0;
    sprintf(msg, "Lost mass = %g [kg] (from %g [kg])\n", C->fluidMass - fluidMass, C->fluidMass);
    S->addMessageF(1, msg);
	// Closing data data
	if(!_log_file_id)
	    return;
	fprintf(_log_file_id, "<br><hr>\n");
	struct timeval now_time;
	gettimeofday(&now_time, NULL);
	const time_t seconds = now_time.tv_sec;
	fprintf(_log_file_id, "<b><font color=\"#000000\">[INFO (End of simulation)] Lost mass = %g [kg] (from %g [kg])</font></b><br>\n", C->fluidMass - fluidMass, C->fluidMass);
	fprintf(_log_file_id, "<p align=\"left\">%s</p>\n", ctime(&seconds));

	fprintf(_log_file_id, "</body>\n");
	fprintf(_log_file_id, "</html>\n");
}

void FileManager::initEnFile()
{
	if(!_en_file_id)
	    return;

	fprintf(_en_file_id,"#########################################################\n");
	fprintf(_en_file_id,"#                                                       #\n");
	fprintf(_en_file_id,"#    #    ##   #  #   #                           #     #\n");
	fprintf(_en_file_id,"#   # #  #  #  #  #  # #                          #     #\n");
	fprintf(_en_file_id,"#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
	fprintf(_en_file_id,"#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
	fprintf(_en_file_id,"#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
	fprintf(_en_file_id,"#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
	fprintf(_en_file_id,"#                            # #             #          #\n");
	fprintf(_en_file_id,"#                          ##  #             #          #\n");
	fprintf(_en_file_id,"#                                                       #\n");
	fprintf(_en_file_id,"#########################################################\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"# Another QUAlity GPU-SPH, by CEHINAV.\n");
	fprintf(_en_file_id,"#\thttp://canal.etsin.upm.es/\n");
	fprintf(_en_file_id,"# Authors:\n");
	fprintf(_en_file_id,"#\tCercós Pita, Jose Luis\n");
	fprintf(_en_file_id,"#\tMiguel Gonzalez, Leo\n");
	fprintf(_en_file_id,"#\tSouto Iglesias, Antonio\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"#########################################################\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"# Energy report file autogenerated during the simulation.\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"# E.Pot. = Potential energy (almost unuseful).\n");
	fprintf(_en_file_id,"# E.Kin. = Kinetic energy.\n");
	fprintf(_en_file_id,"# U = Internal energy (H + TS).\n");
	fprintf(_en_file_id,"# H = Enthalpy.\n");
	fprintf(_en_file_id,"# TS = Entropy (multiplied by the temperature).\n");
	fprintf(_en_file_id,"# E = Energy (U + E.kin.).\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"# The entropy generated by the viscous terms has not been\n");
	fprintf(_en_file_id,"# implemented yet\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"# Columns:\n");
	fprintf(_en_file_id,"# Time, E.Pot., E.Kin., U, H, TS, E\n");
	fprintf(_en_file_id,"#\n");
	fprintf(_en_file_id,"#########################################################\n");
	fprintf(_en_file_id,"#\n");
	fflush(_en_file_id);
}

void FileManager::endEnFile()
{
	if(!_en_file_id)
	    return;
	fprintf(_en_file_id,"# Simulation finished");
}

void FileManager::initBoundsFile()
{
	if(!_bounds_file_id)
	    return;

	fprintf(_bounds_file_id,"#########################################################\n");
	fprintf(_bounds_file_id,"#                                                       #\n");
	fprintf(_bounds_file_id,"#    #    ##   #  #   #                           #     #\n");
	fprintf(_bounds_file_id,"#   # #  #  #  #  #  # #                          #     #\n");
	fprintf(_bounds_file_id,"#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
	fprintf(_bounds_file_id,"#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
	fprintf(_bounds_file_id,"#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
	fprintf(_bounds_file_id,"#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
	fprintf(_bounds_file_id,"#                            # #             #          #\n");
	fprintf(_bounds_file_id,"#                          ##  #             #          #\n");
	fprintf(_bounds_file_id,"#                                                       #\n");
	fprintf(_bounds_file_id,"#########################################################\n");
	fprintf(_bounds_file_id,"#\n");
	fprintf(_bounds_file_id,"# Another QUAlity GPU-SPH, by CEHINAV.\n");
	fprintf(_bounds_file_id,"#\thttp://canal.etsin.upm.es/\n");
	fprintf(_bounds_file_id,"# Authors:\n");
	fprintf(_bounds_file_id,"#\tCercós Pita, Jose Luis\n");
	fprintf(_bounds_file_id,"#\tMiguel Gonzalez, Leo\n");
	fprintf(_bounds_file_id,"#\tSouto Iglesias, Antonio\n");
	fprintf(_bounds_file_id,"#\n");
	fprintf(_bounds_file_id,"#########################################################\n");
	fprintf(_bounds_file_id,"#\n");
	fprintf(_bounds_file_id,"# Fluid bounds report file.\n");
	fprintf(_bounds_file_id,"# Columns:\n");
	fprintf(_bounds_file_id,"# Time, X_Min, Y_Min, ");
	#ifdef HAVE_3D
	    fprintf(_bounds_file_id,"Z_Min, ");
	#endif
	fprintf(_bounds_file_id,"X_Max, Y_Max, ");
	#ifdef HAVE_3D
	    fprintf(_bounds_file_id,"Z_Max, ");
	#endif
	fprintf(_bounds_file_id,"V_Min, V_Max\n");
	fprintf(_bounds_file_id,"#\n");
	fprintf(_bounds_file_id,"#########################################################\n");
	fprintf(_bounds_file_id,"#\n");
	fflush(_bounds_file_id);
}

void FileManager::endBoundsFile()
{
	if(!_bounds_file_id)
	    return;
	fprintf(_bounds_file_id,"# Simulation finished");
}

bool FileManager::closeFiles()
{
	ProblemSetup *P = P->singleton();
	ArgumentsManager *A = A->singleton();
	endLogFile();
	endEnFile();
	endBoundsFile();
	if(_log_file_id)
		fclose( _log_file_id ); _log_file_id=0;
	if(_en_file_id)
		fclose( _en_file_id ); _en_file_id=0;
	if(_bounds_file_id)
		fclose( _bounds_file_id ); _bounds_file_id=0;
	#ifdef HAVE_H5PART
	    if(P->time_opts.output_mode > __NO_OUTPUT_MODE__) {
	        if(P->time_opts.output_format & __H5Part__) {
	            H5PartCloseFile(_H5Part_file_id);
	            if(A->mustReassembly()){
	                assemblyH5Part();
	            }
	        }
	    }
	    #ifdef HAVE_MPI
	        MPI_Finalize();
	    #endif
	#endif

	return false;
}

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
//			 ||									 ||
//			 ||		Auxiliar methods			 ||
//			\  /								\  /
//			 \/									 \/
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

bool FileManager::parseSettings(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Settings"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Looking for verbose
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Verbose"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        P->settings.verbose_level = atoi(xmlAttribute(s_elem, "level"));
	    }
	    // Looking for starting mode
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Start"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        P->settings.start_mode = atoi(xmlAttribute(s_elem, "mode"));
	    }
	    // Looking for OpenCL device
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Device"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
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
	            sprintf(msg, "Unknow \"%s\" type of device\n", xmlAttribute(s_elem, "type"));
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

bool FileManager::parseOpenCL(DOMElement *root)
{
	ProblemSetup *P = ProblemSetup::singleton();
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("OpenCL"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Looking for predictor
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Predictor"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.predictor, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for link-list
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("LinkList"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.link_list, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for rates
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Rates"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.rates, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for corrector
	    s_nodes = root->getElementsByTagName(XMLString::transcode("Corrector"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.corrector, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for time step
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("TimeStep"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.time_step, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for reduction
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Reduction"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.reduction, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for radix sort
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("RadixSort"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.radix_sort, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for ElasticBounce boundary condition
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("ElasticBounce"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.elastic_bounce, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for DeLeffe
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("DeLeffe"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.de_Leffe, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Ghost particles
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("GhostParticles"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.ghost, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for 0th order correction
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Shepard"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.shepard, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Density interpolation
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("DensInterpolation"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.dens_int, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Torque
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Torque"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.torque, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Energy
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Energy"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.energy, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Bounds
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Bounds"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.bounds, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Domain
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Domain"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.domain, xmlAttribute(s_elem, "file"));
	    }
	    // Looking for Portals
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Portal"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        strcpy(P->OpenCL_kernels.portal, xmlAttribute(s_elem, "file"));
	    }
	}
	return false;
}

bool FileManager::parseTiming(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Timing"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Get options
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Option"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
			if(!strcmp(xmlAttribute(s_elem, "name"), "SimulationStop")) {
				if(!strcmp(xmlAttribute(s_elem, "type"), "Time") || !strcmp(xmlAttribute(s_elem, "type"), "T")) {
					P->time_opts.sim_end_mode = P->time_opts.sim_end_mode | __TIME_MODE__;
					P->time_opts.sim_end_time = atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "Steps") || !strcmp(xmlAttribute(s_elem, "type"), "S")) {
					P->time_opts.sim_end_mode = P->time_opts.sim_end_mode | __ITER_MODE__;
					P->time_opts.sim_end_step = atoi(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "Frames") || !strcmp(xmlAttribute(s_elem, "type"), "F")) {
					P->time_opts.sim_end_mode = P->time_opts.sim_end_mode | __FRAME_MODE__;
					P->time_opts.sim_end_frame = atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg, "Unknow simulation stop criteria \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tTime\n");
	                S->addMessage(0, "\t\tSteps\n");
	                S->addMessage(0, "\t\tFrames\n");
	                return true;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "LogFiles")) {
				if(!strcmp(xmlAttribute(s_elem, "type"), "No")) {
					P->time_opts.log_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "FPS")) {
					P->time_opts.log_mode = P->time_opts.log_mode | __FPS_MODE__;
					P->time_opts.log_fps = atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "IPF")) {
					P->time_opts.log_mode = P->time_opts.log_mode | __IPF_MODE__;
					P->time_opts.log_ipf = atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg, "Unknow log file print criteria \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "EnFiles")) {
				if(!strcmp(xmlAttribute(s_elem, "type"), "No")) {
					P->time_opts.energy_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "FPS")) {
					P->time_opts.energy_mode = P->time_opts.energy_mode | __FPS_MODE__;
					P->time_opts.energy_fps = atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "IPF")) {
					P->time_opts.energy_mode = P->time_opts.energy_mode | __IPF_MODE__;
					P->time_opts.energy_ipf = atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg, "Unknow energy report file print criteria \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "BoundsFiles")) {
				if(!strcmp(xmlAttribute(s_elem, "type"), "No")) {
					P->time_opts.bounds_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "FPS")) {
					P->time_opts.bounds_mode = P->time_opts.bounds_mode | __FPS_MODE__;
					P->time_opts.bounds_fps = atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "IPF")) {
					P->time_opts.bounds_mode = P->time_opts.bounds_mode | __IPF_MODE__;
					P->time_opts.bounds_ipf = atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg, "Unknow fluid bounds file print criteria \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Output")) {
				if(!strcmp(xmlAttribute(s_elem, "type"), "No")) {
					P->time_opts.output_mode = __NO_OUTPUT_MODE__;
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "FPS")) {
					P->time_opts.output_mode = P->time_opts.output_mode | __FPS_MODE__;
					P->time_opts.output_fps = atof(xmlAttribute(s_elem, "value"));
				}
				else if(!strcmp(xmlAttribute(s_elem, "type"), "IPF")) {
					P->time_opts.output_mode = P->time_opts.output_mode | __IPF_MODE__;
					P->time_opts.output_ipf = atoi(xmlAttribute(s_elem, "value"));
				}
				else {
	                sprintf(msg, "Unknow output file print criteria \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid options are:\n");
	                S->addMessage(0, "\t\tNo\n");
	                S->addMessage(0, "\t\tFPS\n");
	                S->addMessage(0, "\t\tIPF\n");
	                return true;
				}
				if(!strcmp(xmlAttribute(s_elem, "format"), "H5Part")) {
					#ifdef HAVE_H5PART
	                    P->time_opts.output_format = P->time_opts.output_format | __H5Part__;
					#else
	                    sprintf(msg, "The H5Part format is not supported (due to the compilation).\n");
	                    S->addMessageF(3, msg);
	                    return true;
					#endif
				}
				else if(!strcmp(xmlAttribute(s_elem, "format"), "VTK")) {
					#ifdef HAVE_VTK
						P->time_opts.output_format = P->time_opts.output_format | __VTK__;
					#else
	                    sprintf(msg, "The VTK format is not supported (due to the compilation).\n");
	                    S->addMessageF(3, msg);
	                    return true;
					#endif
				}
				else if(!strcmp(xmlAttribute(s_elem, "format"), "Tecplot")) {
					#ifdef HAVE_TECPLOT
						P->time_opts.output_format = P->time_opts.output_format | __TECPLOT__;
	                    sprintf(msg, "The Tecplot format is outdated.\n");
	                    S->addMessageF(2, msg);
					#else
	                    sprintf(msg, "The Tecplot format is not supported (due to the compilation).\n");
	                    S->addMessageF(3, msg);
	                    return true;
					#endif
				}
				else {
	                sprintf(msg, "Unknow output file format \"%s\"\n", xmlAttribute(s_elem, "type"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tThe valid formats are:\n");
	                S->addMessage(0, "\t\tH5Part\n");
	                S->addMessage(0, "\t\tVTK\n");
	                S->addMessage(0, "\t\tTecplot\n");
	                return true;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "ExcludeField")) {
	            excludeField(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "TimeStep")) {
				if(   !strcmp(xmlAttribute(s_elem, "value"), "Courant")
	               || !strcmp(xmlAttribute(s_elem, "value"), "Variable") ) {
					P->time_opts.dt_mode = __DT_VARIABLE__;
				}
				else if(   !strcmp(xmlAttribute(s_elem, "value"), "Fix")
	                    || !strcmp(xmlAttribute(s_elem, "value"), "Fixed") ) {
					P->time_opts.dt_mode = __DT_FIXCALCULATED__;
				}
				else {
					P->time_opts.dt_mode = __DT_FIX__;
					P->time_opts.dt = atof(xmlAttribute(s_elem, "value"));
					if(P->time_opts.dt <= 0.f){
	                    sprintf(msg, "%g s is not a valid time step\n", P->time_opts.dt);
	                    S->addMessageF(2, msg);
	                    S->addMessage(0, "\tAutomatic calculated time step will used\n");
	                    P->time_opts.dt_mode = __DT_FIXCALCULATED__;
					}
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "MinTimeStep")) {
	            P->time_opts.dt_min = atof(xmlAttribute(s_elem, "value"));
	            if(P->time_opts.dt_min < 0.f){
	                sprintf(msg, "The minimum time step is lower than 0 s, so it will not take effect\n");
	                S->addMessageF(2, msg);
	            }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "ClampVel")) {
	            if(!strcmp(xmlAttribute(s_elem, "value"), "True") ||
	               !strcmp(xmlAttribute(s_elem, "value"), "true"))
	            {
	                P->time_opts.velocity_clamp = true;
	                if(P->time_opts.dt_min <= 0.f){
	                    sprintf(msg, "Velocity clamping has been activated, but minimum time step is %g s\n", P->time_opts.dt_min);
	                    S->addMessageF(2, msg);
	                    S->addMessage(0, "\tVelocity clamping will not take effect until minimum time step is not greather than 0 s\n");
	                }
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "False") ||
	                    !strcmp(xmlAttribute(s_elem, "value"), "false"))
	            {
	                P->time_opts.velocity_clamp = false;
	            }
	            else{
	                sprintf(msg, "%s is not a valid velocity clamping option\n", xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\ttrue\n");
	                S->addMessage(0, "\t\tfalse\n");
	                return true;
	            }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Stabilization")) {
	            P->time_opts.stabilization_time = atof(xmlAttribute(s_elem, "value"));
	            if(P->time_opts.stabilization_time < 0.f){
	                sprintf(msg, "%g s is not a valid stabilization time\n", P->time_opts.stabilization_time);
	                S->addMessageF(2, msg);
	                S->addMessage(0, "\tNo stabilization time will set\n");
	                P->time_opts.stabilization_time = 0.f;
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

bool FileManager::parseSPH(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("SPH"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Get options
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Option"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
			if(!strcmp(xmlAttribute(s_elem, "name"), "gamma")) {
				P->SPH_opts.gamma = atof(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "g")) {
				P->SPH_opts.g.x = atof(xmlAttribute(s_elem, "x"));
				P->SPH_opts.g.y = atof(xmlAttribute(s_elem, "y"));
	            #ifdef HAVE_3D
	                P->SPH_opts.g.z = atof(xmlAttribute(s_elem, "z"));
	            #endif
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "hfac")) {
				P->SPH_opts.hfac = atof(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "deltar")) {
				P->SPH_opts.deltar.x = atof(xmlAttribute(s_elem, "x"));
	            P->SPH_opts.deltar.y = atof(xmlAttribute(s_elem, "y"));
	            #ifdef HAVE_3D
	                P->SPH_opts.deltar.z = atof(xmlAttribute(s_elem, "z"));
	            #endif
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "cs")) {
				P->SPH_opts.cs = atof(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "DivDt")) {
				P->SPH_opts.dt_divisor = atof(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "LLSteps")) {
				P->SPH_opts.link_list_steps = atoi(xmlAttribute(s_elem, "value"));
				if(P->SPH_opts.link_list_steps < 1){
	                S->addMessageF(2, "LLSteps minor than 1 (value will be ignored).\n");
	                P->SPH_opts.link_list_steps = 1;
				}
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "DensSteps")) {
				P->SPH_opts.dens_int_steps = atoi(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "DensBounds")) {
			    if(xmlHasAttribute(s_elem,"min"))
                    P->SPH_opts.rho_min = atof(xmlAttribute(s_elem, "min"));
			    if(xmlHasAttribute(s_elem,"max"))
                    P->SPH_opts.rho_max = atof(xmlAttribute(s_elem, "max"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Boundary")) {
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
	                sprintf(msg, "Unknow boundary condition \"%s\"\n", xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tElasticBounce\n");
	                S->addMessage(0, "\t\tFixedParticles\n");
	                S->addMessage(0, "\t\tDeLeffe\n");
	                return true;
	            }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "SlipCondition")) {
	            if(!strcmp(xmlAttribute(s_elem, "value"), "NoSlip")){
	                P->SPH_opts.slip_condition = 0;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "value"), "FreeSlip")){
	                P->SPH_opts.slip_condition = 1;
	            }
	            else{
	                sprintf(msg, "Unknow slip condition \"%s\".\n", xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tNoSlip\n");
	                S->addMessage(0, "\t\tFreeSlip\n");
	                return true;
	            }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "BoundElasticFactor")) {
			    P->SPH_opts.elastic_factor = atof(xmlAttribute(s_elem, "value"));
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "BoundDist")) {
			    P->SPH_opts.elastic_dist = fabs(atof(xmlAttribute(s_elem, "value")));
			    if(!strcmp(xmlAttribute(s_elem, "force"), "true") || !strcmp(xmlAttribute(s_elem, "force"), "True")){
                    P->SPH_opts.elastic_dist = - P->SPH_opts.elastic_dist;
			    }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Shepard")) {
			    if(!strcmp(xmlAttribute(s_elem, "value"), "None")) {
	                P->SPH_opts.has_shepard = 0;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "Force")) {
	                P->SPH_opts.has_shepard = 1;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "Dens")) {
	                P->SPH_opts.has_shepard = 2;
			    }
			    else if (!strcmp(xmlAttribute(s_elem, "value"), "ForceDens")) {
	                P->SPH_opts.has_shepard = 3;
			    }
			    else {
	                sprintf(msg, "Unknow shepard application environment \"%s\".\n", xmlAttribute(s_elem, "value"));
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\tValid options are:\n");
	                S->addMessage(0, "\t\tNone\n");
	                S->addMessage(0, "\t\tForce\n");
	                S->addMessage(0, "\t\tDens\n");
	                S->addMessage(0, "\t\tForceDens\n");
	                return true;
			    }
			}
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Domain")) {
			    P->SPH_opts.has_domain = true;
			    P->SPH_opts.domain_min.x = atof(xmlAttribute(s_elem, "x"));
			    P->SPH_opts.domain_min.y = atof(xmlAttribute(s_elem, "y"));
			    P->SPH_opts.domain_max.x = atof(xmlAttribute(s_elem, "l")) + P->SPH_opts.domain_min.x;
	            #ifdef HAVE_3D
	                P->SPH_opts.domain_min.z = atof(xmlAttribute(s_elem, "z"));
	                P->SPH_opts.domain_max.y = atof(xmlAttribute(s_elem, "d")) + P->SPH_opts.domain_min.y;
	                P->SPH_opts.domain_max.z = atof(xmlAttribute(s_elem, "h")) + P->SPH_opts.domain_min.z;
	            #else
	                P->SPH_opts.domain_max.y = atof(xmlAttribute(s_elem, "h")) + P->SPH_opts.domain_min.y;
	            #endif
	            if(!strcmp(xmlAttribute(s_elem, "move"), "true") || !strcmp(xmlAttribute(s_elem, "move"), "True")){
                    P->SPH_opts.domain_motion = true;
	            }
	            else if(!strcmp(xmlAttribute(s_elem, "move"), "false") || !strcmp(xmlAttribute(s_elem, "move"), "False")){
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
	            sprintf(msg, "Unknow SPH option\n\t\t%s\n", xmlAttribute(s_elem, "name"));
	            S->addMessageF(3, msg);
	            return true;
			}
	    }
	}
	return false;
}

bool FileManager::parseFluid(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Fluid"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Add a new fluid
	    P->addFluid();
	    // Reads the number of particles
	    P->fluids[P->n_fluids-1].n = atoi(xmlAttribute(elem, "n"));
	    // Get options
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Option"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
			if(!strcmp(xmlAttribute(s_elem, "name"), "gamma"))
				P->fluids[P->n_fluids-1].gamma   = atof(xmlAttribute(s_elem, "value"));
			else if(!strcmp(xmlAttribute(s_elem, "name"), "refd"))
				P->fluids[P->n_fluids-1].refd    = atof(xmlAttribute(s_elem, "value"));
			else if(!strcmp(xmlAttribute(s_elem, "name"), "Viscdyn"))
				P->fluids[P->n_fluids-1].visc_dyn = atof(xmlAttribute(s_elem, "value"));
			else if(!strcmp(xmlAttribute(s_elem, "name"), "alpha"))
				P->fluids[P->n_fluids-1].alpha   = atof(xmlAttribute(s_elem, "value"));
			else if(!strcmp(xmlAttribute(s_elem, "name"), "delta"))
				P->fluids[P->n_fluids-1].delta   = atof(xmlAttribute(s_elem, "value"));
			else{
	            sprintf(msg, "Unknow Fluid option \"%s\"\n", xmlAttribute(s_elem, "name"));
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
	    // Get load methods
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Script"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
			strcpy(P->fluids[P->n_fluids-1].Path, xmlAttribute(s_elem, "path"));
			strcpy(P->fluids[P->n_fluids-1].Script, xmlAttribute(s_elem, "script"));
	    }
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Load"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
			strcpy(P->fluids[P->n_fluids-1].path, xmlAttribute(s_elem, "file"));
	    }
	}
	return false;
}

bool FileManager::parseSensors(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Sensors"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Get output frecuency
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("FPS"));
	    bool haveFPS = false;
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveFPS = true;
	        P->SensorsParameters.fps = atof(xmlAttribute(s_elem, "value"));
	    }
		if(!haveFPS) {
	        sprintf(msg, "Sensors section without FPS output frecuency.\n");
	        S->addMessageF(2, msg);
		}
		// Look for script
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Script"));
	    bool haveScript = false;
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        haveScript = true;
	        strcpy(P->SensorsParameters.script, xmlAttribute(s_elem, "file"));
	    }
		if(!haveScript) {
	        sprintf(msg, "Sensors section without script.\n");
	        S->addMessageF(2, msg);
		}
	    // Look for sensors
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Sensor"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
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
	            sprintf(msg, "%s type of sensor unknow.\n", xmlAttribute(s_elem, "type"));
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

bool FileManager::parseMotions(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Movements"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Look for movements
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("Movement"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        // Build new movement
	        ProblemSetup::sphMoveParameters *Move = new ProblemSetup::sphMoveParameters();
	        // Get movement type
	        if(!strcmp(xmlAttribute(s_elem, "type"),"Quaternion")){
	            Move->MoveType = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"LIQuaternion")){
	            Move->MoveType = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"C1Quaternion")){
	            Move->MoveType = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "type"),"ScriptQuaternion")){
	            Move->MoveType = 3;
	        }
	        else{
	            sprintf(msg, "Unknow type of movement \"%s\"\n", xmlAttribute(s_elem, "type"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tLIQuaternion\n");
	            S->addMessage(0, "\t\tScriptQuaternion\n");
	            return true;
	        }
	        // Get definition file
	        strcpy(Move->defFile, xmlAttribute(s_elem, "file"));
	        // Add movement to list
	        P->MoveParameters.push_back(Move);
	    }
	}
	return false;
}

bool FileManager::parsePortals(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
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
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("Portal"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Look for in/out planes
	    DOMNodeList* ss_nodes = elem->getElementsByTagName(XMLString::transcode("Plane"));
	    if(!ss_nodes->getLength()){
	        continue;
	    }
	    if(ss_nodes->getLength() != 2){
	        sprintf(msg, "Invalid planes data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\t2 planes by portal must be stablished, but %u found.\n", ss_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    // -----------------------------------------------------
	    // In portal
	    // -----------------------------------------------------
	    hasNormal = false;
	    ss_node = ss_nodes->item(0);
	    ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Corner"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\t1 corner must be provided, but %u found.\n", s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    in.corner.x = atof(xmlAttribute(s_elem, "x"));
	    in.corner.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        in.corner.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Up"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\t1 up vector must be provided, but %u found.\n", s_nodes->getLength());
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
	        s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Side"));
	        if(s_nodes->getLength() != 1){
	            sprintf(msg, "Invalid first plane data.\n");
	            S->addMessageF(3, msg);
	            sprintf(msg, "\t1 side vector must be provided, but %u found.\n", s_nodes->getLength());
	            S->addMessage(0, msg);
	            continue;
	        }
	        s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	        in.side.x = atof(xmlAttribute(s_elem, "x"));
	        in.side.y = atof(xmlAttribute(s_elem, "y"));
	        in.side.z = atof(xmlAttribute(s_elem, "z"));
	    #else
	        in.side.x = 0.f;
	        in.side.y = 0.f;
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Normal"));
	    if(s_nodes->getLength() > 1) {
	        sprintf(msg, "Invalid first plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\tOnly 1 normal vector can be provided, but %u found.\n", s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    else if (s_nodes->getLength() == 1) {
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
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Corner"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid second plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\t1 corner must be provided, but %u found.\n", s_nodes->getLength());
	        S->addMessage(0, msg);
	        continue;
	    }
	    s_elem = dynamic_cast< xercesc::DOMElement* >( s_nodes->item(0) );
	    out.corner.x = atof(xmlAttribute(s_elem, "x"));
	    out.corner.y = atof(xmlAttribute(s_elem, "y"));
	    #ifdef HAVE_3D
	        out.corner.z = atof(xmlAttribute(s_elem, "z"));
	    #endif
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Up"));
	    if(s_nodes->getLength() != 1){
	        sprintf(msg, "Invalid second plane data.\n");
	        S->addMessageF(3, msg);
	        sprintf(msg, "\t1 up vector must be provided, but %u found.\n", s_nodes->getLength());
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
	        s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Side"));
	        if(s_nodes->getLength() != 1){
                sprintf(msg, "Invalid second plane data.\n");
                S->addMessageF(3, msg);
                sprintf(msg, "\t1 side vector must be provided, but %u found.\n", s_nodes->getLength());
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
	    s_nodes = ss_elem->getElementsByTagName(XMLString::transcode("Normal"));
	    if(s_nodes->getLength() > 1) {
            sprintf(msg, "Invalid second plane data.\n");
            S->addMessageF(3, msg);
            sprintf(msg, "\tOnly 1 normal vector can be provided, but %u found.\n", s_nodes->getLength());
            S->addMessage(0, msg);
	        continue;
	    }
	    else if (s_nodes->getLength() == 1) {
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
	    // -----------------------------------------------------
	    // Append new portal
	    // -----------------------------------------------------
	    ProblemSetup::sphPortal *Portal = new ProblemSetup::sphPortal();
	    Portal->in  = in;
	    Portal->out = out;
	    P->Portals.push_back(Portal);
	}
	return false;
}

bool FileManager::parseGhostParticles(DOMElement *root)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	char msg[1024]; strcpy(msg, "");
	DOMNodeList* nodes = root->getElementsByTagName(XMLString::transcode("GhostParticles"));
	for( XMLSize_t i=0; i<nodes->getLength();i++ ){
	    DOMNode* node = nodes->item(i);
	    if( node->getNodeType() != DOMNode::ELEMENT_NODE )
	        continue;
	    DOMElement* elem = dynamic_cast< xercesc::DOMElement* >( node );
	    // Get pressure model
	    DOMNodeList* s_nodes = elem->getElementsByTagName(XMLString::transcode("PressModel"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->GhostParticles.pressModel = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->GhostParticles.pressModel = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->GhostParticles.pressModel = 2;
	        }
	        else{
	            sprintf(msg, "Invalid %s pressure extension model.\n", xmlAttribute(s_elem, "value"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tASM\n");
	            S->addMessage(0, "\t\tSSM\n");
	            S->addMessage(0, "\t\tTakeda\n");
	            return true;
	        }
	    }
	    // Get normal velocity model
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("NormalUModel"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->GhostParticles.nVelModel = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->GhostParticles.nVelModel = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->GhostParticles.nVelModel = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "U0M")){
	            P->GhostParticles.nVelModel = 3;
	        }
	        else{
	            sprintf(msg, "Invalid %s normal velocity extension model.\n", xmlAttribute(s_elem, "value"));
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
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("TangentUModel"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        if(!strcmp(xmlAttribute(s_elem, "value"), "ASM")){
	            P->GhostParticles.tVelModel = 0;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "SSM")){
	            P->GhostParticles.tVelModel = 1;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "Takeda")){
	            P->GhostParticles.tVelModel = 2;
	        }
	        else if(!strcmp(xmlAttribute(s_elem, "value"), "U0M")){
	            P->GhostParticles.tVelModel = 3;
	        }
	        else{
	            sprintf(msg, "Invalid %s tangent velocity extension model.\n", xmlAttribute(s_elem, "value"));
	            S->addMessageF(3, msg);
	            S->addMessage(0, "\tValid options are:\n");
	            S->addMessage(0, "\t\tASM\n");
	            S->addMessage(0, "\t\tSSM\n");
	            S->addMessage(0, "\t\tTakeda\n");
	            S->addMessage(0, "\t\tU0M\n");
	            return true;
	        }
	    }
	    // Look for walls
	    s_nodes = elem->getElementsByTagName(XMLString::transcode("Wall"));
	    for( XMLSize_t j=0; j<s_nodes->getLength();j++ ){
	        DOMNode* s_node = s_nodes->item(j);
	        if( s_node->getNodeType() != DOMNode::ELEMENT_NODE )
	            continue;
	        DOMElement* s_elem = dynamic_cast< xercesc::DOMElement* >( s_node );
	        DOMNodeList* ss_nodes = s_elem->getElementsByTagName(XMLString::transcode("Vertex"));
	        #ifdef HAVE_3D
	            if( (ss_nodes->getLength() != 3) && (ss_nodes->getLength() != 4) ){
	                sprintf(msg, "Found wall with %u vertexes.\n", ss_nodes->getLength());
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\t3D simulations accepts walls of 3 or 4 vertexes\n");
	                return true;
	            }
	            if( ss_nodes->getLength() == 3 ){
	                vec p[3];
	                for( XMLSize_t k=0; k<ss_nodes->getLength();k++ ){
	                    DOMNode* ss_node = ss_nodes->item(k);
	                    if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                        continue;
	                    DOMElement* ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	                    p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                    p[k].y = atof(xmlAttribute(ss_elem, "y"));
	                    p[k].z = atof(xmlAttribute(ss_elem, "z"));
	                    p[k].w = 0.f;
	                }
	                P->GhostParticles.add(p[0], p[1], p[2]);
	            }
	            if( ss_nodes->getLength() == 4 ){
	                vec p[4];
	                for( XMLSize_t k=0; k<ss_nodes->getLength();k++ ){
	                    DOMNode* ss_node = ss_nodes->item(k);
	                    if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                        continue;
	                    DOMElement* ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	                    p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                    p[k].y = atof(xmlAttribute(ss_elem, "y"));
	                    p[k].z = atof(xmlAttribute(ss_elem, "z"));
	                    p[k].w = 0.f;
	                }
	                P->GhostParticles.add(p[0], p[1], p[2], p[3]);
	            }
	        #else
	            if(ss_nodes->getLength() != 2){
	                sprintf(msg, "Found wall with %u vertexes.\n", ss_nodes->getLength());
	                S->addMessageF(3, msg);
	                S->addMessage(0, "\t2D simulations accepts walls of 2 vertexes\n");
	                return true;
	            }
	            vec p[2];
	            for( XMLSize_t k=0; k<ss_nodes->getLength();k++ ){
	                DOMNode* ss_node = ss_nodes->item(k);
	                if( ss_node->getNodeType() != DOMNode::ELEMENT_NODE )
	                    continue;
	                DOMElement* ss_elem = dynamic_cast< xercesc::DOMElement* >( ss_node );
	                p[k].x = atof(xmlAttribute(ss_elem, "x"));
	                p[k].y = atof(xmlAttribute(ss_elem, "y"));
	            }
	            P->GhostParticles.add(p[0], p[1]);
	        #endif
	    }
	}
	return false;
}

}}  // namespace
