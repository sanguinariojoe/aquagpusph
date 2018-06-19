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
 * @brief Simulation configuration files manager.
 * (See Aqua::InputOutput::State for details)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fnmatch.h>
#include <map>
#include <limits>

#include <InputOutput/State.h>
#include <ScreenManager.h>
#include <TimeManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

#include <vector>
#include <deque>
#include <algorithm>
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
using namespace std;

namespace Aqua{ namespace InputOutput{

State::State()
    : _output_file(NULL)
{
    unsigned int i;
    struct lconv *lc;
    char *s, msg[256];
    ScreenManager *S = ScreenManager::singleton();

    // Set the decimal-point character (which is depending on the locale)
    s = setlocale(LC_NUMERIC, NULL);
    if(strcmp(s, "C")){
        sprintf(msg, "\"%s\" numeric locale found\n", s);
        S->addMessageF(L_INFO, msg);
        S->addMessage(L_DEBUG, "\tIt is replaced by \"C\"\n");
        setlocale(LC_NUMERIC, "C");
    }
    lc = localeconv();
    s = lc->decimal_point;
    if(strcmp(s, ".")){
        sprintf(msg, "\"%s\" decimal point character found\n", s);
        S->addMessageF(L_WARNING, msg);
        S->addMessage(L_DEBUG, "\tIt is replaced by \".\"\n");
        lc->decimal_point = ".";
    }
    s = lc->thousands_sep;
    if(strcmp(s, "")){
        sprintf(msg, "\"%s\" thousands separator character found\n", s);
        S->addMessageF(L_WARNING, msg);
        S->addMessage(L_DEBUG, "\tIt is removed\n");
        lc->thousands_sep = "";
    }

    // Start the XML parser
    try {
        XMLPlatformUtils::Initialize();
    }
    catch( XMLException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(L_ERROR, "XML toolkit initialization error.\n");
        char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        exit(EXIT_FAILURE);
    }
    catch( ... ){
        S->addMessageF(L_ERROR, "XML toolkit initialization error.\n");
        S->addMessage(L_DEBUG, "\tUnhandled exception\n");
        xmlClear();
        exit(EXIT_FAILURE);
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
        char* message = xmlS( e.getMessage() );
        S->addMessageF(L_ERROR, "XML toolkit exit error.\n");
        char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        exit(EXIT_FAILURE);
    }
    catch( ... ){
        S->addMessageF(L_ERROR, "XML toolkit exit error.\n");
        S->addMessage(L_DEBUG, "\tUnhandled exception\n");
        xmlClear();
        exit(EXIT_FAILURE);
    }

    delete[] _output_file;
}

bool State::save(ProblemSetup sim_data, std::vector<Particles*> savers)
{
    return write(_output_file, sim_data, savers);
}

bool State::load(std::string input_file, ProblemSetup &sim_data)
{
    return parse(input_file.c_str(), sim_data);
}

bool State::parse(const char* filepath,
                  ProblemSetup &sim_data,
                  const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    DOMNodeList* nodes = NULL;
    char msg[1024];
    strcpy(msg, "");

    sprintf(msg,
            "Parsing the XML file \"%s\" with prefix \"%s\"\n",
            filepath,
            prefix);
    S->addMessageF(L_INFO, msg);
    // Try to open as ascii file, just to know if the file already exist
    FILE *dummy=0; dummy = fopen(filepath, "r");
    if(!dummy){
        S->addMessageF(L_ERROR, "File inaccessible!\n");
        return true;
    }
    fclose(dummy);
    XercesDOMParser *parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Never);
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(filepath);
     DOMDocument* doc = parser->getDocument();
    DOMElement* root = doc->getDocumentElement();
    if( !root ){
        S->addMessageF(L_ERROR, "Empty XML file\n");
        return true;
    }
    // Look for XML files included to parse them before
    nodes = root->getElementsByTagName(xmlS("Include"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        // By default, the include statements are parsed at the very beginning.
        if(xmlHasAttribute(elem, "when")){
            if(strcmp(xmlAttribute(elem, "when"), "begin"))
                continue;
        }
        const char* included_file = xmlAttribute(elem, "file");
        char* included_prefix = (char*)prefix;
        if(xmlHasAttribute(elem, "prefix")){
            included_prefix = xmlAttribute(elem, "prefix");
        }
        if(parse(included_file, sim_data, (const char*)included_prefix)){
            xmlClear();
            return true;
        }
    }

    if(parseSettings(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseVariables(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseDefinitions(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseTools(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseReports(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseTiming(root, sim_data, prefix)){
        xmlClear();
        return true;
    }
    if(parseSets(root, sim_data, prefix)){
        xmlClear();
        return true;
    }

    nodes = root->getElementsByTagName(xmlS("Include"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        if(!xmlHasAttribute(elem, "when"))
            continue;
        if(strcmp(xmlAttribute(elem, "when"), "end"))
            continue;
        const char* included_file = xmlAttribute(elem, "file");
        char* included_prefix = (char*)prefix;
        if(xmlHasAttribute(elem, "prefix")){
            included_prefix = xmlAttribute(elem, "prefix");
        }
        if(parse(included_file, sim_data, (const char*)included_prefix)){
            xmlClear();
            return true;
        }
    }
    
    xmlClear();
    delete parser;
    return false;
}

bool State::parseSettings(DOMElement *root,
                          ProblemSetup &sim_data,
                          const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[1024]; strcpy(msg, "");
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Settings"));
    for(XMLSize_t i=0; i<nodes->getLength();i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes;

        s_nodes = elem->getElementsByTagName(xmlS("Verbose"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            sim_data.settings.verbose_level = atoi(xmlAttribute(s_elem, "level"));
        }

        s_nodes = elem->getElementsByTagName(xmlS("SaveOnFail"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            if(!strcmp(xmlAttribute(s_elem, "value"), "true") ||
               !strcmp(xmlAttribute(s_elem, "value"), "True") ||
               !strcmp(xmlAttribute(s_elem, "value"), "TRUE")){
                sim_data.settings.save_on_fail = true;
            }
            else{
                sim_data.settings.save_on_fail = false;
            }
        }

        s_nodes = elem->getElementsByTagName(xmlS("Device"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            sim_data.settings.platform_id = atoi(xmlAttribute(s_elem, "platform"));
            sim_data.settings.device_id   = atoi(xmlAttribute(s_elem, "device"));
            if(!strcmp("ALL", xmlAttribute(s_elem, "type")))
                sim_data.settings.device_type = CL_DEVICE_TYPE_ALL;
            else if(!strcmp("CPU", xmlAttribute(s_elem, "type")))
                sim_data.settings.device_type = CL_DEVICE_TYPE_CPU;
            else if(!strcmp("GPU", xmlAttribute(s_elem, "type")))
                sim_data.settings.device_type = CL_DEVICE_TYPE_GPU;
            else if(!strcmp("ACCELERATOR", xmlAttribute(s_elem, "type")))
                sim_data.settings.device_type = CL_DEVICE_TYPE_ACCELERATOR;
            else if(!strcmp("DEFAULT", xmlAttribute(s_elem, "type")))
                sim_data.settings.device_type = CL_DEVICE_TYPE_DEFAULT;
            else{
                sprintf(msg,
                        "Unknow \"%s\" type of device\n",
                        xmlAttribute(s_elem, "type"));
                S->addMessageF(L_ERROR, msg);
                S->addMessage(L_DEBUG, "\tThe valid options are:\n");
                S->addMessage(L_DEBUG, "\t\tALL\n");
                S->addMessage(L_DEBUG, "\t\tCPU\n");
                S->addMessage(L_DEBUG, "\t\tGPU\n");
                S->addMessage(L_DEBUG, "\t\tACCELERATOR\n");
                S->addMessage(L_DEBUG, "\t\tDEFAULT\n");
                return true;
            }
        }
        s_nodes = elem->getElementsByTagName(xmlS("RootPath"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            sim_data.settings.base_path = xmlAttribute(s_elem, "path");
        }
    }
    return false;
}

bool State::parseVariables(DOMElement *root,
                           ProblemSetup &sim_data,
                           const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Variables"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Variable"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);

            if(!strstr(xmlAttribute(s_elem, "type"), "*")){
                sim_data.variables.registerVariable(xmlAttribute(s_elem, "name"),
                                              xmlAttribute(s_elem, "type"),
                                              "1",
                                              xmlAttribute(s_elem, "value"));
            }
            else{
                sim_data.variables.registerVariable(xmlAttribute(s_elem, "name"),
                                              xmlAttribute(s_elem, "type"),
                                              xmlAttribute(s_elem, "length"),
                                              "NULL");
            }
        }
    }
    return false;
}

bool State::parseDefinitions(DOMElement *root,
                             ProblemSetup &sim_data,
                             const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Definitions"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Define"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            if(!xmlHasAttribute(s_elem, "name")){
                S->addMessageF(L_ERROR, "Found a definition without name\n");
                return true;
            }
            if(!xmlHasAttribute(s_elem, "value")){
                sim_data.definitions.define(xmlAttribute(s_elem, "name"),
                                      "",
                                      false);
                continue;
            }

            bool evaluate = false;
            if(!strcmp(xmlAttribute(s_elem, "evaluate"), "true") ||
               !strcmp(xmlAttribute(s_elem, "evaluate"), "True") ||
               !strcmp(xmlAttribute(s_elem, "evaluate"), "TRUE")){
                evaluate = true;
            }
            if(!strcmp(xmlAttribute(s_elem, "evaluate"), "yes") ||
               !strcmp(xmlAttribute(s_elem, "evaluate"), "Yes") ||
               !strcmp(xmlAttribute(s_elem, "evaluate"), "YES")){
                evaluate = true;
            }

            sim_data.definitions.define(xmlAttribute(s_elem, "name"),
                                  xmlAttribute(s_elem, "value"),
                                  evaluate);
        }
    }
    return false;
}


/// Helper storage for the functions _toolsList() and _toolsName()
static std::deque<unsigned int> _tool_places;

/** @brief Helper function to get a list of tool placements from a list of names
 * @param list List of tools, separated by commas
 * @param prefix prefix to become inserted at the beggining of the name of each
 * tool of the list
 * @return The positions of the tools
 * @warning This methos is not thread safe
 */
static std::deque<unsigned int> _toolsList(const char* list,
                                           ProblemSetup &sim_data,
                                           const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    // Clear the list of tools (from previous executions)
    _tool_places.clear();
    // Loop along the list items
    char *token = strtok((char*)list, ",");
    while(token != NULL)
    {
        // Insert the prefix at the beggining of the tool name
        char *toolname = (char*)malloc(
            (strlen(prefix) + strlen(token) + 1) * sizeof(char));
        if(!toolname){
            S->addMessageF(L_ERROR, "Failure allocating memory.\n");
            return _tool_places;
        }
        strcpy(toolname, prefix);
        strcat(toolname, token);
        // Look for the tool in the already defined ones
        unsigned int place;
        for(place = 0; place < sim_data.tools.size(); place++){
            if(!strcmp(sim_data.tools.at(place)->get("name"),
                       toolname))
            {
                _tool_places.push_back(place);
            }
        }
        free(toolname);
        toolname = NULL;
        // Move to the next tool of the list
        token = strtok(NULL, ",");
    }

    return _tool_places;
}

/** @brief Helper function to get a list of tool placements from a wildcard
 * @param name Wildcard formatted tool name
 * @param prefix prefix to become inserted at the beggining of the name of each
 * tool of the list
 * @return The positions of the tools
 * @warning This methos is not thread safe
 */
static std::deque<unsigned int> _toolsName(const char* name,
                                           ProblemSetup &sim_data,
                                           const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    // Clear the list of tools (from previous executions)
    _tool_places.clear();
    // Insert the prefix at the beggining of the tool name
    char *toolname = (char*)malloc(
        (strlen(prefix) + strlen(name) + 1) * sizeof(char));
    if(!toolname){
        S->addMessageF(L_ERROR, "Failure allocating memory.\n");
        return _tool_places;
    }
    strcpy(toolname, prefix);
    strcat(toolname, name);
    // Look for the patterns in the already defined tool names
    unsigned int place;
    for(place = 0; place < sim_data.tools.size(); place++){
        if(!fnmatch(toolname,
                    sim_data.tools.at(place)->get("name"),
                    0))
        {
            _tool_places.push_back(place);
        }
    }
    free(toolname);
    toolname = NULL;
    return _tool_places;
}

bool State::parseTools(DOMElement *root,
                       ProblemSetup &sim_data,
                       const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[1024]; strcpy(msg, "");
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Tools"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Tool"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            if(!xmlHasAttribute(s_elem, "name")){
                S->addMessageF(L_ERROR, "Found a tool without name\n");
                return true;
            }
            if(!xmlHasAttribute(s_elem, "type")){
                S->addMessageF(L_ERROR, "Found a tool without type\n");
                return true;
            }

            // Create the tool
            ProblemSetup::sphTool *tool = new ProblemSetup::sphTool();

            // Set the name (with prefix)
            size_t name_len = strlen(prefix) +
                              strlen(xmlAttribute(s_elem, "name")) + 1;
            char *name = (char*)malloc(name_len * sizeof(char));
            if(!name){
                S->addMessageF(L_ERROR, "Failure allocating memory for the name\n");
                return true;
            }
            strcpy(name, prefix);
            strcat(name, xmlAttribute(s_elem, "name"));
            tool->set("name", name);
            free(name);
            name = NULL;

            tool->set("type", xmlAttribute(s_elem, "type"));

            // Check if the conditions to add the tool are fulfilled
            if(xmlHasAttribute(s_elem, "ifdef")){
                if(!sim_data.definitions.isDefined(xmlAttribute(s_elem, "ifdef"))){
                    sprintf(msg,
                            "Ignoring the tool \"%s\" because \"%s\" has not been defined.\n",
                            tool->get("name"),
                            xmlAttribute(s_elem, "ifdef"));
                    S->addMessageF(L_WARNING, msg);
                    delete tool;
                    continue;
                }
            }
            else if(xmlHasAttribute(s_elem, "ifndef")){
                if(sim_data.definitions.isDefined(xmlAttribute(s_elem, "ifndef"))){
                    sprintf(msg,
                            "Ignoring the tool \"%s\" because \"%s\" has been defined.\n",
                            tool->get("name"),
                            xmlAttribute(s_elem, "ifndef"));
                    S->addMessageF(L_WARNING, msg);
                    delete tool;
                    continue;
                }
            }

            // Place the tool
            if(!xmlHasAttribute(s_elem, "action") ||
               !strcmp(xmlAttribute(s_elem, "action"), "add")){
                sim_data.tools.push_back(tool);
            }
            else if(!strcmp(xmlAttribute(s_elem, "action"), "insert") ||
                    !strcmp(xmlAttribute(s_elem, "action"), "try_insert")){
                unsigned int place;
                std::deque<unsigned int> places;
                std::deque<unsigned int> all_places;
                std::deque<unsigned int>::iterator it;

                bool try_insert = !strcmp(xmlAttribute(s_elem, "action"),
                                          "try_insert");

                if(xmlHasAttribute(s_elem, "at")){
                    places.push_back(atoi(xmlAttribute(s_elem, "at")));
                }
                else if(xmlHasAttribute(s_elem, "before")){
                    char *att_str = xmlAttribute(s_elem, "before");
                    if(strchr((const char*)att_str, ',')){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList((const char*)att_str, sim_data, "");
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tools cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Get just the first one
                        it = std::min_element(all_places.begin(), all_places.end());
                        places.push_back(*it);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName((const char*)att_str, sim_data, "");
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tool cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Deep copy the places
                        for(it = all_places.begin(); it != all_places.end(); ++it)
                            places.push_back(*it);
                    }
                }
                else if(xmlHasAttribute(s_elem, "after")){
                    const char *att_str = xmlAttribute(s_elem, "after");
                    if(strchr((const char*)att_str, ',')){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList((const char*)att_str, sim_data, "");
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tools cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Get just the last one (and insert after that)
                        it = std::max_element(all_places.begin(), all_places.end());
                        places.push_back(*it + 1);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName((const char*)att_str, sim_data, "");
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tool cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Deep copy the places (adding 1 to insert after that)
                        for(it = all_places.begin(); it != all_places.end(); ++it)
                            places.push_back(*it + 1);
                    }
                }
                else if(xmlHasAttribute(s_elem, "before_prefix")){
                    const char *att_str = xmlAttribute(s_elem, "before_prefix");
                    if(strchr((const char*)att_str, ',')){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList((const char*)att_str, sim_data, prefix);
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tools cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Get just the first one
                        it = std::min_element(all_places.begin(), all_places.end());
                        places.push_back(*it);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName((const char*)att_str, sim_data, prefix);
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tool cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Deep copy the places
                        for(it = all_places.begin(); it != all_places.end(); ++it)
                            places.push_back(*it);
                    }
                }
                else if(xmlHasAttribute(s_elem, "after_prefix")){
                    const char *att_str = xmlAttribute(s_elem, "after_prefix");
                    if(strchr((const char*)att_str, ',')){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList((const char*)att_str, sim_data, prefix);
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tools cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Get just the last one (and insert after that)
                        it = std::max_element(all_places.begin(), all_places.end());
                        places.push_back(*it + 1);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName((const char*)att_str, sim_data, prefix);
                        if(!all_places.size()){
                            sprintf(msg,
                                    "The tool \"%s\" must be inserted before \"%s\", but such tool cannot be found.\n",
                                    tool->get("name"),
                                    att_str);
                            delete tool;
                            if(try_insert){
                                S->addMessageF(L_WARNING, msg);
                                continue;
                            }
                            else{
                                S->addMessageF(L_ERROR, msg);
                                return true;
                            }
                        }
                        // Deep copy the places (adding 1 to insert after that)
                        for(it = all_places.begin(); it != all_places.end(); ++it)
                            places.push_back(*it + 1);
                    }
                }
                else{
                    sprintf(msg,
                            "Missed the place where the tool \"%s\" should be inserted.\n",
                            tool->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    S->addMessage(L_DEBUG, "Please set one of the following attributes:\n");
                    S->addMessage(L_DEBUG, "\t\"in\"\n");
                    S->addMessage(L_DEBUG, "\t\"before\"\n");
                    S->addMessage(L_DEBUG, "\t\"after\"\n");
                    S->addMessage(L_DEBUG, "\t\"before_prefix\"\n");
                    S->addMessage(L_DEBUG, "\t\"after_prefix\"\n");
                    return true;
                }
                // We cannot directly insert the tools, because the places would
                // change meanwhile, so better backward adding them
                for(place = places.size(); place > 0; place--){
                    sim_data.tools.insert(sim_data.tools.begin() + places.at(place - 1),
                                    tool);
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "action"), "remove") ||
                    !strcmp(xmlAttribute(s_elem, "action"), "try_remove")){
                bool try_remove = !strcmp(xmlAttribute(s_elem, "action"),
                                          "try_remove");
                unsigned int place;
                std::deque<unsigned int> places;
                // Get the places of the tools selected
                places = _toolsName((const char*)(tool->get("name")), sim_data, prefix);
                if(!places.size()){
                    sprintf(msg,
                            "Failure removing the tool \"%s\" (tool not found).\n",
                            tool->get("name"));
                    delete tool;
                    if(try_remove){
                        S->addMessageF(L_WARNING, msg);
                        continue;
                    }
                    else{
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                }
                // Delete the new tool (which is useless)
                delete tool;
                // Delete the tools in backward order
                for(place = places.size(); place > 0; place--){
                    if(sim_data.toolInstances(sim_data.tools.at(places.at(place - 1))) == 1){
                        // This is the last instance
                        delete sim_data.tools.at(places.at(place - 1));
                    }
                    // Drop the tool from the list
                    sim_data.tools.erase(sim_data.tools.begin() + places.at(place - 1));
                }
                continue;
            }
            else if(!strcmp(xmlAttribute(s_elem, "action"), "replace") ||
                    !strcmp(xmlAttribute(s_elem, "action"), "try_replace")){
                bool try_replace = !strcmp(xmlAttribute(s_elem, "action"),
                                          "try_replace");
                unsigned int place;
                std::deque<unsigned int> places;
                // Get the places
                places = _toolsName((const char*)(tool->get("name")), sim_data, prefix);
                if(!places.size()){
                    sprintf(msg,
                            "Failure replacing the tool \"%s\" (tool not found).\n",
                            tool->get("name"));
                    delete tool;
                    if(try_replace){
                        S->addMessageF(L_WARNING, msg);
                        continue;
                    }
                    else{
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                }
                // Replace the tools
                for(place = 0; place < places.size(); place++){
                    if(sim_data.toolInstances(sim_data.tools.at(places.at(place))) == 1){
                        // This is the last instance
                        delete sim_data.tools.at(places.at(place));
                    }
                    // Set the new tool
                    sim_data.tools.at(places.at(place)) = tool;
                }
            }
            else{
                sprintf(msg,
                        "Unknown \"action\" for the tool \"%s\".\n",
                        tool->get("name"));
                S->addMessageF(L_ERROR, msg);
                S->addMessage(L_DEBUG, "\tThe valid actions are:\n");
                S->addMessage(L_DEBUG, "\t\tadd\n");
                S->addMessage(L_DEBUG, "\t\tinsert\n");
                S->addMessage(L_DEBUG, "\t\treplace\n");
                S->addMessage(L_DEBUG, "\t\tremove\n");
                return true;
            }

            // Configure the tool
            if(!strcmp(xmlAttribute(s_elem, "type"), "kernel")){
                if(!xmlHasAttribute(s_elem, "path")){
                    sprintf(msg,
                            "Tool \"%s\" is of type \"kernel\", but \"path\" is not defined.\n",
                            tool->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                tool->set("path", xmlAttribute(s_elem, "path"));
                if(!xmlHasAttribute(s_elem, "entry_point")){
                    tool->set("entry_point", "entry");
                }
                else{
                    tool->set("entry_point", xmlAttribute(s_elem, "entry_point"));
                }
                if(!xmlHasAttribute(s_elem, "n")){
                    tool->set("n", "N");
                }
                else{
                    tool->set("n", xmlAttribute(s_elem, "n"));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "copy")){
                const char *atts[2] = {"in", "out"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        sprintf(msg,
                                "Tool \"%s\" is of type \"set\", but \"%s\" is not defined.\n",
                                tool->get("name"),
                                atts[k]);
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "python")){
                if(!xmlHasAttribute(s_elem, "path")){
                    sprintf(msg,
                            "Tool \"%s\" is of type \"python\", but \"path\" is not defined.\n",
                            tool->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                tool->set("path", xmlAttribute(s_elem, "path"));
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "set")){
                const char *atts[2] = {"in", "value"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        sprintf(msg,
                                "Tool \"%s\" is of type \"set\", but \"%s\" is not defined.\n",
                                tool->get("name"),
                                atts[k]);
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "set_scalar")){
                const char *atts[2] = {"in", "value"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        sprintf(msg,
                                "Tool \"%s\" is of type \"set\", but \"%s\" is not defined.\n",
                                tool->get("name"),
                                atts[k]);
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "reduction")){
                const char *atts[3] = {"in", "out", "null"};
                for(unsigned int k = 0; k < 3; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        sprintf(msg,
                                "Tool \"%s\" is of type \"reduction\", but \"%s\" is not defined.\n",
                                tool->get("name"),
                                atts[k]);
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
                if(!strcmp(xmlS(s_elem->getTextContent()), "")){
                    sprintf(msg,
                            "No operation specified for the reduction \"%s\".\n",
                            tool->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                tool->set("operation", xmlS(s_elem->getTextContent()));
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "link-list")){
                if(!xmlHasAttribute(s_elem, "in")){
                    tool->set("in", "r");
                    continue;
                }
                tool->set("in", xmlAttribute(s_elem, "in"));
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "radix-sort")){
                const char *atts[3] = {"in", "perm", "inv_perm"};
                for(unsigned int k = 0; k < 3; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        sprintf(msg,
                                "Tool \"%s\" is of type \"set\", but \"%s\" is not defined.\n",
                                tool->get("name"),
                                atts[k]);
                        S->addMessageF(L_ERROR, msg);
                        return true;
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "assert")){
                if(!xmlHasAttribute(s_elem, "condition")){
                    sprintf(msg,
                            "Tool \"%s\" is of type \"assert\", but \"condition\" is not defined.\n",
                            tool->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                tool->set("condition", xmlAttribute(s_elem, "condition"));
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "dummy")){
            	// Without options
            }
            else{
                sprintf(msg,
                        "Unknown \"type\" for the tool \"%s\".\n",
                        tool->get("name"));
                S->addMessageF(L_ERROR, msg);
                S->addMessage(L_DEBUG, "\tThe valid types are:\n");
                S->addMessage(L_DEBUG, "\t\tkernel\n");
                S->addMessage(L_DEBUG, "\t\tcopy\n");
                S->addMessage(L_DEBUG, "\t\tpython\n");
                S->addMessage(L_DEBUG, "\t\tset\n");
                S->addMessage(L_DEBUG, "\t\tset_scalar\n");
                S->addMessage(L_DEBUG, "\t\treduction\n");
                S->addMessage(L_DEBUG, "\t\tlink-list\n");
                S->addMessage(L_DEBUG, "\t\tradix-sort\n");
                S->addMessage(L_DEBUG, "\t\tdummy\n");
                return true;
            }
        }
    }
    return false;
}

bool State::parseTiming(DOMElement *root,
                        ProblemSetup &sim_data,
                        const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
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
            if(!strcmp(name, "End") || !strcmp(name, "SimulationStop")){
                const char *type = xmlAttribute(s_elem, "type");
                if(!strcmp(type, "Time") || !strcmp(type, "T")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __TIME_MODE__;
                    sim_data.time_opts.sim_end_time =
                        atof(xmlAttribute(s_elem, "value"));
                }
                else if(!strcmp(type, "Steps") || !strcmp(type, "S")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __ITER_MODE__;
                    sim_data.time_opts.sim_end_step =
                        atoi(xmlAttribute(s_elem, "value"));
                }
                else if(!strcmp(type, "Frames") || !strcmp(type, "F")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __FRAME_MODE__;
                    sim_data.time_opts.sim_end_frame =
                        atoi(xmlAttribute(s_elem, "value"));
                }
                else {
                    sprintf(msg,
                            "Unknow simulation stop criteria \"%s\"\n",
                            type);
                    S->addMessageF(L_ERROR, msg);
                    S->addMessage(L_DEBUG, "\tThe valid options are:\n");
                    S->addMessage(L_DEBUG, "\t\tTime\n");
                    S->addMessage(L_DEBUG, "\t\tSteps\n");
                    S->addMessage(L_DEBUG, "\t\tFrames\n");
                    return true;
                }
            }
            else if(!strcmp(name, "LogFile")){
                const char *type = xmlAttribute(s_elem, "type");
                if(!strcmp(type, "No")){
                    sim_data.time_opts.log_mode = __NO_OUTPUT_MODE__;
                }
                else if(!strcmp(type, "FPS")){
                    sim_data.time_opts.log_mode =
                        sim_data.time_opts.log_mode | __FPS_MODE__;
                    sim_data.time_opts.log_fps =
                        atof(xmlAttribute(s_elem, "value"));
                }
                else if(!strcmp(type, "IPF")){
                    sim_data.time_opts.log_mode =
                        sim_data.time_opts.log_mode | __IPF_MODE__;
                    sim_data.time_opts.log_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
                }
                else {
                    sprintf(msg,
                            "Unknow log file print criteria \"%s\"\n",
                            xmlAttribute(s_elem, "type"));
                    S->addMessageF(L_ERROR, msg);
                    S->addMessage(L_DEBUG, "\tThe valid options are:\n");
                    S->addMessage(L_DEBUG, "\t\tNo\n");
                    S->addMessage(L_DEBUG, "\t\tFPS\n");
                    S->addMessage(L_DEBUG, "\t\tIPF\n");
                    return true;
                }
            }
            else if(!strcmp(name, "Output")){
                const char *type = xmlAttribute(s_elem, "type");
                if(!strcmp(type, "No")){
                    sim_data.time_opts.output_mode = __NO_OUTPUT_MODE__;
                }
                else if(!strcmp(type, "FPS")){
                    sim_data.time_opts.output_mode =
                        sim_data.time_opts.output_mode | __FPS_MODE__;
                    sim_data.time_opts.output_fps =
                        atof(xmlAttribute(s_elem, "value"));
                }
                else if(!strcmp(type, "IPF")){
                    sim_data.time_opts.output_mode =
                        sim_data.time_opts.output_mode | __IPF_MODE__;
                    sim_data.time_opts.output_ipf =
                        atoi(xmlAttribute(s_elem, "value"));
                }
                else {
                    sprintf(msg,
                            "Unknow output file print criteria \"%s\"\n",
                            type);
                    S->addMessageF(L_ERROR, msg);
                    S->addMessage(L_DEBUG, "\tThe valid options are:\n");
                    S->addMessage(L_DEBUG, "\t\tNo\n");
                    S->addMessage(L_DEBUG, "\t\tFPS\n");
                    S->addMessage(L_DEBUG, "\t\tIPF\n");
                    return true;
                }
            }

            else {
                sprintf(msg, "Unknow timing option \"%s\"\n", xmlAttribute(s_elem, "name"));
                S->addMessageF(L_ERROR, msg);
                return true;
            }
        }
    }
    return false;
}

bool State::parseSets(DOMElement *root,
                      ProblemSetup &sim_data,
                      const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[1024]; strcpy(msg, "");
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("ParticlesSet"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        if(!xmlHasAttribute(elem, "n")){
            S->addMessageF(L_ERROR, "Found a particles set without \"n\" attribute.\n");
            return true;
        }

        ProblemSetup::sphParticlesSet *Set = new ProblemSetup::sphParticlesSet();
        Set->n(atoi(xmlAttribute(elem, "n")));

        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Scalar"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);

            const char *name = xmlAttribute(s_elem, "name");
            const char *value = xmlAttribute(s_elem, "value");
            Set->addScalar(name, value);
        }

        s_nodes = elem->getElementsByTagName(xmlS("Load"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            const char *path = xmlAttribute(s_elem, "file");
            const char *format = xmlAttribute(s_elem, "format");
            const char *fields = xmlAttribute(s_elem, "fields");
            Set->input(path, format, fields);
        }

        s_nodes = elem->getElementsByTagName(xmlS("Save"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            const char *path = xmlAttribute(s_elem, "file");
            const char *format = xmlAttribute(s_elem, "format");
            const char *fields = xmlAttribute(s_elem, "fields");
            Set->output(path, format, fields);
        }
        sim_data.sets.push_back(Set);
    }
    return false;
}

bool State::parseReports(DOMElement *root,
                         ProblemSetup &sim_data,
                         const char* prefix)
{
    ScreenManager *S = ScreenManager::singleton();
    char msg[1024]; strcpy(msg, "");
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Reports"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Report"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            if(!xmlHasAttribute(s_elem, "name")){
                S->addMessageF(L_ERROR, "Found a report without name\n");
                return true;
            }
            if(!xmlHasAttribute(s_elem, "type")){
                S->addMessageF(L_ERROR, "Found a report without type\n");
                return true;
            }

            // Create the report
            ProblemSetup::sphTool *report = new ProblemSetup::sphTool();

            // Set the name (with prefix)
            size_t name_len = strlen(prefix) +
                              strlen(xmlAttribute(s_elem, "name")) + 1;
            char *name = (char*)malloc(name_len * sizeof(char));
            if(!name){
                S->addMessageF(L_ERROR, "Failure allocating memory for the name\n");
                return true;
            }
            strcpy(name, prefix);
            strcat(name, xmlAttribute(s_elem, "name"));
            report->set("name", name);
            free(name);
            name = NULL;

            report->set("type", xmlAttribute(s_elem, "type"));
            sim_data.reports.push_back(report);

            // Configure the report
            if(!strcmp(xmlAttribute(s_elem, "type"), "screen")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    S->addMessageF(L_ERROR, "Found a \"screen\" report without fields\n");
                    return true;
                }
                report->set("fields", xmlAttribute(s_elem, "fields"));
                if(xmlHasAttribute(s_elem, "bold")){
                    report->set("bold", xmlAttribute(s_elem, "bold"));
                }
                else{
                    report->set("bold", "false");
                }
                if(xmlHasAttribute(s_elem, "color")){
                    report->set("color", xmlAttribute(s_elem, "color"));
                }
                else{
                    report->set("color", "white");
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "file")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    S->addMessageF(L_ERROR, "Found a \"file\" report without fields\n");
                    return true;
                }
                report->set("fields", xmlAttribute(s_elem, "fields"));
                if(!xmlHasAttribute(s_elem, "path")){
                    sprintf(msg,
                            "Report \"%s\" is of type \"file\", but the output \"path\" is not defined.\n",
                            report->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                report->set("path", xmlAttribute(s_elem, "path"));
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "particles")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    S->addMessageF(L_ERROR, "Found a \"particles\" report without fields\n");
                    return true;
                }
                report->set("fields", xmlAttribute(s_elem, "fields"));
                if(!xmlHasAttribute(s_elem, "path")){
                    sprintf(msg,
                            "Report \"%s\" is of type \"particles\", but the output \"path\" is not defined.\n",
                            report->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                report->set("path", xmlAttribute(s_elem, "path"));

                if(!xmlHasAttribute(s_elem, "set")){
                    sprintf(msg,
                            "Report \"%s\" is of type \"particles\", but the \"set\" is not defined.\n",
                            report->get("name"));
                    S->addMessageF(L_ERROR, msg);
                    return true;
                }
                report->set("set", xmlAttribute(s_elem, "set"));

                if(!xmlHasAttribute(s_elem, "ipf")){
                    report->set("ipf", "1");
                }
                else{
                    report->set("ipf", xmlAttribute(s_elem, "ipf"));
                }
                if(!xmlHasAttribute(s_elem, "fps")){
                    report->set("fps", "0.0");
                }
                else{
                    report->set("fps", xmlAttribute(s_elem, "fps"));
                }
            }
            else if(!strcmp(xmlAttribute(s_elem, "type"), "performance")){
                if(xmlHasAttribute(s_elem, "bold")){
                    report->set("bold", xmlAttribute(s_elem, "bold"));
                }
                else{
                    report->set("bold", "false");
                }
                if(xmlHasAttribute(s_elem, "color")){
                    report->set("color", xmlAttribute(s_elem, "color"));
                }
                else{
                    report->set("color", "white");
                }
                if(xmlHasAttribute(s_elem, "path")){
                    report->set("path", xmlAttribute(s_elem, "path"));
                }
                else{
                    report->set("path", "");
                }
            }
            else{
                sprintf(msg,
                        "Unknown \"type\" for the report \"%s\".\n",
                        report->get("name"));
                S->addMessageF(L_ERROR, msg);
                S->addMessage(L_DEBUG, "\tThe valid types are:\n");
                S->addMessage(L_DEBUG, "\t\tscreen\n");
                S->addMessage(L_DEBUG, "\t\tfile\n");
                S->addMessage(L_DEBUG, "\t\tparticles\n");
                S->addMessage(L_DEBUG, "\t\tperformance\n");
                return true;
            }
        }
    }
    return false;
}

bool State::write(const char* filepath,
                  ProblemSetup sim_data,
                  std::vector<Particles*> savers)
{
    DOMImplementation* impl;
    char msg[64 + strlen(filepath)];
    ScreenManager *S = ScreenManager::singleton();
    sprintf(msg, "Writing \"%s\" SPH state file...\n", filepath);
    S->addMessageF(L_INFO, msg);

    impl = DOMImplementationRegistry::getDOMImplementation(xmlS("Range"));
    DOMDocument* doc = impl->createDocument(
        NULL,
        xmlS("sphInput"),
        NULL);
    DOMElement* root = doc->getDocumentElement();

    if(writeSettings(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeVariables(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeDefinitions(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeTools(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeReports(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeTiming(doc, root, sim_data)){
        xmlClear();
        return true;
    }
    if(writeSets(doc, root, sim_data, savers)){
        xmlClear();
        return true;
    }

    // Save the XML document to a file
    impl = DOMImplementationRegistry::getDOMImplementation(xmlS("LS"));
    DOMLSSerializer *saver = ((DOMImplementationLS*)impl)->createLSSerializer();

    if(saver->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
        saver->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
    saver->setNewLine(xmlS("\r\n"));

    XMLFormatTarget *target = new LocalFileFormatTarget(filepath);
    // XMLFormatTarget *target = new StdOutFormatTarget();
    DOMLSOutput *output = ((DOMImplementationLS*)impl)->createLSOutput();
    output->setByteStream(target);
    output->setEncoding(xmlS("UTF-8"));

    try {
        saver->write(doc, output);
    }
    catch( XMLException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(L_ERROR, "XML toolkit writing error.\n");
        char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        exit(EXIT_FAILURE);
    }
    catch( DOMException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(L_ERROR, "XML DOM writing error.\n");
        char msg[strlen(message) + 3];
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        exit(EXIT_FAILURE);
    }
    catch( ... ){
        S->addMessageF(L_ERROR, "Writing error.\n");
        S->addMessage(L_DEBUG, "\tUnhandled exception\n");
        exit(EXIT_FAILURE);
    }

    target->flush();

    delete target;
    saver->release();
    output->release();
    doc->release();
    xmlClear();

    return false;
}

bool State::writeSettings(xercesc::DOMDocument* doc,
                          xercesc::DOMElement *root,
                          ProblemSetup sim_data)
{
    DOMElement *elem, *s_elem;
    char att[1024];

    elem = doc->createElement(xmlS("Settings"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("Verbose"));
    sprintf(att, "%d", sim_data.settings.verbose_level);
    s_elem->setAttribute(xmlS("level"), xmlS(att));
    elem->appendChild(s_elem);
    s_elem = doc->createElement(xmlS("Device"));
    sprintf(att, "%u", sim_data.settings.platform_id);
    s_elem->setAttribute(xmlS("platform"), xmlS(att));
    sprintf(att, "%u", sim_data.settings.device_id);
    s_elem->setAttribute(xmlS("device"), xmlS(att));
    switch(sim_data.settings.device_type){
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
    s_elem->setAttribute(xmlS("type"), xmlS(att));
    elem->appendChild(s_elem);
    return false;
}

bool State::writeVariables(xercesc::DOMDocument* doc,
                           xercesc::DOMElement *root,
                           ProblemSetup sim_data)
{
    unsigned int i;
    DOMElement *elem, *s_elem;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    elem = doc->createElement(xmlS("Variables"));
    root->appendChild(elem);

    std::deque<Variable*> vars = C->variables()->getAll();

    for(i = 0; i < vars.size(); i++){
        s_elem = doc->createElement(xmlS("Variable"));
        elem->appendChild(s_elem);
        Variable* var = vars.at(i);

        s_elem->setAttribute(xmlS("name"), xmlS(var->name()));
        const char* type = var->type();
        s_elem->setAttribute(xmlS("type"), xmlS(type));

        // Array variable
        if(strstr(type, "*")){
            size_t length =
                ((ArrayVariable*)var)->size() / Variables::typeToBytes(var->type());
            char length_txt[16];
            sprintf(length_txt, "%lu", length);
            s_elem->setAttribute(xmlS("length"), xmlS(length_txt));
            continue;
        }
        // Scalar variable
        char* value_txt = (char*)var->asString();
        if(value_txt[0] == '('){
            value_txt[0] = ' ';
        }
        if(value_txt[strlen(value_txt) - 1] == ')'){
            value_txt[strlen(value_txt) - 1] = ' ';
        }
        s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
    }

    return false;
}

bool State::writeDefinitions(xercesc::DOMDocument* doc,
                             xercesc::DOMElement *root,
                             ProblemSetup sim_data)
{
    unsigned int i;
    DOMElement *elem, *s_elem;

    elem = doc->createElement(xmlS("Definitions"));
    root->appendChild(elem);

    ProblemSetup::sphDefinitions defs = sim_data.definitions;

    for(i = 0; i < defs.names.size(); i++){
        s_elem = doc->createElement(xmlS("Define"));
        elem->appendChild(s_elem);

        s_elem->setAttribute(xmlS("name"),
                           xmlS(defs.names.at(i)));
        s_elem->setAttribute(xmlS("value"),
                           xmlS(defs.values.at(i)));
        if(defs.evaluations.at(i)){
            s_elem->setAttribute(xmlS("evaluate"),
                               xmlS("true"));
        }
        else{
            s_elem->setAttribute(xmlS("evaluate"),
                               xmlS("false"));
        }
    }

    return false;
}

bool State::writeTools(xercesc::DOMDocument* doc,
                       xercesc::DOMElement *root,
                       ProblemSetup sim_data)
{
    unsigned int i, j;
    DOMElement *elem, *s_elem;

    elem = doc->createElement(xmlS("Tools"));
    root->appendChild(elem);

    std::deque<ProblemSetup::sphTool*> tools = sim_data.tools;

    for(i = 0; i < tools.size(); i++){
        s_elem = doc->createElement(xmlS("Tool"));
        elem->appendChild(s_elem);

        ProblemSetup::sphTool* tool = tools.at(i);
        for(j = 0; j < tool->n(); j++){
            const char* name = tool->getName(j);
            const char* value = tool->get(j);
            if(!strcmp(name, "operation")){
                // The reduction operation is not an attribute, but a text
                s_elem->setTextContent(xmlS(value));
                continue;
            }
            s_elem->setAttribute(xmlS(name), xmlS(value));
        }
    }

    return false;
}

bool State::writeReports(xercesc::DOMDocument* doc,
                         xercesc::DOMElement *root,
                         ProblemSetup sim_data)
{
    unsigned int i, j;
    DOMElement *elem, *s_elem;

    elem = doc->createElement(xmlS("Reports"));
    root->appendChild(elem);

    std::deque<ProblemSetup::sphTool*> reports = sim_data.reports;

    for(i = 0; i < reports.size(); i++){
        s_elem = doc->createElement(xmlS("Report"));
        elem->appendChild(s_elem);

        ProblemSetup::sphTool* report = reports.at(i);
        for(j = 0; j < report->n(); j++){
            const char* name = report->getName(j);
            const char* value = report->get(j);
            s_elem->setAttribute(xmlS(name), xmlS(value));
        }
    }

    return false;
}

bool State::writeTiming(xercesc::DOMDocument* doc,
                        xercesc::DOMElement *root,
                        ProblemSetup sim_data)
{
    DOMElement *elem, *s_elem;
    char att[1024];
    TimeManager *T = TimeManager::singleton();

    elem = doc->createElement(xmlS("Timing"));
    root->appendChild(elem);

    if(sim_data.time_opts.sim_end_mode & __TIME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Time"));
        sprintf(att, "%g", sim_data.time_opts.sim_end_time);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.sim_end_mode & __ITER_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Steps"));
        sprintf(att, "%d", sim_data.time_opts.sim_end_step);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.sim_end_mode & __FRAME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Frames"));
        sprintf(att, "%d", sim_data.time_opts.sim_end_frame);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    if(sim_data.time_opts.log_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("LogFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", sim_data.time_opts.log_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.log_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("LogFile"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", sim_data.time_opts.log_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    if(sim_data.time_opts.output_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        sprintf(att, "%g", sim_data.time_opts.output_fps);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.output_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        sprintf(att, "%d", sim_data.time_opts.output_ipf);
        s_elem->setAttribute(xmlS("value"), xmlS(att));
        elem->appendChild(s_elem);
    }

    return false;
}

bool State::writeSets(xercesc::DOMDocument* doc,
                      xercesc::DOMElement *root,
                      ProblemSetup sim_data,
                      std::vector<Particles*> savers)
{
    unsigned int i, j;
    char att[16];
    DOMElement *elem, *s_elem;

    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    Variables* vars = C->variables();
    std::deque<ProblemSetup::sphParticlesSet*> sets = sim_data.sets;

    for(i = 0; i < sets.size(); i++){
        elem = doc->createElement(xmlS("ParticlesSet"));
        sprintf(att, "%u", sim_data.sets.at(i)->n());
        elem->setAttribute(xmlS("n"), xmlS(att));
        root->appendChild(elem);

        for(j = 0; j < sets.at(i)->scalarNames().size(); j++){
            const char* name = sets.at(i)->scalarNames().at(j);
            s_elem = doc->createElement(xmlS("Scalar"));
            s_elem->setAttribute(xmlS("name"), xmlS(name));

            ArrayVariable* var = (ArrayVariable*)vars->get(name);
            char* value_txt = (char*)var->asString(i);
            if(!value_txt){
                return true;
            }
            if(value_txt[0] == '('){
                value_txt[0] = ' ';
            }
            if(value_txt[strlen(value_txt) - 1] == ')'){
                value_txt[strlen(value_txt) - 1] = ' ';
            }
            s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
            elem->appendChild(s_elem);
        }

        unsigned int n = 1;
        char *fields = new char[n];
        strcpy(fields, "");
        for(j = 0; j < sets.at(i)->outputFields().size(); j++){
            const char* field = sets.at(i)->outputFields().at(j);
            n += strlen(field) + 1;
            char *backup = fields;
            fields = new char[n];
            strcpy(fields, backup);
            delete[] backup; backup = NULL;
            strcat(fields, field);
            if(j < sets.at(i)->outputFields().size() - 1)
                strcat(fields, ",");
        }
        s_elem = doc->createElement(xmlS("Load"));
        s_elem->setAttribute(xmlS("file"), xmlS(savers.at(i)->file()));
        s_elem->setAttribute(xmlS("format"), xmlS(sets.at(i)->outputFormat()));
        s_elem->setAttribute(xmlS("fields"), xmlS(fields));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Save"));
        s_elem->setAttribute(xmlS("file"), xmlS(sets.at(i)->outputPath()));
        s_elem->setAttribute(xmlS("format"), xmlS(sets.at(i)->outputFormat()));
        s_elem->setAttribute(xmlS("fields"), xmlS(fields));
        elem->appendChild(s_elem);
        delete[] fields; fields=NULL;
    }

    return false;
}

}}  // namespace
