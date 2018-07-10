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

#include <fnmatch.h>
#include <map>
#include <limits>
#include <deque>
#include <algorithm>

#include <InputOutput/State.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>


static std::deque<std::string> cpp_str;
static std::deque<XMLCh*> xml_str;

static std::string xmlTranscode(const XMLCh *txt)
{
    std::string str = std::string(xercesc::XMLString::transcode(txt));
    cpp_str.push_back(str);
    return str;
}

static XMLCh *xmlTranscode(std::string txt)
{
    XMLCh *str = xercesc::XMLString::transcode(txt.c_str());
    xml_str.push_back(str);
    return str;
}

static void xmlClear()
{
    unsigned int i;
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
{
    unsigned int i;
    struct lconv *lc;
    char *s;

    // Set the decimal-point character (which is depending on the locale)
    s = setlocale(LC_NUMERIC, NULL);
    if(strcmp(s, "C")){
        std::ostringstream msg;
        msg << "\"" << s << "\" numeric locale found" << std::endl;
        LOG(L_INFO, msg.str());
        LOG0(L_DEBUG, "\tIt is replaced by \"C\"\n");
        setlocale(LC_NUMERIC, "C");
    }
    lc = localeconv();
    s = lc->decimal_point;
    if(strcmp(s, ".")){
        std::ostringstream msg;
        msg << "\"" << s << "\" decimal point character found" << std::endl;
        LOG(L_WARNING, msg.str());
        LOG0(L_DEBUG, "\tIt is replaced by \".\"\n");
        lc->decimal_point = ".";
    }
    s = lc->thousands_sep;
    if(strcmp(s, "")){
        std::ostringstream msg;
        msg << "\"" << s << "\" thousands separator character found" << std::endl;
        LOG(L_WARNING, msg.str());
        LOG0(L_DEBUG, "\tIt is removed\n");
        lc->thousands_sep = "";
    }

    // Start the XML parser
    try {
        XMLPlatformUtils::Initialize();
    }
    catch( XMLException& e ){
        std::ostringstream msg;
        LOG(L_ERROR, "XML toolkit initialization error.\n");
        msg << "\t" << xmlS(e.getMessage()) << std::endl;
        LOG0(L_DEBUG, msg.str());
        xmlClear();
        throw;
    }
    catch( ... ){
        LOG(L_ERROR, "XML toolkit initialization error.\n");
        LOG0(L_DEBUG, "\tUnhandled exception\n");
        xmlClear();
        throw;
    }

    // Look ofr the first available file place
    i = 0;
    std::ostringstream file_name;
    while(true){
        file_name.str("");
        file_name << "AQUAgpusph.save." << i << ".xml";
        std::ifstream f(file_name.str());
        if(f.is_open()){
            // The file already exist, look for another one
            i++;
            f.close();
            continue;
        }
        break;
    }
    _output_file = file_name.str();
}

State::~State()
{
    unsigned int i;
    // Terminate Xerces
    try{
        XMLPlatformUtils::Terminate();
    }
    catch( xercesc::XMLException& e ){
        std::ostringstream msg;
        LOG(L_ERROR, "XML toolkit exit error.\n");
        msg << "\t" << xmlS(e.getMessage()) << std::endl;
        LOG0(L_DEBUG, msg.str());
        xmlClear();
        throw;
    }
    catch( ... ){
        LOG(L_ERROR, "XML toolkit exit error.\n");
        LOG0(L_DEBUG, "\tUnhandled exception\n");
        xmlClear();
        throw;
    }
}

void State::save(ProblemSetup &sim_data, std::vector<Particles*> savers)
{
    return write(_output_file, sim_data, savers);
}

void State::load(std::string input_file, ProblemSetup &sim_data)
{
    return parse(input_file, sim_data);
}

void State::parse(std::string filepath,
                  ProblemSetup &sim_data,
                  std::string prefix)
{
    DOMNodeList* nodes = NULL;
    std::ostringstream msg;
    msg << "Parsing the XML file \"" << filepath
        << "\" with prefix \"" << prefix << "\"" << std::endl;
    LOG(L_INFO, msg.str());

    // Try to open as ascii file, just to know if the file already exist
    std::ifstream f(filepath);
    if(!f.is_open()){
        LOG(L_ERROR, "File inaccessible!\n");
        throw std::runtime_error("File inaccessible!");
    }
    f.close();

    // Now we can proceed to properly parse the XML file
    XercesDOMParser *parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Never);
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(filepath.c_str());
    DOMDocument* doc = parser->getDocument();
    DOMElement* root = doc->getDocumentElement();
    if( !root ){
        LOG(L_ERROR, "Empty XML file\n");
        throw std::runtime_error("Empty XML file");
    }

    // Parse <Include> tags to recursively load linked XML files
    nodes = root->getElementsByTagName(xmlS("Include"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        // By default, the include statements are parsed at the very beginning.
        if(xmlHasAttribute(elem, "when")){
            if(xmlAttribute(elem, "when").compare("begin"))
                continue;
        }
        std::string included_file = xmlAttribute(elem, "file");
        std::string included_prefix = prefix;
        if(xmlHasAttribute(elem, "prefix")){
            included_prefix = xmlAttribute(elem, "prefix");
        }
        try {
            parse(included_file, sim_data, included_prefix);
        } catch (...) {
            xmlClear();
            throw;
        }
    }

    // Parse the file itself
    try {
        parseSettings(root, sim_data, prefix);
        parseVariables(root, sim_data, prefix);
        parseDefinitions(root, sim_data, prefix);
        parseTools(root, sim_data, prefix);
        parseReports(root, sim_data, prefix);
        parseTiming(root, sim_data, prefix);
        parseSets(root, sim_data, prefix);
    } catch (...) {
        xmlClear();
        throw;
    }

    // Parse <Include> tags to recursively load linked XML files.
    nodes = root->getElementsByTagName(xmlS("Include"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        if(!xmlHasAttribute(elem, "when"))
            continue;
        if(xmlAttribute(elem, "when").compare("end"))
            continue;
        std::string included_file = xmlAttribute(elem, "file");
        std::string included_prefix = prefix;
        if(xmlHasAttribute(elem, "prefix")){
            included_prefix = xmlAttribute(elem, "prefix");
        }
        try {
            parse(included_file, sim_data, included_prefix);
        } catch (...) {
            xmlClear();
            throw;
        }
    }
    
    xmlClear();
    delete parser;
}

void State::parseSettings(DOMElement *root,
                          ProblemSetup &sim_data,
                          std::string prefix)
{
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Settings"));
    for(XMLSize_t i=0; i<nodes->getLength();i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        DOMNodeList* s_nodes;

        s_nodes = elem->getElementsByTagName(xmlS("SaveOnFail"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            if(!xmlAttribute(s_elem, "value").compare("true") ||
               !xmlAttribute(s_elem, "value").compare("True") ||
               !xmlAttribute(s_elem, "value").compare("TRUE")){
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
            sim_data.settings.platform_id = std::stoi(xmlAttribute(s_elem, "platform"));
            sim_data.settings.device_id   = std::stoi(xmlAttribute(s_elem, "device"));
            if(!xmlAttribute(s_elem, "type").compare("ALL"))
                sim_data.settings.device_type = CL_DEVICE_TYPE_ALL;
            else if(!xmlAttribute(s_elem, "type").compare("CPU"))
                sim_data.settings.device_type = CL_DEVICE_TYPE_CPU;
            else if(!xmlAttribute(s_elem, "type").compare("GPU"))
                sim_data.settings.device_type = CL_DEVICE_TYPE_GPU;
            else if(!xmlAttribute(s_elem, "type").compare("ACCELERATOR"))
                sim_data.settings.device_type = CL_DEVICE_TYPE_ACCELERATOR;
            else if(!xmlAttribute(s_elem, "type").compare("DEFAULT"))
                sim_data.settings.device_type = CL_DEVICE_TYPE_DEFAULT;
            else{
                std::ostringstream msg;
                msg << "Unknow \"" << xmlAttribute(s_elem, "type")
                    << "\" type of device" << std::endl;
                LOG(L_ERROR, msg.str());
                LOG0(L_DEBUG, "\tThe valid options are:\n");
                LOG0(L_DEBUG, "\t\tALL\n");
                LOG0(L_DEBUG, "\t\tCPU\n");
                LOG0(L_DEBUG, "\t\tGPU\n");
                LOG0(L_DEBUG, "\t\tACCELERATOR\n");
                LOG0(L_DEBUG, "\t\tDEFAULT\n");
                throw std::runtime_error("Invalid device type");
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
}

void State::parseVariables(DOMElement *root,
                           ProblemSetup &sim_data,
                           std::string prefix)
{
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

            if(xmlAttribute(s_elem, "type").find('*') == std::string::npos){
                sim_data.variables.registerVariable(xmlAttribute(s_elem, "name"),
                                                    xmlAttribute(s_elem, "type"),
                                                    "1",
                                                    xmlAttribute(s_elem, "value"));
            }
            else{
                sim_data.variables.registerVariable(xmlAttribute(s_elem, "name"),
                                                    xmlAttribute(s_elem, "type"),
                                                    xmlAttribute(s_elem, "length"),
                                                    "");
            }
        }
    }
    return;
}

void State::parseDefinitions(DOMElement *root,
                             ProblemSetup &sim_data,
                             std::string prefix)
{
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
                LOG(L_ERROR, "Found a definition without name\n");
                throw std::runtime_error("Name shall be specified for definitions");
            }
            if(!xmlHasAttribute(s_elem, "value")){
                sim_data.definitions.define(xmlAttribute(s_elem, "name"),
                                            "",
                                            false);
                continue;
            }

            bool evaluate = false;
            if(!xmlAttribute(s_elem, "evaluate").compare("true") ||
               !xmlAttribute(s_elem, "evaluate").compare("True") ||
               !xmlAttribute(s_elem, "evaluate").compare("TRUE")){
                evaluate = true;
            }
            if(!xmlAttribute(s_elem, "evaluate").compare("yes") ||
               !xmlAttribute(s_elem, "evaluate").compare("Yes") ||
               !xmlAttribute(s_elem, "evaluate").compare("YES")){
                evaluate = true;
            }

            sim_data.definitions.define(xmlAttribute(s_elem, "name"),
                                        xmlAttribute(s_elem, "value"),
                                        evaluate);
        }
    }
    return;
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
static std::deque<unsigned int> _toolsList(std::string list,
                                           ProblemSetup &sim_data,
                                           std::string prefix)
{
    _tool_places.clear();

    std::istringstream f(list);
    std::string s;
    while (getline(f, s, ',')) {
        std::ostringstream toolname;
        toolname << prefix << s;
        // Look for the tool in the already defined ones
        unsigned int place;
        for(place = 0; place < sim_data.tools.size(); place++){
            if(!toolname.str().compare(sim_data.tools.at(place)->get("name")))
            {
                _tool_places.push_back(place);
            }
        }        
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
static std::deque<unsigned int> _toolsName(std::string name,
                                           ProblemSetup &sim_data,
                                           std::string prefix)
{
    _tool_places.clear();
    std::ostringstream toolname;
    toolname << prefix << name;

    // Look for the patterns in the already defined tool names
    unsigned int place;
    for(place = 0; place < sim_data.tools.size(); place++){
        if(!fnmatch(toolname.str().c_str(),
                    sim_data.tools.at(place)->get("name").c_str(),
                    0))
        {
            _tool_places.push_back(place);
        }
    }

    return _tool_places;
}

void State::parseTools(DOMElement *root,
                       ProblemSetup &sim_data,
                       std::string prefix)
{
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
                LOG(L_ERROR, "Found a tool without name\n");
                throw std::runtime_error("Name shall be defined for tools");
            }
            if(!xmlHasAttribute(s_elem, "type")){
                LOG(L_ERROR, "Found a tool without type\n");
                throw std::runtime_error("Type shall be defined for tools");
            }

            // Create the tool
            ProblemSetup::sphTool *tool = new ProblemSetup::sphTool();

            // Set the name (with prefix)
            std::ostringstream name;
            name << prefix << xmlAttribute(s_elem, "name");
            tool->set("name", name.str());
            tool->set("type", xmlAttribute(s_elem, "type"));

            // Check if the conditions to add the tool are fulfilled
            if(xmlHasAttribute(s_elem, "ifdef")){
                if(!sim_data.definitions.isDefined(xmlAttribute(s_elem, "ifdef"))){
                    std::ostringstream msg;
                    msg << "Ignoring the tool \"" << tool->get("name")
                        << "\" because \"" << xmlAttribute(s_elem, "ifdef")
                        << "\" has not been defined." << std::endl;
                    LOG(L_WARNING, msg.str());
                    delete tool;
                    continue;
                }
            }
            else if(xmlHasAttribute(s_elem, "ifndef")){
                if(sim_data.definitions.isDefined(xmlAttribute(s_elem, "ifndef"))){
                    std::ostringstream msg;
                    msg << "Ignoring the tool \"" << tool->get("name")
                        << "\" because \"" << xmlAttribute(s_elem, "ifndef")
                        << "\" has been defined." << std::endl;
                    LOG(L_WARNING, msg.str());
                    delete tool;
                    continue;
                }
            }

            // Place the tool
            if(!xmlHasAttribute(s_elem, "action") ||
               !xmlAttribute(s_elem, "action").compare("add")){
                sim_data.tools.push_back(tool);
            }
            else if(!xmlAttribute(s_elem, "action").compare("insert") ||
                    !xmlAttribute(s_elem, "action").compare("try_insert")){
                unsigned int place;
                std::deque<unsigned int> places;
                std::deque<unsigned int> all_places;
                std::deque<unsigned int>::iterator it;

                bool try_insert = !xmlAttribute(s_elem, "action").compare(
                    "try_insert");

                if(xmlHasAttribute(s_elem, "at")){
                    places.push_back(std::stoi(xmlAttribute(s_elem, "at")));
                }
                else if(xmlHasAttribute(s_elem, "before") ||
                        xmlHasAttribute(s_elem, "before_prefix")){
                    std::string att_str, att_prefix;
                    if(xmlHasAttribute(s_elem, "before")) {
                        att_str = xmlAttribute(s_elem, "before");
                        att_prefix = "";
                    }
                    else {
                        att_str = xmlAttribute(s_elem, "before_prefix");
                        att_prefix = prefix;
                    }
                    if(att_str.find(',') != std::string::npos){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList(att_str, sim_data, att_prefix);
                        if(!all_places.size()){
                            std::ostringstream msg;
                            msg << "The tool \"" << tool->get("name")
                                << "\" must be inserted before \"" << att_str
                                << "\", but such tools cannot be found." << std::endl;
                            delete tool;
                            if(try_insert){
                                LOG(L_WARNING, msg.str());
                                continue;
                            }
                            else{
                                LOG(L_ERROR, msg.str());
                                throw std::runtime_error("Refered tool cannot be found");
                            }
                        }
                        // Get just the first one
                        it = std::min_element(all_places.begin(), all_places.end());
                        places.push_back(*it);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName(att_str, sim_data, att_prefix);
                        if(!all_places.size()){
                            std::ostringstream msg;
                            msg << "The tool \"" << tool->get("name")
                                << "\" must be inserted before \"" << att_str
                                << "\", but such tools cannot be found." << std::endl;
                            delete tool;
                            if(try_insert){
                                LOG(L_WARNING, msg.str());
                                continue;
                            }
                            else{
                                LOG(L_ERROR, msg.str());
                                throw std::runtime_error("Refered tool cannot be found");
                            }
                        }
                        // Deep copy the places
                        for(auto place : all_places)
                            places.push_back(place);
                    }
                }
                else if(xmlHasAttribute(s_elem, "after") ||
                        xmlHasAttribute(s_elem, "after_prefix")){
                    std::string att_str, att_prefix;
                    if(xmlHasAttribute(s_elem, "after")) {
                        att_str = xmlAttribute(s_elem, "after");
                        att_prefix = "";
                    }
                    else {
                        att_str = xmlAttribute(s_elem, "after_prefix");
                        att_prefix = prefix;
                    }
                    if(att_str.find(',') != std::string::npos){
                        // It is a list of names. We must get all the matching
                        // places and select the most convenient one
                        all_places = _toolsList(att_str, sim_data, att_prefix);
                        if(!all_places.size()){
                            std::ostringstream msg;
                            msg << "The tool \"" << tool->get("name")
                                << "\" must be inserted after \"" << att_str
                                << "\", but such tools cannot be found." << std::endl;
                            delete tool;
                            if(try_insert){
                                LOG(L_WARNING, msg.str());
                                continue;
                            }
                            else{
                                LOG(L_ERROR, msg.str());
                                throw std::runtime_error("Refered tool cannot be found");
                            }
                        }
                        // Get just the last one (and insert after that)
                        it = std::max_element(all_places.begin(), all_places.end());
                        places.push_back(*it + 1);
                    }
                    else{
                        // We can treat the string as a wildcard. Right now we
                        // wanna get all the places matching the pattern
                        all_places = _toolsName(att_str, sim_data, att_prefix);
                        if(!all_places.size()){
                            std::ostringstream msg;
                            msg << "The tool \"" << tool->get("name")
                                << "\" must be inserted after \"" << att_str
                                << "\", but such tools cannot be found." << std::endl;
                            delete tool;
                            if(try_insert){
                                LOG(L_WARNING, msg.str());
                                continue;
                            }
                            else{
                                LOG(L_ERROR, msg.str());
                                throw std::runtime_error("Refered tool cannot be found");
                            }
                        }
                        // Deep copy the places (adding 1 to insert after that)
                        for(auto place : all_places)
                            places.push_back(place + 1);
                    }
                }
                else{
                    std::ostringstream msg;
                    msg << "Missed the place where the tool \"" << tool->get("name")
                        << "\" should be inserted." << std::endl;
                    LOG(L_ERROR, msg.str());
                    LOG0(L_DEBUG, "Please set one of the following attributes:\n");
                    LOG0(L_DEBUG, "\t\"in\"\n");
                    LOG0(L_DEBUG, "\t\"before\"\n");
                    LOG0(L_DEBUG, "\t\"after\"\n");
                    LOG0(L_DEBUG, "\t\"before_prefix\"\n");
                    LOG0(L_DEBUG, "\t\"after_prefix\"\n");
                    throw std::runtime_error("A way to insert the tool shall be specified");
                }
                // We cannot directly insert the tools, because the places would
                // change meanwhile, so better backward adding them
                for(place = places.size(); place > 0; place--){
                    sim_data.tools.insert(sim_data.tools.begin() + places.at(place - 1),
                                    tool);
                }
            }
            else if(!xmlAttribute(s_elem, "action").compare("remove") ||
                    !xmlAttribute(s_elem, "action").compare("try_remove")){
                bool try_remove = !xmlAttribute(s_elem, "action").compare(
                    "try_remove");
                unsigned int place;
                std::deque<unsigned int> places;
                // Get the places of the tools selected
                places = _toolsName(tool->get("name"), sim_data, prefix);
                if(!places.size()){
                    std::ostringstream msg;
                    msg << "Failure removing the tool \"" << tool->get("name")
                        << "\". No such tool." << std::endl;
                    delete tool;
                    if(try_remove){
                        LOG(L_WARNING, msg.str());
                        continue;
                    }
                    else{
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("No such tool");
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
            else if(!xmlAttribute(s_elem, "action").compare("replace") ||
                    !xmlAttribute(s_elem, "action").compare("try_replace")){
                bool try_replace = !xmlAttribute(s_elem, "action").compare(
                    "try_replace");
                unsigned int place;
                std::deque<unsigned int> places;
                // Get the places
                places = _toolsName(tool->get("name"), sim_data, prefix);
                if(!places.size()){
                    std::ostringstream msg;
                    msg << "Failure replacing the tool \"" << tool->get("name")
                        << "\". No such tool." << std::endl;
                    delete tool;
                    if(try_replace){
                        LOG(L_WARNING, msg.str());
                        continue;
                    }
                    else{
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("No such tool");
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
                std::ostringstream msg;
                msg << "Unknown \"action\" for the tool \"" << tool->get("name")
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                LOG0(L_DEBUG, "\tThe valid actions are:\n");
                LOG0(L_DEBUG, "\t\tadd\n");
                LOG0(L_DEBUG, "\t\tinsert\n");
                LOG0(L_DEBUG, "\t\treplace\n");
                LOG0(L_DEBUG, "\t\tremove\n");
                throw std::runtime_error("Invalid action");
            }

            // Configure the tool
            if(!xmlAttribute(s_elem, "type").compare("kernel")){
                if(!xmlHasAttribute(s_elem, "path")){
                    std::ostringstream msg;
                    msg << "Tool \"" << tool->get("name")
                        << "\" is of type \"kernel\", but \"path\" is not defined."
                        << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Undefined OpenCL script path");
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
            else if(!xmlAttribute(s_elem, "type").compare("copy")){
                const char *atts[2] = {"in", "out"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        std::ostringstream msg;
                        msg << "Tool \"" << tool->get("name")
                            << "\" is of type \"copy\", but \"" << atts[k]
                            << "\" is not defined." << std::endl;
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("Missing attributes");
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!xmlAttribute(s_elem, "type").compare("python")){
                if(!xmlHasAttribute(s_elem, "path")){
                    std::ostringstream msg;
                    msg << "Tool \"" << tool->get("name")
                        << "\" is of type \"python\", but \"" << "path"
                        << "\" is not defined." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Undefined Python script path");
                }
                tool->set("path", xmlAttribute(s_elem, "path"));
            }
            else if(!xmlAttribute(s_elem, "type").compare("set")){
                const char *atts[2] = {"in", "value"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        std::ostringstream msg;
                        msg << "Tool \"" << tool->get("name")
                            << "\" is of type \"set\", but \"" << atts[k]
                            << "\" is not defined." << std::endl;
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("Missing attributes");
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!xmlAttribute(s_elem, "type").compare("set_scalar")){
                const char *atts[2] = {"in", "value"};
                for(unsigned int k = 0; k < 2; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        std::ostringstream msg;
                        msg << "Tool \"" << tool->get("name")
                            << "\" is of type \"set\", but \"" << atts[k]
                            << "\" is not defined." << std::endl;
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("Missing attributes");
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!xmlAttribute(s_elem, "type").compare("reduction")){
                const char *atts[3] = {"in", "out", "null"};
                for(unsigned int k = 0; k < 3; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        std::ostringstream msg;
                        msg << "Tool \"" << tool->get("name")
                            << "\" is of type \"reduction\", but \"" << atts[k]
                            << "\" is not defined." << std::endl;
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("Missing attributes");
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
                if(!xmlS(s_elem->getTextContent()).compare("")){
                    std::ostringstream msg;
                    msg << "No operation specified for the reduction \"" << tool->get("name")
                        << "\"." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Missing reduction operation");
                }
                tool->set("operation", xmlS(s_elem->getTextContent()));
            }
            else if(!xmlAttribute(s_elem, "type").compare("link-list")){
                if(!xmlHasAttribute(s_elem, "in")){
                    tool->set("in", "r");
                    continue;
                }
                tool->set("in", xmlAttribute(s_elem, "in"));
            }
            else if(!xmlAttribute(s_elem, "type").compare("radix-sort")){
                const char *atts[3] = {"in", "perm", "inv_perm"};
                for(unsigned int k = 0; k < 3; k++){
                    if(!xmlHasAttribute(s_elem, atts[k])){
                        std::ostringstream msg;
                        msg << "Tool \"" << tool->get("name")
                            << "\" is of type \"radix-sort\", but \"" << atts[k]
                            << "\" is not defined." << std::endl;
                        LOG(L_ERROR, msg.str());
                        throw std::runtime_error("Missing attribute");
                    }
                    tool->set(atts[k], xmlAttribute(s_elem, atts[k]));
                }
            }
            else if(!xmlAttribute(s_elem, "type").compare("assert")){
                if(!xmlHasAttribute(s_elem, "condition")){
                    std::ostringstream msg;
                    msg << "Tool \"" << tool->get("name")
                        << "\" is of type \"assert\", but \"" << "condition"
                        << "\" is not defined." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Missing attribute");
                }
                tool->set("condition", xmlAttribute(s_elem, "condition"));
            }
            else if(!xmlAttribute(s_elem, "type").compare("dummy")){
                // Without options
            }
            else{
                std::ostringstream msg;
                msg << "Unknown \"type\" for the tool \"" << tool->get("name")
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                LOG0(L_DEBUG, "\tThe valid types are:\n");
                LOG0(L_DEBUG, "\t\tkernel\n");
                LOG0(L_DEBUG, "\t\tcopy\n");
                LOG0(L_DEBUG, "\t\tpython\n");
                LOG0(L_DEBUG, "\t\tset\n");
                LOG0(L_DEBUG, "\t\tset_scalar\n");
                LOG0(L_DEBUG, "\t\treduction\n");
                LOG0(L_DEBUG, "\t\tlink-list\n");
                LOG0(L_DEBUG, "\t\tradix-sort\n");
                LOG0(L_DEBUG, "\t\tdummy\n");
                throw std::runtime_error("Unknown tool type");
            }
        }
    }
}

void State::parseTiming(DOMElement *root,
                        ProblemSetup &sim_data,
                        std::string prefix)
{
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
            std::string name = xmlAttribute(s_elem, "name");
            if(!name.compare("End") || !name.compare("SimulationStop")){
                std::string type = xmlAttribute(s_elem, "type");
                if(!type.compare("Time") || !type.compare("T")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __TIME_MODE__;
                    sim_data.time_opts.sim_end_time =
                        std::stof(xmlAttribute(s_elem, "value"));
                }
                else if(!type.compare("Steps") || !type.compare("S")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __ITER_MODE__;
                    sim_data.time_opts.sim_end_step =
                        std::stoi(xmlAttribute(s_elem, "value"));
                }
                else if(!type.compare("Frames") || !type.compare("F")){
                    sim_data.time_opts.sim_end_mode =
                        sim_data.time_opts.sim_end_mode | __FRAME_MODE__;
                    sim_data.time_opts.sim_end_frame =
                        std::stoi(xmlAttribute(s_elem, "value"));
                }
                else {
                    std::ostringstream msg;
                    msg << "Unknown simulation stop criteria \"" << type
                        << "\"." << std::endl;
                    LOG(L_ERROR, msg.str());
                    LOG0(L_DEBUG, "\tThe valid options are:\n");
                    LOG0(L_DEBUG, "\t\tTime\n");
                    LOG0(L_DEBUG, "\t\tSteps\n");
                    LOG0(L_DEBUG, "\t\tFrames\n");
                    throw std::runtime_error("Unknown simulation stop criteria");
                }
            }
            else if(!name.compare("Output")){
                std::string type = xmlAttribute(s_elem, "type");
                if(!type.compare("No")){
                    sim_data.time_opts.output_mode = __NO_OUTPUT_MODE__;
                }
                else if(!type.compare("FPS")){
                    sim_data.time_opts.output_mode =
                        sim_data.time_opts.output_mode | __FPS_MODE__;
                    sim_data.time_opts.output_fps =
                        std::stof(xmlAttribute(s_elem, "value"));
                }
                else if(!type.compare("IPF")){
                    sim_data.time_opts.output_mode =
                        sim_data.time_opts.output_mode | __IPF_MODE__;
                    sim_data.time_opts.output_ipf =
                        std::stoi(xmlAttribute(s_elem, "value"));
                }
                else {
                    std::ostringstream msg;
                    msg << "Unknow output file print criteria \"" << type
                        << "\"." << std::endl;
                    LOG(L_ERROR, msg.str());
                    LOG0(L_DEBUG, "\tThe valid options are:\n");
                    LOG0(L_DEBUG, "\t\tNo\n");
                    LOG0(L_DEBUG, "\t\tFPS\n");
                    LOG0(L_DEBUG, "\t\tIPF\n");
                    throw std::runtime_error("Unknow output file print criteria");
                }
            }

            else {
                std::ostringstream msg;
                msg << "Unknow timing option \"" << xmlAttribute(s_elem, "name")
                    << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                throw std::runtime_error("Unknow timing option");
            }
        }
    }
}

void State::parseSets(DOMElement *root,
                      ProblemSetup &sim_data,
                      std::string prefix)
{
    DOMNodeList* nodes = root->getElementsByTagName(xmlS("ParticlesSet"));
    for(XMLSize_t i=0; i<nodes->getLength(); i++){
        DOMNode* node = nodes->item(i);
        if(node->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
        DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
        if(!xmlHasAttribute(elem, "n")){
            LOG(L_ERROR, "Found a particles set without \"n\" attribute.\n");
            throw std::runtime_error("Missing number of particles in a set");
        }

        ProblemSetup::sphParticlesSet *set = new ProblemSetup::sphParticlesSet();
        set->n(std::stoi(xmlAttribute(elem, "n")));

        DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Scalar"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);

            std::string name = xmlAttribute(s_elem, "name");
            std::string value = xmlAttribute(s_elem, "value");
            set->addScalar(name, value);
        }

        s_nodes = elem->getElementsByTagName(xmlS("Load"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            std::string path = xmlAttribute(s_elem, "file");
            std::string format = xmlAttribute(s_elem, "format");
            std::string fields = xmlAttribute(s_elem, "fields");
            set->input(path, format, fields);
        }

        s_nodes = elem->getElementsByTagName(xmlS("Save"));
        for(XMLSize_t j=0; j<s_nodes->getLength(); j++){
            DOMNode* s_node = s_nodes->item(j);
            if(s_node->getNodeType() != DOMNode::ELEMENT_NODE)
                continue;
            DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
            std::string path = xmlAttribute(s_elem, "file");
            std::string format = xmlAttribute(s_elem, "format");
            std::string fields = xmlAttribute(s_elem, "fields");
            set->output(path, format, fields);
        }
        sim_data.sets.push_back(set);
    }
}

void State::parseReports(DOMElement *root,
                         ProblemSetup &sim_data,
                         std::string prefix)
{
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
                LOG(L_ERROR, "Found a report without name\n");
                throw std::runtime_error("Missing report name");
            }
            if(!xmlHasAttribute(s_elem, "type")){
                LOG(L_ERROR, "Found a report without type\n");
                throw std::runtime_error("Missing report type");
            }

            // Create the report
            ProblemSetup::sphTool *report = new ProblemSetup::sphTool();

            // Set the name (with prefix)
            std::ostringstream name;
            name << prefix << xmlAttribute(s_elem, "name");
            report->set("name", name.str());

            report->set("type", xmlAttribute(s_elem, "type"));
            sim_data.reports.push_back(report);

            // Configure the report
            if(!xmlAttribute(s_elem, "type").compare("screen")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    LOG(L_ERROR, "Found a \"screen\" report without fields\n");
                    throw std::runtime_error("Missing report fields");
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
            else if(!xmlAttribute(s_elem, "type").compare("file")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    LOG(L_ERROR, "Found a \"file\" report without fields\n");
                    throw std::runtime_error("Missing report fields");
                }
                report->set("fields", xmlAttribute(s_elem, "fields"));
                if(!xmlHasAttribute(s_elem, "path")){
                    std::ostringstream msg;
                    msg << "Report \"" << report->get("name")
                        << "\" is of type \"file\", but the output \"path\" is not defined." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Missing report file path");
                }
                report->set("path", xmlAttribute(s_elem, "path"));
            }
            else if(!xmlAttribute(s_elem, "type").compare("particles")){
                if(!xmlHasAttribute(s_elem, "fields")){
                    LOG(L_ERROR, "Found a \"particles\" report without fields\n");
                    throw std::runtime_error("Missing report fields");
                }
                report->set("fields", xmlAttribute(s_elem, "fields"));
                if(!xmlHasAttribute(s_elem, "path")){
                    std::ostringstream msg;
                    msg << "Report \"" << report->get("name")
                        << "\" is of type \"particles\", but the output \"path\" is not defined." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Missing report file path");
                }
                report->set("path", xmlAttribute(s_elem, "path"));

                if(!xmlHasAttribute(s_elem, "set")){
                    std::ostringstream msg;
                    msg << "Report \"" << report->get("name")
                        << "\" is of type \"particles\", but the output \"set\" is not defined." << std::endl;
                    LOG(L_ERROR, msg.str());
                    throw std::runtime_error("Missing report particles set");
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
            else if(!xmlAttribute(s_elem, "type").compare("performance")){
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
                std::ostringstream msg;
                msg << "Unknown \"type\" for the report \""
                    << report->get("name") << "\"." << std::endl;
                LOG(L_ERROR, msg.str());
                LOG0(L_DEBUG, "\tThe valid types are:\n");
                LOG0(L_DEBUG, "\t\tscreen\n");
                LOG0(L_DEBUG, "\t\tfile\n");
                LOG0(L_DEBUG, "\t\tparticles\n");
                LOG0(L_DEBUG, "\t\tperformance\n");
                throw std::runtime_error("Invalid report type");
            }
        }
    }
}

void State::write(std::string filepath,
                  ProblemSetup &sim_data,
                  std::vector<Particles*> savers)
{
    DOMImplementation* impl;
    std::ostringstream msg;
    msg << "Writing \"" << filepath
        << "\" SPH state file..." << std::endl;
    LOG(L_INFO, msg.str());

    impl = DOMImplementationRegistry::getDOMImplementation(xmlS("Range"));
    DOMDocument* doc = impl->createDocument(
        NULL,
        xmlS("sphInput"),
        NULL);
    DOMElement* root = doc->getDocumentElement();

    try {
        writeSettings(doc, root, sim_data);
        writeVariables(doc, root, sim_data);
        writeDefinitions(doc, root, sim_data);
        writeTools(doc, root, sim_data);
        writeReports(doc, root, sim_data);
        writeTiming(doc, root, sim_data);
        writeSets(doc, root, sim_data, savers);
    } catch (...) {
        xmlClear();
        throw;
    }

    // Save the XML document to a file
    impl = DOMImplementationRegistry::getDOMImplementation(xmlS("LS"));
    DOMLSSerializer *saver = ((DOMImplementationLS*)impl)->createLSSerializer();

    if(saver->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
        saver->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
    saver->setNewLine(xmlS("\r\n"));

    XMLFormatTarget *target = new LocalFileFormatTarget(filepath.c_str());
    // XMLFormatTarget *target = new StdOutFormatTarget();
    DOMLSOutput *output = ((DOMImplementationLS*)impl)->createLSOutput();
    output->setByteStream(target);
    output->setEncoding(xmlS("UTF-8"));

    try {
        saver->write(doc, output);
    }
    catch( XMLException& e ){
        LOG(L_ERROR, "XML toolkit writing error.\n");
        msg.str("");
        msg << "\t" << xmlS(e.getMessage()) << std::endl;
        LOG0(L_DEBUG, msg.str());
        xmlClear();
        throw std::runtime_error("Failure writing XML file");
    }
    catch( DOMException& e ){
        LOG(L_ERROR, "XML DOM writing error.\n");
        msg.str("");
        msg << "\t" << xmlS(e.getMessage()) << std::endl;
        LOG0(L_DEBUG, msg.str());
        xmlClear();
        throw std::runtime_error("XML DOM error while writing XML");
    }
    catch( ... ){
        LOG(L_ERROR, "Writing error.\n");
        LOG0(L_DEBUG, "\tUnhandled exception\n");
        throw std::runtime_error("Unhandled exception while writing XML");
    }

    target->flush();

    delete target;
    saver->release();
    output->release();
    doc->release();
    xmlClear();
}

void State::writeSettings(xercesc::DOMDocument* doc,
                          xercesc::DOMElement *root,
                          ProblemSetup &sim_data)
{
    DOMElement *elem, *s_elem;
    std::ostringstream att;

    elem = doc->createElement(xmlS("Settings"));
    root->appendChild(elem);

    s_elem = doc->createElement(xmlS("Device"));
    att << sim_data.settings.platform_id;
    s_elem->setAttribute(xmlS("platform"), xmlS(att.str()));
    att.str(""); att << sim_data.settings.device_id;
    s_elem->setAttribute(xmlS("device"), xmlS(att.str()));
    att.str("");
    switch(sim_data.settings.device_type){
        case CL_DEVICE_TYPE_ALL:
            att << "ALL";
            break;
        case CL_DEVICE_TYPE_CPU:
            att << "CPU";
            break;
        case CL_DEVICE_TYPE_GPU:
            att << "GPU";
            break;
        case CL_DEVICE_TYPE_ACCELERATOR:
            att << "ACCELERATOR";
            break;
        case CL_DEVICE_TYPE_DEFAULT:
            att << "DEFAULT";
            break;
    }
    s_elem->setAttribute(xmlS("type"), xmlS(att.str()));
    elem->appendChild(s_elem);
}

void State::writeVariables(xercesc::DOMDocument* doc,
                           xercesc::DOMElement *root,
                           ProblemSetup &sim_data)
{
    unsigned int i;
    DOMElement *elem, *s_elem;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    elem = doc->createElement(xmlS("Variables"));
    root->appendChild(elem);

    std::vector<Variable*> vars = C->variables()->getAll();
    for (auto var : vars) {
        s_elem = doc->createElement(xmlS("Variable"));
        elem->appendChild(s_elem);

        s_elem->setAttribute(xmlS("name"), xmlS(var->name()));
        std::string type = var->type();
        s_elem->setAttribute(xmlS("type"), xmlS(type));

        // Array variable
        if(type.find("*") != std::string::npos) {
            size_t length =
                ((ArrayVariable*)var)->size() / Variables::typeToBytes(var->type());
            std::ostringstream length_txt;
            length_txt << length;
            s_elem->setAttribute(xmlS("length"), xmlS(length_txt.str()));
            continue;
        }
        // Scalar variable
        std::string value_txt = var->asString();
        if(value_txt.at(0) == '('){
            value_txt.at(0) = ' ';
        }
        if(value_txt.back() == ')'){
            value_txt.back() = ' ';
        }
        s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
    }
}

void State::writeDefinitions(xercesc::DOMDocument* doc,
                             xercesc::DOMElement *root,
                             ProblemSetup &sim_data)
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
}

void State::writeTools(xercesc::DOMDocument* doc,
                       xercesc::DOMElement *root,
                       ProblemSetup &sim_data)
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
            std::string name = tool->getName(j);
            std::string value = tool->get(j);
            if(!name.compare("operation")){
                // The reduction operation is not an attribute, but a text
                s_elem->setTextContent(xmlS(value));
                continue;
            }
            s_elem->setAttribute(xmlS(name), xmlS(value));
        }
    }
}

void State::writeReports(xercesc::DOMDocument* doc,
                         xercesc::DOMElement *root,
                         ProblemSetup &sim_data)
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
            std::string name = report->getName(j);
            std::string value = report->get(j);
            s_elem->setAttribute(xmlS(name), xmlS(value));
        }
    }
}

void State::writeTiming(xercesc::DOMDocument* doc,
                        xercesc::DOMElement *root,
                        ProblemSetup &sim_data)
{
    DOMElement *elem, *s_elem;
    std::ostringstream att;

    elem = doc->createElement(xmlS("Timing"));
    root->appendChild(elem);

    if(sim_data.time_opts.sim_end_mode & __TIME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Time"));
        att.str(""); att << sim_data.time_opts.sim_end_time;
        s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.sim_end_mode & __ITER_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Steps"));
        att.str(""); att << sim_data.time_opts.sim_end_step;
        s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.sim_end_mode & __FRAME_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("End"));
        s_elem->setAttribute(xmlS("type"), xmlS("Frames"));
        att.str(""); att << sim_data.time_opts.sim_end_frame;
        s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
        elem->appendChild(s_elem);
    }

    if(sim_data.time_opts.output_mode & __FPS_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
        att.str(""); att << sim_data.time_opts.output_fps;
        s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
        elem->appendChild(s_elem);
    }
    if(sim_data.time_opts.output_mode & __IPF_MODE__){
        s_elem = doc->createElement(xmlS("Option"));
        s_elem->setAttribute(xmlS("name"), xmlS("Output"));
        s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
        att.str(""); att << sim_data.time_opts.output_ipf;
        s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
        elem->appendChild(s_elem);
    }
}

void State::writeSets(xercesc::DOMDocument* doc,
                      xercesc::DOMElement *root,
                      ProblemSetup &sim_data,
                      std::vector<Particles*> savers)
{
    unsigned int i, j;
    std::ostringstream att;
    DOMElement *elem, *s_elem;

    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    for(i = 0; i < sim_data.sets.size(); i++){
        elem = doc->createElement(xmlS("ParticlesSet"));
        att.str(""); att << sim_data.sets.at(i)->n();
        elem->setAttribute(xmlS("n"), xmlS(att.str()));
        root->appendChild(elem);

        for(j = 0; j < sim_data.sets.at(i)->scalarNames().size(); j++){
            std::string name = sim_data.sets.at(i)->scalarNames().at(j);
            s_elem = doc->createElement(xmlS("Scalar"));
            s_elem->setAttribute(xmlS("name"), xmlS(name));

            ArrayVariable* var = (ArrayVariable*)C->variables()->get(name);
            std::string value_txt = var->asString(i);
            if(value_txt.at(0) == '('){
                value_txt.at(0) = ' ';
            }
            if(value_txt.back() == ')'){
                value_txt.back() = ' ';
            }
            s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
            elem->appendChild(s_elem);
        }

        std::ostringstream fields;
        for(j = 0; j < sim_data.sets.at(i)->outputFields().size(); j++){
            std::string field = sim_data.sets.at(i)->outputFields().at(j);
            fields << field;
            if(j < sim_data.sets.at(i)->outputFields().size() - 1)
                fields << ",";            
        }
        s_elem = doc->createElement(xmlS("Load"));
        s_elem->setAttribute(xmlS("file"), xmlS(savers.at(i)->file()));
        s_elem->setAttribute(xmlS("format"), xmlS(sim_data.sets.at(i)->outputFormat()));
        s_elem->setAttribute(xmlS("fields"), xmlS(fields.str()));
        elem->appendChild(s_elem);

        s_elem = doc->createElement(xmlS("Save"));
        s_elem->setAttribute(xmlS("file"), xmlS(sim_data.sets.at(i)->outputPath()));
        s_elem->setAttribute(xmlS("format"), xmlS(sim_data.sets.at(i)->outputFormat()));
        s_elem->setAttribute(xmlS("fields"), xmlS(fields.str()));
        elem->appendChild(s_elem);
    }
}

}}  // namespace
