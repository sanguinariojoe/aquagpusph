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

#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <sstream>
#include <filesystem>
#include <system_error>

#include "aquagpusph/sphPrerequisites.hpp"
#include "State.hpp"
#include "Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

/** @brief Checks if two given strings match.
 * 
 * The first string may contain wildcard characters.
 * @param first First string, which may contain wildcards
 * @param second Second string
 * @return true if they match, false otherwise
 * @see https://www.geeksforgeeks.org/wildcard-character-matching/
 */
bool __match(const char* first, const char* second) 
{ 
    // If we reach at the end of both strings, we are done 
    if (*first == '\0' && *second == '\0') 
        return true; 
  
    // Make sure to eliminate consecutive '*' 
    if (*first == '*') { 
        while (*(first + 1) == '*') 
            first++; 
    } 
  
    // Make sure that the characters after '*' are present 
    // in second string. This function assumes that the 
    // first string will not contain two consecutive '*' 
    if (*first == '*' && *(first + 1) != '\0'
        && *second == '\0') 
        return false; 
  
    // If the first string contains '?', or current 
    // characters of both strings match 
    if (*first == '?' || *first == *second) 
        return __match(first + 1, second + 1); 
  
    // If there is *, then there are two possibilities 
    // a) We consider current character of second string 
    // b) We ignore current character of second string. 
    if (*first == '*') 
        return __match(first + 1, second) 
               || __match(first, second + 1); 
    return false; 
}

/** @brief Wrapper for __match().
 * @param first First string, which may contain wildcards
 * @param second Second string
 * @return true if they match, false otherwise
 */
bool match(const std::string first, const std::string second) 
{
	if (first == second)
		return true;
	return __match(first.c_str(), second.c_str());
}

static std::vector<std::string> cpp_str;
static std::vector<XMLCh*> xml_str;

static std::string
xmlTranscode(const XMLCh* txt)
{
	char* temp = xercesc::XMLString::transcode(txt);
	std::string str(temp);
	xercesc::XMLString::release(&temp);
	cpp_str.push_back(str);
	return str;
}

static XMLCh*
xmlTranscode(std::string txt)
{
	XMLCh* str = xercesc::XMLString::transcode(txt.c_str());
	xml_str.push_back(str);
	return str;
}

static void
xmlClear()
{
	unsigned int i;
	for (i = 0; i < xml_str.size(); i++) {
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
#define xmlAttribute(elem, att) xmlS(elem->getAttribute(xmlS(att)))

#ifdef xmlHasAttribute
#undef xmlHasAttribute
#endif
#define xmlHasAttribute(elem, att) elem->hasAttribute(xmlS(att))

using namespace xercesc;
using namespace std;

namespace Aqua {
namespace InputOutput {

State::State()
{
	unsigned int i;
	struct lconv* lc;
	char* s;

	// Set the decimal-point character (which is depending on the locale)
	s = setlocale(LC_NUMERIC, NULL);
	if (strcmp(s, "C")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" numeric locale found" << std::endl;
		LOG(L_INFO, msg.str());
		LOG0(L_DEBUG, "\tIt is replaced by \"C\"\n");
		setlocale(LC_NUMERIC, "C");
	}
	lc = localeconv();
	s = lc->decimal_point;
	if (strcmp(s, ".")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" decimal point character found" << std::endl;
		LOG(L_WARNING, msg.str());
		LOG0(L_DEBUG, "\tIt is replaced by \".\"\n");
		lc->decimal_point = ".";
	}
	s = lc->thousands_sep;
	if (strcmp(s, "")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" thousands separator character found"
		    << std::endl;
		LOG(L_WARNING, msg.str());
		LOG0(L_DEBUG, "\tIt is removed\n");
		lc->thousands_sep = "";
	}

	// Start the XML parser
	try {
		XMLPlatformUtils::Initialize();
	} catch (XMLException& e) {
		std::ostringstream msg;
		LOG(L_ERROR, "XML toolkit initialization error.\n");
		msg << "\t" << xmlS(e.getMessage()) << std::endl;
		LOG0(L_DEBUG, msg.str());
		xmlClear();
		throw;
	} catch (...) {
		LOG(L_ERROR, "XML toolkit initialization error.\n");
		LOG0(L_DEBUG, "\tUnhandled exception\n");
		xmlClear();
		throw;
	}

	// Look ofr the first available file place
	i = 0;
	std::ostringstream file_name;
	while (true) {
		file_name.str("");
		file_name << "AQUAgpusph.save." << i << ".xml";
		if (isFile(file_name.str())) {
			// The file already exist, look for another one
			i++;
			continue;
		}
		break;
	}
	_output_file = file_name.str();
}

State::~State()
{
	// Terminate Xerces
	try {
		XMLPlatformUtils::Terminate();
	} catch (xercesc::XMLException& e) {
		std::ostringstream msg;
		LOG(L_ERROR, "XML toolkit exit error.\n");
		msg << "\t" << xmlS(e.getMessage()) << std::endl;
		LOG0(L_DEBUG, msg.str());
		xmlClear();
	} catch (...) {
		LOG(L_ERROR, "XML toolkit exit error.\n");
		LOG0(L_DEBUG, "\tUnhandled exception\n");
		xmlClear();
	}
}

void
State::save(ProblemSetup& sim_data, std::vector<Particles*> savers)
{
	return write(_output_file, sim_data, savers);
}

void
State::load(std::string input_file, ProblemSetup& sim_data)
{
	return parse(input_file, sim_data);
}

const std::string
State::findPath(const std::string& filepath, ProblemSetup& sim_data)
{
	std::vector<std::string> candidates;
	auto fp = std::filesystem::path(filepath);
	if (std::filesystem::exists(fp)) {
		return std::filesystem::canonical(fp).string();
	}
	candidates.push_back((std::filesystem::current_path() / fp).string());
	if (fp.is_relative()) {
		if (sim_data.settings.base_path != "") {
			auto f =
				std::filesystem::path(sim_data.settings.base_path) / fp;
			if (std::filesystem::exists(f)) {
				return std::filesystem::canonical(f).string();
			}
			candidates.push_back(f.string());
		}
		auto f =
			std::filesystem::path(_xml_paths.back()) / fp;
		if (std::filesystem::exists(f)) {
			return std::filesystem::canonical(f).string();
		}
		candidates.push_back(f.string());
	}

	std::ostringstream msg;
	msg << "No such file or directory '" << fp.string() << "'" << std::endl
	    << "The following paths were checked:" << std::endl;
	for (auto candidate : candidates) {
		msg << "  '" << candidate << "'" << std::endl;
	}
	throw std::filesystem::filesystem_error(msg.str(), std::error_code());
}

void
State::parse(std::string filepath, ProblemSetup& sim_data, std::string prefix)
{
	DOMNodeList* nodes = NULL;
	std::ostringstream msg;
	msg << "Parsing the XML file \"" << filepath << "\" with prefix \""
	    << prefix << "\"" << std::endl;
	LOG(L_INFO, msg.str());

	// Try to open as ascii file, just to know if the file already exist
	std::ifstream f(filepath);
	if (!f) {
		LOG(L_ERROR, "File inaccessible!\n");
		throw std::ifstream::failure("File inaccessible!");
	}
	f.close();

	auto fp = std::filesystem::path(filepath);
	_xml_paths.push_back(fp.parent_path().string());

	// Now we can proceed to properly parse the XML file
	XercesDOMParser* parser = new XercesDOMParser();
	parser->setValidationScheme(XercesDOMParser::Val_Never);
	parser->setDoNamespaces(false);
	parser->setDoSchema(false);
	parser->setLoadExternalDTD(false);
	parser->parse(filepath.c_str());
	DOMDocument* doc = parser->getDocument();
	DOMElement* root = doc->getDocumentElement();
	if (!root) {
		LOG(L_ERROR, "Empty XML file\n");
		throw std::runtime_error("Empty XML file");
	}

	// Parse <Include> tags to recursively load linked XML files
	nodes = root->getElementsByTagName(xmlS("Include"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		// By default, the include statements are parsed at the very beginning.
		if (xmlHasAttribute(elem, "when")) {
			if (xmlAttribute(elem, "when").compare("begin"))
				continue;
		}
		std::string included_file = trimCopy(xmlAttribute(elem, "file"));
		try {
			included_file = findPath(included_file, sim_data);
		} catch (std::filesystem::filesystem_error const& e) {
			LOG(L_ERROR,
			    std::string("Failure including XML:\n") + e.what() + "\n");
			throw;
		}
		std::string included_prefix = prefix;
		if (xmlHasAttribute(elem, "prefix")) {
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
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		if (!xmlHasAttribute(elem, "when"))
			continue;
		if (xmlAttribute(elem, "when").compare("end"))
			continue;
		std::string included_file = trimCopy(xmlAttribute(elem, "file"));
		try {
			included_file = findPath(included_file, sim_data);
		} catch (std::filesystem::filesystem_error const& e) {
			LOG(L_ERROR,
			    std::string("Failure including XML:\n") + e.what() + "\n");
			throw;
		}
		std::string included_prefix = prefix;
		if (xmlHasAttribute(elem, "prefix")) {
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
	_xml_paths.pop_back();
	delete parser;
}

void
State::parseSettings(DOMElement* root,
                     ProblemSetup& sim_data,
                     std::string UNUSED_PARAM prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Settings"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		DOMNodeList* s_nodes;

		s_nodes = elem->getElementsByTagName(xmlS("SaveOnFail"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!toLowerCopy(xmlAttribute(s_elem, "value")).compare("true")) {
				sim_data.settings.save_on_fail = true;
			} else {
				sim_data.settings.save_on_fail = false;
			}
		}
		s_nodes = elem->getElementsByTagName(xmlS("DebugTools"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			std::vector<std::pair<std::string, ProblemSetup::sphSettings::debug_opts>> attrs = {
				{"events", ProblemSetup::sphSettings::DEBUG_EVENTS},
				{"args", ProblemSetup::sphSettings::DEBUG_ARGS},
				{"sync", ProblemSetup::sphSettings::DEBUG_DEPS_SYNC},
				{"deps", ProblemSetup::sphSettings::DEBUG_SYNC}};
			for (const auto& attr : attrs) {
				const std::string s = attr.first;
				const ProblemSetup::sphSettings::debug_opts o = attr.second;
				if (!toLowerCopy(xmlAttribute(s_elem, s)).compare("true")) {
					sim_data.settings.debug_tools |= o;
				} else if (!toLowerCopy(xmlAttribute(s_elem, s)).compare("")){
					if(sim_data.settings.debug_tools & o)
						sim_data.settings.debug_tools -= o;
				}
			}
		}

		s_nodes = elem->getElementsByTagName(xmlS("RootPath"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!xmlHasAttribute(s_elem, "path")) {
				LOG(L_ERROR, "RootPath requires a \"path\" attribute\n");
				throw std::runtime_error("Missing attribute");
			}
			std::string p;
			try {
				p = findPath(xmlAttribute(s_elem, "path"), sim_data);
			} catch (std::filesystem::filesystem_error const& e) {
				std::ostringstream msg;
				msg << "Invalid RootPath \""
					<< xmlAttribute(s_elem, "path") << "\"" << std::endl
					<< e.what() << std::endl;
				LOG(L_WARNING, msg.str());
				continue;
			}
			LOG(L_INFO, std::string("RootPath='") + p + "'\n");
			sim_data.settings.base_path = p;
		}

		s_nodes = elem->getElementsByTagName(xmlS("Device"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			unsigned int platform_id =
			    std::stoi(xmlAttribute(s_elem, "platform"));
			unsigned int device_id = std::stoi(xmlAttribute(s_elem, "device"));
			cl_device_type device_type = CL_DEVICE_TYPE_ALL;
			if (xmlHasAttribute(s_elem, "type")) {
				if (!xmlAttribute(s_elem, "type").compare("ALL"))
					device_type = CL_DEVICE_TYPE_ALL;
				else if (!xmlAttribute(s_elem, "type").compare("CPU"))
					device_type = CL_DEVICE_TYPE_CPU;
				else if (!xmlAttribute(s_elem, "type").compare("GPU"))
					device_type = CL_DEVICE_TYPE_GPU;
				else if (!xmlAttribute(s_elem, "type").compare("ACCELERATOR"))
					device_type = CL_DEVICE_TYPE_ACCELERATOR;
				else if (!xmlAttribute(s_elem, "type").compare("DEFAULT"))
					device_type = CL_DEVICE_TYPE_DEFAULT;
				else {
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
			unsigned int addr_bits = 32;
			if (xmlHasAttribute(s_elem, "addr_bits")) {
				addr_bits = std::stoi(xmlAttribute(s_elem, "addr_bits"));
			}
			std::string compile_flags = "";
			if (xmlHasAttribute(s_elem, "compile_flags")) {
				compile_flags = xmlAttribute(s_elem, "compile_flags");
			}
			sim_data.settings.devices.push_back(
				ProblemSetup::sphSettings::device(platform_id,
				                                  device_id,
				                                  device_type,
				                                  addr_bits,
				                                  compile_flags));
			std::vector<std::pair<
				std::string,
				ProblemSetup::sphSettings::device::patch_state>> attrs = {
					{"enable_patch",
					 ProblemSetup::sphSettings::device::patch_state::ENABLED},
					{"disable_patch",
					 ProblemSetup::sphSettings::device::patch_state::DISABLED}
				};
			for (const auto& attr : attrs) {
				const std::string s = attr.first;
				const auto o = attr.second;
				if (xmlHasAttribute(s_elem, s)) {
					for (auto name : split(xmlAttribute(s_elem, s))) {
						sim_data.settings.devices.back().patches[name] = o;
					}
				}
			}

		}
	}
}

void
State::parseVariables(DOMElement* root,
                      ProblemSetup& sim_data,
                      std::string UNUSED_PARAM prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Variables"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Variable"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!xmlHasAttribute(s_elem, "name")) {
				LOG(L_ERROR, "Unnamed variable");
				throw std::runtime_error("Invalid variable name");
			}
			std::string var_name = xmlAttribute(s_elem, "name");
			if (hasPrefix(var_name, "__")) {
				std::ostringstream msg;
				msg << "Invalid variable name \"" << var_name
				    << "\", since prefix \"__\" is reserved" << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("Invalid variable name");
			}
			for (auto suffix : { "_x", "_y", "_z", "_w" }) {
				if (hasSuffix(var_name, suffix)) {
					std::ostringstream msg;
					msg << "Invalid variable name \"" << var_name
					    << "\", since suffix \"" << suffix << "\" is reserved"
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Invalid variable name");
				}
			}

			if (xmlAttribute(s_elem, "type").find('*') == std::string::npos) {
				sim_data.variables.registerVariable(
				    var_name,
				    xmlAttribute(s_elem, "type"),
				    "1",
				    xmlAttribute(s_elem, "value"));
			} else {
				sim_data.variables.registerVariable(
				    var_name,
				    xmlAttribute(s_elem, "type"),
				    xmlAttribute(s_elem, "length"),
				    "");
			}
		}
	}
	return;
}

void
State::parseDefinitions(DOMElement* root,
                        ProblemSetup& sim_data,
                        std::string UNUSED_PARAM prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Definitions"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Define"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!xmlHasAttribute(s_elem, "name")) {
				LOG(L_ERROR, "Found a definition without name\n");
				throw std::runtime_error(
				    "Name shall be specified for definitions");
			}
			if (!xmlHasAttribute(s_elem, "value")) {
				sim_data.definitions.define(
				    xmlAttribute(s_elem, "name"), "", false);
				continue;
			}

			bool evaluate = false;
			if (!toLowerCopy(xmlAttribute(s_elem, "evaluate"))
			         .compare("true") ||
			    !toLowerCopy(xmlAttribute(s_elem, "evaluate")).compare("yes")) {
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
static std::vector<unsigned int> _tool_places;

/** @brief Helper function to get a list of tool placements from a list of names
 * @param list List of tools, separated by commas
 * @param prefix prefix to become inserted at the beggining of the name of each
 * tool of the list
 * @return The positions of the tools
 * @warning This methos is not thread safe
 */
static std::vector<unsigned int>
_toolsList(std::string list, ProblemSetup& sim_data, std::string prefix)
{
	_tool_places.clear();

	std::istringstream f(list);
	std::string s;
	while (getline(f, s, ',')) {
		std::ostringstream toolname;
		toolname << prefix << s;
		// Look for the tool in the already defined ones
		unsigned int place;
		for (place = 0; place < sim_data.tools.size(); place++) {
			if (!toolname.str().compare(
			        sim_data.tools.at(place)->get("name"))) {
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
static std::vector<unsigned int>
_toolsName(std::string name, ProblemSetup& sim_data, std::string prefix)
{
	_tool_places.clear();
	std::ostringstream toolname;
	toolname << prefix << name;

	// Look for the patterns in the already defined tool names
	unsigned int place;
	for (place = 0; place < sim_data.tools.size(); place++) {
		if (match(toolname.str(), sim_data.tools.at(place)->get("name"))) {
			_tool_places.push_back(place);
		}
	}

	return _tool_places;
}

/** @brief Helper function to get and XML attribute for a tool
 *
 * This function will set the default value if the attribute is not found
 * @param tool The tool to edit
 * @param elem The XML element
 * @param attr The XML element attribute
 * @param def_val The default string to be set if it cannot be found
 */
static void
_toolAttr(ProblemSetup::sphTool* tool,
          DOMElement* elem,
          std::string attr,
          std::string def_val)
{
	if (!xmlHasAttribute(elem, attr.c_str())) {
		tool->set(attr, def_val);
		return;
	}
	tool->set(attr, xmlAttribute(elem, attr.c_str()));
}

/** @brief Helper function to get and XML attribute for a tool
 *
 * This function will throw an exception if the attribute is missing
 * @param tool The tool to edit
 * @param elem The XML element
 * @param attr The XML element attribute
 * @throws std::runtime_error If @p attr is not an attribute of the @p elem XML
 * element
 */
static void
_toolAttr(ProblemSetup::sphTool* tool, DOMElement* elem, std::string attr)
{
	if (!xmlHasAttribute(elem, attr.c_str())) {
		std::ostringstream msg;
		msg << "Tool \"" << tool->get("name")
		    << "\" requires the missing attribute \"" << attr << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Missing attribute");
	}
	tool->set(attr, xmlAttribute(elem, attr.c_str()));
}

void
State::parseTools(DOMElement* root, ProblemSetup& sim_data, std::string prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Tools"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Tool"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!xmlHasAttribute(s_elem, "name")) {
				LOG(L_ERROR, "Found a tool without name\n");
				throw std::runtime_error("Name shall be defined for tools");
			}
			if (!xmlHasAttribute(s_elem, "type")) {
				LOG(L_ERROR, "Found a tool without type\n");
				throw std::runtime_error("Type shall be defined for tools");
			}

			// Create the tool
			ProblemSetup::sphTool* tool = new ProblemSetup::sphTool();

			// Set the name (with prefix)
			std::ostringstream name;
			name << prefix << xmlAttribute(s_elem, "name");
			tool->set("name", name.str());

			tool->set("type", xmlAttribute(s_elem, "type"));
			if (xmlHasAttribute(s_elem, "once")) {
				tool->set("once", toLowerCopy(xmlAttribute(s_elem, "once")));
			} else {
				tool->set("once", "false");
			}

			// Check if the conditions to add the tool are fulfilled
			if (xmlHasAttribute(s_elem, "ifdef")) {
				if (!sim_data.definitions.isDefined(
				        xmlAttribute(s_elem, "ifdef"))) {
					std::ostringstream msg;
					msg << "Ignoring the tool \"" << tool->get("name")
					    << "\" because \"" << xmlAttribute(s_elem, "ifdef")
					    << "\" has not been defined." << std::endl;
					LOG(L_WARNING, msg.str());
					delete tool;
					continue;
				}
			} else if (xmlHasAttribute(s_elem, "ifndef")) {
				if (sim_data.definitions.isDefined(
				        xmlAttribute(s_elem, "ifndef"))) {
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
			if (!xmlHasAttribute(s_elem, "action") ||
			    !xmlAttribute(s_elem, "action").compare("add")) {
				sim_data.tools.push_back(tool);
			} else if (!xmlAttribute(s_elem, "action").compare("insert") ||
			           !xmlAttribute(s_elem, "action").compare("try_insert")) {
				unsigned int place;
				std::vector<unsigned int> places;
				std::vector<unsigned int> all_places;
				std::vector<unsigned int>::iterator it;

				bool try_insert =
				    !xmlAttribute(s_elem, "action").compare("try_insert");

				if (xmlHasAttribute(s_elem, "at")) {
					places.push_back(std::stoi(xmlAttribute(s_elem, "at")));
				} else if (xmlHasAttribute(s_elem, "before") ||
				           xmlHasAttribute(s_elem, "before_prefix")) {
					std::string att_str, att_prefix;
					if (xmlHasAttribute(s_elem, "before")) {
						att_str = xmlAttribute(s_elem, "before");
						att_prefix = "";
					} else {
						att_str = xmlAttribute(s_elem, "before_prefix");
						att_prefix = prefix;
					}
					if (att_str.find(',') != std::string::npos) {
						// It is a list of names. We must get all the matching
						// places and select the most convenient one
						all_places = _toolsList(att_str, sim_data, att_prefix);
						if (!all_places.size()) {
							std::ostringstream msg;
							msg << "The tool \"" << tool->get("name")
							    << "\" must be inserted before \"" << att_str
							    << "\", but such tools cannot be found."
							    << std::endl;
							delete tool;
							if (try_insert) {
								LOG(L_WARNING, msg.str());
								continue;
							} else {
								LOG(L_ERROR, msg.str());
								throw std::runtime_error(
								    "Refered tool cannot be found");
							}
						}
						// Get just the first one
						it = std::min_element(all_places.begin(),
						                      all_places.end());
						places.push_back(*it);
					} else {
						// We can treat the string as a wildcard. Right now we
						// wanna get all the places matching the pattern
						all_places = _toolsName(att_str, sim_data, att_prefix);
						if (!all_places.size()) {
							std::ostringstream msg;
							msg << "The tool \"" << tool->get("name")
							    << "\" must be inserted before \"" << att_str
							    << "\", but such tools cannot be found."
							    << std::endl;
							delete tool;
							if (try_insert) {
								LOG(L_WARNING, msg.str());
								continue;
							} else {
								LOG(L_ERROR, msg.str());
								throw std::runtime_error(
								    "Refered tool cannot be found");
							}
						}
						// Deep copy the places
						for (auto place : all_places)
							places.push_back(place);
					}
				} else if (xmlHasAttribute(s_elem, "after") ||
				           xmlHasAttribute(s_elem, "after_prefix")) {
					std::string att_str, att_prefix;
					if (xmlHasAttribute(s_elem, "after")) {
						att_str = xmlAttribute(s_elem, "after");
						att_prefix = "";
					} else {
						att_str = xmlAttribute(s_elem, "after_prefix");
						att_prefix = prefix;
					}
					if (att_str.find(',') != std::string::npos) {
						// It is a list of names. We must get all the matching
						// places and select the most convenient one
						all_places = _toolsList(att_str, sim_data, att_prefix);
						if (!all_places.size()) {
							std::ostringstream msg;
							msg << "The tool \"" << tool->get("name")
							    << "\" must be inserted after \"" << att_str
							    << "\", but such tools cannot be found."
							    << std::endl;
							delete tool;
							if (try_insert) {
								LOG(L_WARNING, msg.str());
								continue;
							} else {
								LOG(L_ERROR, msg.str());
								throw std::runtime_error(
								    "Refered tool cannot be found");
							}
						}
						// Get just the last one (and insert after that)
						it = std::max_element(all_places.begin(),
						                      all_places.end());
						places.push_back(*it + 1);
					} else {
						// We can treat the string as a wildcard. Right now we
						// wanna get all the places matching the pattern
						all_places = _toolsName(att_str, sim_data, att_prefix);
						if (!all_places.size()) {
							std::ostringstream msg;
							msg << "The tool \"" << tool->get("name")
							    << "\" must be inserted after \"" << att_str
							    << "\", but such tools cannot be found."
							    << std::endl;
							delete tool;
							if (try_insert) {
								LOG(L_WARNING, msg.str());
								continue;
							} else {
								LOG(L_ERROR, msg.str());
								throw std::runtime_error(
								    "Refered tool cannot be found");
							}
						}
						// Deep copy the places (adding 1 to insert after that)
						for (auto place : all_places)
							places.push_back(place + 1);
					}
				} else {
					std::ostringstream msg;
					msg << "Missed the place where the tool \""
					    << tool->get("name") << "\" should be inserted."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					LOG0(L_DEBUG,
					     "Please set one of the following attributes:\n");
					LOG0(L_DEBUG, "\t\"in\"\n");
					LOG0(L_DEBUG, "\t\"before\"\n");
					LOG0(L_DEBUG, "\t\"after\"\n");
					LOG0(L_DEBUG, "\t\"before_prefix\"\n");
					LOG0(L_DEBUG, "\t\"after_prefix\"\n");
					throw std::runtime_error(
					    "A way to insert the tool shall be specified");
				}
				// We cannot directly insert the tools, because the places would
				// change meanwhile, so better backward adding them
				for (place = places.size(); place > 0; place--) {
					sim_data.tools.insert(
					    sim_data.tools.begin() + places.at(place - 1), tool);
				}
			} else if (!xmlAttribute(s_elem, "action").compare("remove") ||
			           !xmlAttribute(s_elem, "action").compare("try_remove")) {
				bool try_remove =
				    !xmlAttribute(s_elem, "action").compare("try_remove");
				unsigned int place;
				std::vector<unsigned int> places;
				// Get the places of the tools selected
				places = _toolsName(tool->get("name"), sim_data, prefix);
				if (!places.size()) {
					std::ostringstream msg;
					msg << "Failure removing the tool \"" << tool->get("name")
					    << "\". No such tool." << std::endl;
					delete tool;
					if (try_remove) {
						LOG(L_WARNING, msg.str());
						continue;
					} else {
						LOG(L_ERROR, msg.str());
						throw std::runtime_error("No such tool");
					}
				}
				// Delete the new tool (which is useless)
				delete tool;
				// Delete the tools in backward order
				for (place = places.size(); place > 0; place--) {
					if (sim_data.toolInstances(
					        sim_data.tools.at(places.at(place - 1))) == 1) {
						// This is the last instance
						delete sim_data.tools.at(places.at(place - 1));
					}
					// Drop the tool from the list
					sim_data.tools.erase(sim_data.tools.begin() +
					                     places.at(place - 1));
				}
				continue;
			} else if (!xmlAttribute(s_elem, "action").compare("replace") ||
			           !xmlAttribute(s_elem, "action").compare("try_replace")) {
				bool try_replace =
				    !xmlAttribute(s_elem, "action").compare("try_replace");
				unsigned int place;
				std::vector<unsigned int> places;
				// Get the places
				places = _toolsName(tool->get("name"), sim_data, prefix);
				if (!places.size()) {
					std::ostringstream msg;
					msg << "Failure replacing the tool \"" << tool->get("name")
					    << "\". No such tool." << std::endl;
					delete tool;
					if (try_replace) {
						LOG(L_WARNING, msg.str());
						continue;
					} else {
						LOG(L_ERROR, msg.str());
						throw std::runtime_error("No such tool");
					}
				}
				// Replace the tools
				for (place = 0; place < places.size(); place++) {
					if (sim_data.toolInstances(
					        sim_data.tools.at(places.at(place))) == 1) {
						// This is the last instance
						delete sim_data.tools.at(places.at(place));
					}
					// Set the new tool
					sim_data.tools.at(places.at(place)) = tool;
				}
			} else {
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
			if (!xmlAttribute(s_elem, "type").compare("kernel")) {
				if (!xmlHasAttribute(s_elem, "path")) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" requires the missing attribute \"path\"."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing attribute");
				}
				std::string p;
				try {
					p = findPath(xmlAttribute(s_elem, "path"), sim_data);
				} catch (std::filesystem::filesystem_error const& e) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" asks for the inaccessible file \""
					    << xmlAttribute(s_elem, "path") << "\"" << std::endl
					    << e.what() << std::endl;
					LOG(L_ERROR, msg.str());
					throw;
				}
				tool->set("path", p);
				_toolAttr(tool, s_elem, "entry_point", "entry");
				_toolAttr(tool, s_elem, "n", "");
		std::string included_file = trimCopy(xmlAttribute(elem, "file"));

			} else if (!xmlAttribute(s_elem, "type").compare("copy")) {
				for (auto attr : { "in", "out" })
					_toolAttr(tool, s_elem, attr);
			} else if (!xmlAttribute(s_elem, "type").compare("python")) {
				if (!xmlHasAttribute(s_elem, "path")) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" requires the missing attribute \"path\"."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing attribute");
				}
				std::string p;
				try {
					p = findPath(xmlAttribute(s_elem, "path"), sim_data);
				} catch (std::filesystem::filesystem_error const& e) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" asks for the inaccessible file \""
					    << xmlAttribute(s_elem, "path") << "\"" << std::endl
					    << e.what() << std::endl;
					LOG(L_ERROR, msg.str());
					throw;
				}
				tool->set("path", p);
			} else if (!xmlAttribute(s_elem, "type").compare("set")) {
				for (auto attr : { "in", "value" })
					_toolAttr(tool, s_elem, attr);
			} else if (!xmlAttribute(s_elem, "type").compare("set_scalar")) {
				for (auto attr : { "in", "value" })
					_toolAttr(tool, s_elem, attr);
			} else if (!xmlAttribute(s_elem, "type").compare("reduction")) {
				for (auto attr : { "in", "out", "null" })
					_toolAttr(tool, s_elem, attr);
				if (!xmlS(s_elem->getTextContent()).compare("")) {
					std::ostringstream msg;
					msg << "No operation specified for the reduction \""
					    << tool->get("name") << "\"." << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing reduction operation");
				}
				tool->set("operation", xmlS(s_elem->getTextContent()));
			} else if (!xmlAttribute(s_elem, "type").compare("link-list")) {
				_toolAttr(tool, s_elem, "in", "r");
				_toolAttr(tool, s_elem, "min", "r_min");
				_toolAttr(tool, s_elem, "max", "r_max");
				_toolAttr(tool, s_elem, "ihoc", "ihoc");
				_toolAttr(tool, s_elem, "icell", "icell");
				_toolAttr(tool, s_elem, "n_cells", "n_cells");
				_toolAttr(tool, s_elem, "perm", "id_unsorted");
				_toolAttr(tool, s_elem, "inv_perm", "id_sorted");
				_toolAttr(tool, s_elem, "recompute_grid", "true");
				_toolAttr(tool, s_elem, "sorter", "radix-sort");
			} else if (!xmlAttribute(s_elem, "type").compare("radix-sort")) {
				for (auto attr : { "in", "perm", "inv_perm" })
					_toolAttr(tool, s_elem, attr);
			} else if (!xmlAttribute(s_elem, "type").compare("sort")) {
				for (auto attr : { "in", "perm", "inv_perm" })
					_toolAttr(tool, s_elem, attr);
			} else if (!xmlAttribute(s_elem, "type").compare("assert")) {
				_toolAttr(tool, s_elem, "condition");
			} else if (!xmlAttribute(s_elem, "type").compare("if")) {
				_toolAttr(tool, s_elem, "condition");
			} else if (!xmlAttribute(s_elem, "type").compare("while")) {
				_toolAttr(tool, s_elem, "condition");
			} else if (!xmlAttribute(s_elem, "type").compare("endif")) {
			} else if (!xmlAttribute(s_elem, "type").compare("end")) {
			}
#ifdef HAVE_MPI
			else if (!xmlAttribute(s_elem, "type").compare("mpi-sync")) {
				for (auto attr : { "mask", "fields" })
					_toolAttr(tool, s_elem, attr);
				_toolAttr(tool, s_elem, "processes", "");
			}
#endif
			else if (!xmlAttribute(s_elem, "type").compare("installable")) {
				if (!xmlHasAttribute(s_elem, "path")) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" requires the missing attribute \"path\"."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing attribute");
				}
				std::string p;
				try {
					p = findPath(xmlAttribute(s_elem, "path"), sim_data);
				} catch (std::filesystem::filesystem_error const& e) {
					std::ostringstream msg;
					msg << "Tool \"" << tool->get("name")
					    << "\" asks for the inaccessible file \""
					    << xmlAttribute(s_elem, "path") << "\"" << std::endl
					    << e.what() << std::endl;
					LOG(L_ERROR, msg.str());
					throw;
				}
				tool->set("path", p);
			} else if (!xmlAttribute(s_elem, "type").compare("dummy")) {
				// Without options
			}
			// Reports directly dealt as tools. Suited for the non-linear
			// simulations editor
			else if (!xmlAttribute(s_elem, "type").compare("report_screen")) {
				_toolAttr(tool, s_elem, "fields");
				_toolAttr(tool, s_elem, "bold", "false");
				_toolAttr(tool, s_elem, "color", "white");
			} else if (!xmlAttribute(s_elem, "type").compare("report_file")) {
				_toolAttr(tool, s_elem, "fields");
				_toolAttr(tool, s_elem, "path");
			} else if (!xmlAttribute(s_elem, "type")
			                .compare("report_particles")) {
				_toolAttr(tool, s_elem, "fields");
				_toolAttr(tool, s_elem, "path");
				_toolAttr(tool, s_elem, "set");
				_toolAttr(tool, s_elem, "ipf", "1");
				_toolAttr(tool, s_elem, "fps", "0.0");
			} else if (!xmlAttribute(s_elem, "type")
			                .compare("report_dump")) {
				_toolAttr(tool, s_elem, "fields");
				_toolAttr(tool, s_elem, "path");
				_toolAttr(tool, s_elem, "binary", "false");
			} else if (!xmlAttribute(s_elem, "type")
			                .compare("report_performance")) {
				_toolAttr(tool, s_elem, "bold", "false");
				_toolAttr(tool, s_elem, "color", "white");
				_toolAttr(tool, s_elem, "path", "");
			} else {
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
				LOG0(L_DEBUG, "\t\tsort\n");
				LOG0(L_DEBUG, "\t\tassert\n");
				LOG0(L_DEBUG, "\t\tif\n");
				LOG0(L_DEBUG, "\t\twhile\n");
				LOG0(L_DEBUG, "\t\tend\n");
				LOG0(L_DEBUG, "\t\tinstallable\n");
				LOG0(L_DEBUG, "\t\tdummy\n");
				LOG0(L_DEBUG, "\t\treport_screen\n");
				LOG0(L_DEBUG, "\t\treport_file\n");
				LOG0(L_DEBUG, "\t\treport_particles\n");
				LOG0(L_DEBUG, "\t\treport_performance\n");
				throw std::runtime_error("Unknown tool type");
			}
		}
	}
}

void
State::parseTiming(DOMElement* root,
                   ProblemSetup& sim_data,
                   std::string UNUSED_PARAM prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Timing"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		// Get options
		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Option"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			std::string name = xmlAttribute(s_elem, "name");
			if (!name.compare("End") || !name.compare("SimulationStop")) {
				std::string type = xmlAttribute(s_elem, "type");
				if (!type.compare("Time") || !type.compare("T")) {
					sim_data.time_opts.sim_end_mode =
					    sim_data.time_opts.sim_end_mode | __TIME_MODE__;
					sim_data.time_opts.sim_end_time =
					    std::stof(xmlAttribute(s_elem, "value"));
				} else if (!type.compare("Steps") || !type.compare("S")) {
					sim_data.time_opts.sim_end_mode =
					    sim_data.time_opts.sim_end_mode | __ITER_MODE__;
					sim_data.time_opts.sim_end_step =
					    std::stoi(xmlAttribute(s_elem, "value"));
				} else if (!type.compare("Frames") || !type.compare("F")) {
					sim_data.time_opts.sim_end_mode =
					    sim_data.time_opts.sim_end_mode | __FRAME_MODE__;
					sim_data.time_opts.sim_end_frame =
					    std::stoi(xmlAttribute(s_elem, "value"));
				} else {
					std::ostringstream msg;
					msg << "Unknown simulation stop criteria \"" << type
					    << "\"." << std::endl;
					LOG(L_ERROR, msg.str());
					LOG0(L_DEBUG, "\tThe valid options are:\n");
					LOG0(L_DEBUG, "\t\tTime\n");
					LOG0(L_DEBUG, "\t\tSteps\n");
					LOG0(L_DEBUG, "\t\tFrames\n");
					throw std::runtime_error(
					    "Unknown simulation stop criteria");
				}
			} else if (!name.compare("Output")) {
				std::string type = xmlAttribute(s_elem, "type");
				if (!type.compare("No")) {
					sim_data.time_opts.output_mode = __NO_OUTPUT_MODE__;
				} else if (!type.compare("FPS")) {
					sim_data.time_opts.output_mode =
					    sim_data.time_opts.output_mode | __FPS_MODE__;
					sim_data.time_opts.output_fps =
					    std::stof(xmlAttribute(s_elem, "value"));
				} else if (!type.compare("IPF")) {
					sim_data.time_opts.output_mode =
					    sim_data.time_opts.output_mode | __IPF_MODE__;
					sim_data.time_opts.output_ipf =
					    std::stoi(xmlAttribute(s_elem, "value"));
				} else {
					std::ostringstream msg;
					msg << "Unknow output file print criteria \"" << type
					    << "\"." << std::endl;
					LOG(L_ERROR, msg.str());
					LOG0(L_DEBUG, "\tThe valid options are:\n");
					LOG0(L_DEBUG, "\t\tNo\n");
					LOG0(L_DEBUG, "\t\tFPS\n");
					LOG0(L_DEBUG, "\t\tIPF\n");
					throw std::runtime_error(
					    "Unknow output file print criteria");
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

void
State::parseSets(DOMElement* root,
                 ProblemSetup& sim_data,
                 std::string UNUSED_PARAM prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("ParticlesSet"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

		ProblemSetup::sphParticlesSet* set =
		    new ProblemSetup::sphParticlesSet();
		// Now the number of particles can be let unknown, and will be
		// determined later, from the input file. See FileManager::load()
		if (xmlHasAttribute(elem, "n")) {
			unsigned long long n = std::stoll(xmlAttribute(elem, "n"));
			try {
				set->n(narrow_cast<size_t>(n));
			} catch(std::out_of_range&) {
				LOG(L_ERROR,
				    std::string("The particles set ") +
				    std::to_string(sim_data.sets.size()) +
				    " has the attribute n=" + xmlAttribute(elem, "n") +
				    ", which overflows the type size_t\n");
				throw;
			}
		}

		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Scalar"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);

			std::string name = xmlAttribute(s_elem, "name");
			std::string value = xmlAttribute(s_elem, "value");
			set->addScalar(name, value);
		}

		s_nodes = elem->getElementsByTagName(xmlS("Load"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			std::string path = xmlAttribute(s_elem, "file");
			std::string format = xmlAttribute(s_elem, "format");
			std::string fields = xmlAttribute(s_elem, "fields");
			set->input(path, format, fields);
		}

		s_nodes = elem->getElementsByTagName(xmlS("Save"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
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

void
State::parseReports(DOMElement* root,
                    ProblemSetup& sim_data,
                    std::string prefix)
{
	DOMNodeList* nodes = root->getElementsByTagName(xmlS("Reports"));
	for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
		DOMNode* node = nodes->item(i);
		if (node->getNodeType() != DOMNode::ELEMENT_NODE)
			continue;
		DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);
		DOMNodeList* s_nodes = elem->getElementsByTagName(xmlS("Report"));
		for (XMLSize_t j = 0; j < s_nodes->getLength(); j++) {
			DOMNode* s_node = s_nodes->item(j);
			if (s_node->getNodeType() != DOMNode::ELEMENT_NODE)
				continue;
			DOMElement* s_elem = dynamic_cast<xercesc::DOMElement*>(s_node);
			if (!xmlHasAttribute(s_elem, "name")) {
				LOG(L_ERROR, "Found a report without name\n");
				throw std::runtime_error("Missing report name");
			}
			if (!xmlHasAttribute(s_elem, "type")) {
				LOG(L_ERROR, "Found a report without type\n");
				throw std::runtime_error("Missing report type");
			}

			// Create the report
			ProblemSetup::sphTool* report = new ProblemSetup::sphTool();

			// Set the name (with prefix)
			std::ostringstream name;
			name << prefix << xmlAttribute(s_elem, "name");
			report->set("name", name.str());

			report->set("type", xmlAttribute(s_elem, "type"));
			sim_data.reports.push_back(report);

			// Configure the report
			if (!xmlAttribute(s_elem, "type").compare("screen")) {
				if (!xmlHasAttribute(s_elem, "fields")) {
					LOG(L_ERROR, "Found a \"screen\" report without fields\n");
					throw std::runtime_error("Missing report fields");
				}
				report->set("fields", xmlAttribute(s_elem, "fields"));
				if (xmlHasAttribute(s_elem, "bold")) {
					report->set("bold", xmlAttribute(s_elem, "bold"));
				} else {
					report->set("bold", "false");
				}
				if (xmlHasAttribute(s_elem, "color")) {
					report->set("color", xmlAttribute(s_elem, "color"));
				} else {
					report->set("color", "white");
				}
			} else if (!xmlAttribute(s_elem, "type").compare("file")) {
				if (!xmlHasAttribute(s_elem, "fields")) {
					LOG(L_ERROR, "Found a \"file\" report without fields\n");
					throw std::runtime_error("Missing report fields");
				}
				report->set("fields", xmlAttribute(s_elem, "fields"));
				if (!xmlHasAttribute(s_elem, "path")) {
					std::ostringstream msg;
					msg << "Report \"" << report->get("name")
					    << "\" is of type \"file\", but the output \"path\" is "
					       "not defined."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing report file path");
				}
				report->set("path", xmlAttribute(s_elem, "path"));
			} else if (!xmlAttribute(s_elem, "type").compare("particles")) {
				if (!xmlHasAttribute(s_elem, "fields")) {
					LOG(L_ERROR,
					    "Found a \"particles\" report without fields\n");
					throw std::runtime_error("Missing report fields");
				}
				report->set("fields", xmlAttribute(s_elem, "fields"));
				if (!xmlHasAttribute(s_elem, "path")) {
					std::ostringstream msg;
					msg << "Report \"" << report->get("name")
					    << "\" is of type \"particles\", but the output "
					       "\"path\" is not defined."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing report file path");
				}
				report->set("path", xmlAttribute(s_elem, "path"));

				if (!xmlHasAttribute(s_elem, "set")) {
					std::ostringstream msg;
					msg << "Report \"" << report->get("name")
					    << "\" is of type \"particles\", but the output "
					       "\"set\" is not defined."
					    << std::endl;
					LOG(L_ERROR, msg.str());
					throw std::runtime_error("Missing report particles set");
				}
				report->set("set", xmlAttribute(s_elem, "set"));

				if (!xmlHasAttribute(s_elem, "ipf")) {
					report->set("ipf", "1");
				} else {
					report->set("ipf", xmlAttribute(s_elem, "ipf"));
				}
				if (!xmlHasAttribute(s_elem, "fps")) {
					report->set("fps", "0.0");
				} else {
					report->set("fps", xmlAttribute(s_elem, "fps"));
				}
			} else if (!xmlAttribute(s_elem, "type").compare("performance")) {
				if (xmlHasAttribute(s_elem, "bold")) {
					report->set("bold", xmlAttribute(s_elem, "bold"));
				} else {
					report->set("bold", "false");
				}
				if (xmlHasAttribute(s_elem, "color")) {
					report->set("color", xmlAttribute(s_elem, "color"));
				} else {
					report->set("color", "white");
				}
				if (xmlHasAttribute(s_elem, "path")) {
					report->set("path", xmlAttribute(s_elem, "path"));
				} else {
					report->set("path", "");
				}
			} else {
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

void
State::write(std::string filepath,
             ProblemSetup& sim_data,
             std::vector<Particles*> savers)
{
	DOMImplementation* impl;
	std::ostringstream msg;
	msg << "Writing \"" << filepath << "\" SPH state file..." << std::endl;
	LOG(L_INFO, msg.str());

	impl = DOMImplementationRegistry::getDOMImplementation(xmlS("Range"));
	DOMDocument* doc = impl->createDocument(NULL, xmlS("sphInput"), NULL);
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
	DOMLSSerializer* saver = ((DOMImplementationLS*)impl)->createLSSerializer();

	if (saver->getDomConfig()->canSetParameter(
	        XMLUni::fgDOMWRTFormatPrettyPrint, true))
		saver->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint,
		                                    true);
	saver->setNewLine(xmlS("\r\n"));

	XMLFormatTarget* target = new LocalFileFormatTarget(filepath.c_str());
	// XMLFormatTarget *target = new StdOutFormatTarget();
	DOMLSOutput* output = ((DOMImplementationLS*)impl)->createLSOutput();
	output->setByteStream(target);
	output->setEncoding(xmlS("UTF-8"));

	try {
		saver->write(doc, output);
	} catch (XMLException& e) {
		LOG(L_ERROR, "XML toolkit writing error.\n");
		msg.str("");
		msg << "\t" << xmlS(e.getMessage()) << std::endl;
		LOG0(L_DEBUG, msg.str());
		xmlClear();
		throw std::runtime_error("Failure writing XML file");
	} catch (DOMException& e) {
		LOG(L_ERROR, "XML DOM writing error.\n");
		msg.str("");
		msg << "\t" << xmlS(e.getMessage()) << std::endl;
		LOG0(L_DEBUG, msg.str());
		xmlClear();
		throw std::runtime_error("XML DOM error while writing XML");
	} catch (...) {
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

	msg.str("");
	msg << "Wrote \"" << filepath << "\" SPH state file..." << std::endl;
	LOG(L_INFO, msg.str());
}

void
State::writeSettings(xercesc::DOMDocument* doc,
                     xercesc::DOMElement* root,
                     ProblemSetup& sim_data)
{
	DOMElement *elem, *s_elem;
	std::ostringstream att;

	elem = doc->createElement(xmlS("Settings"));
	root->appendChild(elem);

	s_elem = doc->createElement(xmlS("SaveOnFail"));
	if (sim_data.settings.save_on_fail)
		s_elem->setAttribute(xmlS("value"), xmlS("true"));
	else
		s_elem->setAttribute(xmlS("value"), xmlS("false"));
	elem->appendChild(s_elem);

	s_elem = doc->createElement(xmlS("DebugTools"));
	std::vector<std::pair<std::string, ProblemSetup::sphSettings::debug_opts>> attrs = {
		{"events", ProblemSetup::sphSettings::DEBUG_EVENTS},
		{"args", ProblemSetup::sphSettings::DEBUG_ARGS},
		{"sync", ProblemSetup::sphSettings::DEBUG_DEPS_SYNC},
		{"deps", ProblemSetup::sphSettings::DEBUG_SYNC}};
	for (const auto& attr : attrs) {
		const std::string s = attr.first;
		const ProblemSetup::sphSettings::debug_opts o = attr.second;
		if (sim_data.settings.debug_tools & o)
			s_elem->setAttribute(xmlS(s), xmlS("true"));
		else
			s_elem->setAttribute(xmlS(s), xmlS("false"));
	}
	elem->appendChild(s_elem);

	s_elem = doc->createElement(xmlS("RootPath"));
	s_elem->setAttribute(xmlS("path"), xmlS(sim_data.settings.base_path));
	elem->appendChild(s_elem);

	for (auto device : sim_data.settings.devices) {
		s_elem = doc->createElement(xmlS("Device"));
		att.str("");
		att << device.platform_id;
		s_elem->setAttribute(xmlS("platform"), xmlS(att.str()));
		att.str("");
		att << device.device_id;
		s_elem->setAttribute(xmlS("device"), xmlS(att.str()));
		att.str("");
		switch (device.device_type) {
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
		att.str("");
		att << device.addr_bits;
		s_elem->setAttribute(xmlS("addr_bits"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}
}

void
State::writeVariables(xercesc::DOMDocument* doc,
                      xercesc::DOMElement* root,
                      ProblemSetup UNUSED_PARAM &sim_data)
{
	DOMElement *elem, *s_elem;
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();

	elem = doc->createElement(xmlS("Variables"));
	root->appendChild(elem);

	std::vector<Variable*> vars = C->variables()->getAll();
	for (auto var : vars) {
		std::string var_name = var->name();
		if (hasPrefix(var_name, "__"))
			continue;

		s_elem = doc->createElement(xmlS("Variable"));
		elem->appendChild(s_elem);

		s_elem->setAttribute(xmlS("name"), xmlS(var_name));
		std::string type = var->type();
		s_elem->setAttribute(xmlS("type"), xmlS(type));

		// Array variable
		if (type.find("*") != std::string::npos) {
			size_t length = ((ArrayVariable*)var)->size() /
			                Variables::typeToBytes(var->type());
			std::ostringstream length_txt;
			length_txt << length;
			s_elem->setAttribute(xmlS("length"), xmlS(length_txt.str()));
			continue;
		}
		// Scalar variable
		std::string value_txt = var->asString();
		if (value_txt.at(0) == '(') {
			value_txt.at(0) = ' ';
		}
		if (value_txt.back() == ')') {
			value_txt.back() = ' ';
		}
		s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
	}
}

void
State::writeDefinitions(xercesc::DOMDocument* doc,
                        xercesc::DOMElement* root,
                        ProblemSetup& sim_data)
{
	unsigned int i;
	DOMElement *elem, *s_elem;

	elem = doc->createElement(xmlS("Definitions"));
	root->appendChild(elem);

	ProblemSetup::sphDefinitions defs = sim_data.definitions;

	for (i = 0; i < defs.names.size(); i++) {
		s_elem = doc->createElement(xmlS("Define"));
		elem->appendChild(s_elem);

		s_elem->setAttribute(xmlS("name"), xmlS(defs.names.at(i)));
		s_elem->setAttribute(xmlS("value"), xmlS(defs.values.at(i)));
		if (defs.evaluations.at(i)) {
			s_elem->setAttribute(xmlS("evaluate"), xmlS("true"));
		} else {
			s_elem->setAttribute(xmlS("evaluate"), xmlS("false"));
		}
	}
}

void
State::writeTools(xercesc::DOMDocument* doc,
                  xercesc::DOMElement* root,
                  ProblemSetup& sim_data)
{
	unsigned int i;
	DOMElement *elem, *s_elem;

	elem = doc->createElement(xmlS("Tools"));
	root->appendChild(elem);

	std::vector<ProblemSetup::sphTool*> tools = sim_data.tools;

	for (auto tool : tools) {
		s_elem = doc->createElement(xmlS("Tool"));
		elem->appendChild(s_elem);
		for (i = 0; i < tool->n(); i++) {
			std::string name = tool->getName(i);
			std::string value = tool->get(i);
			if (!name.compare("operation")) {
				// The reduction operation is not an attribute, but a text
				s_elem->setTextContent(xmlS(value));
				continue;
			}
			s_elem->setAttribute(xmlS(name), xmlS(value));
		}
	}
}

void
State::writeReports(xercesc::DOMDocument* doc,
                    xercesc::DOMElement* root,
                    ProblemSetup& sim_data)
{
	unsigned int i;
	DOMElement *elem, *s_elem;

	elem = doc->createElement(xmlS("Reports"));
	root->appendChild(elem);

	std::vector<ProblemSetup::sphTool*> reports = sim_data.reports;

	for (auto report : reports) {
		s_elem = doc->createElement(xmlS("Report"));
		elem->appendChild(s_elem);
		for (i = 0; i < report->n(); i++) {
			std::string name = report->getName(i);
			std::string value = report->get(i);
			s_elem->setAttribute(xmlS(name), xmlS(value));
		}
	}
}

void
State::writeTiming(xercesc::DOMDocument* doc,
                   xercesc::DOMElement* root,
                   ProblemSetup& sim_data)
{
	DOMElement *elem, *s_elem;
	std::ostringstream att;

	elem = doc->createElement(xmlS("Timing"));
	root->appendChild(elem);

	if (sim_data.time_opts.sim_end_mode & __TIME_MODE__) {
		s_elem = doc->createElement(xmlS("Option"));
		s_elem->setAttribute(xmlS("name"), xmlS("End"));
		s_elem->setAttribute(xmlS("type"), xmlS("Time"));
		att.str("");
		att << sim_data.time_opts.sim_end_time;
		s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}
	if (sim_data.time_opts.sim_end_mode & __ITER_MODE__) {
		s_elem = doc->createElement(xmlS("Option"));
		s_elem->setAttribute(xmlS("name"), xmlS("End"));
		s_elem->setAttribute(xmlS("type"), xmlS("Steps"));
		att.str("");
		att << sim_data.time_opts.sim_end_step;
		s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}
	if (sim_data.time_opts.sim_end_mode & __FRAME_MODE__) {
		s_elem = doc->createElement(xmlS("Option"));
		s_elem->setAttribute(xmlS("name"), xmlS("End"));
		s_elem->setAttribute(xmlS("type"), xmlS("Frames"));
		att.str("");
		att << sim_data.time_opts.sim_end_frame;
		s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}

	if (sim_data.time_opts.output_mode & __FPS_MODE__) {
		s_elem = doc->createElement(xmlS("Option"));
		s_elem->setAttribute(xmlS("name"), xmlS("Output"));
		s_elem->setAttribute(xmlS("type"), xmlS("FPS"));
		att.str("");
		att << sim_data.time_opts.output_fps;
		s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}
	if (sim_data.time_opts.output_mode & __IPF_MODE__) {
		s_elem = doc->createElement(xmlS("Option"));
		s_elem->setAttribute(xmlS("name"), xmlS("Output"));
		s_elem->setAttribute(xmlS("type"), xmlS("IPF"));
		att.str("");
		att << sim_data.time_opts.output_ipf;
		s_elem->setAttribute(xmlS("value"), xmlS(att.str()));
		elem->appendChild(s_elem);
	}
}

void
State::writeSets(xercesc::DOMDocument* doc,
                 xercesc::DOMElement* root,
                 ProblemSetup& sim_data,
                 std::vector<Particles*> savers)
{
	unsigned int i, j;
	std::ostringstream att;
	DOMElement *elem, *s_elem;

	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();

	for (i = 0; i < sim_data.sets.size(); i++) {
		elem = doc->createElement(xmlS("ParticlesSet"));
		att.str("");
		att << sim_data.sets.at(i)->n();
		elem->setAttribute(xmlS("n"), xmlS(att.str()));
		root->appendChild(elem);

		for (auto name : sim_data.sets.at(i)->scalarNames()) {
			s_elem = doc->createElement(xmlS("Scalar"));
			s_elem->setAttribute(xmlS("name"), xmlS(name));

			ArrayVariable* var = (ArrayVariable*)C->variables()->get(name);
			std::string value_txt = var->asString((size_t)i);
			if (value_txt.at(0) == '(') {
				value_txt.at(0) = ' ';
			}
			if (value_txt.back() == ')') {
				value_txt.back() = ' ';
			}
			s_elem->setAttribute(xmlS("value"), xmlS(value_txt));
			elem->appendChild(s_elem);
		}

		std::ostringstream fields;
		for (j = 0; j < sim_data.sets.at(i)->outputFields().size(); j++) {
			std::string field = sim_data.sets.at(i)->outputFields().at(j);
			fields << field;
			if (j < sim_data.sets.at(i)->outputFields().size() - 1)
				fields << ",";
		}
		s_elem = doc->createElement(xmlS("Load"));
		s_elem->setAttribute(xmlS("file"), xmlS(savers.at(i)->file()));
		s_elem->setAttribute(xmlS("format"),
		                     xmlS(sim_data.sets.at(i)->outputFormat()));
		s_elem->setAttribute(xmlS("fields"), xmlS(fields.str()));
		elem->appendChild(s_elem);

		s_elem = doc->createElement(xmlS("Save"));
		s_elem->setAttribute(xmlS("file"),
		                     xmlS(sim_data.sets.at(i)->outputPath()));
		s_elem->setAttribute(xmlS("format"),
		                     xmlS(sim_data.sets.at(i)->outputFormat()));
		s_elem->setAttribute(xmlS("fields"), xmlS(fields.str()));
		elem->appendChild(s_elem);
	}
}

}
} // namespace
