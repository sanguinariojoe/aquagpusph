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

#ifndef STATE_H_INCLUDED
#define STATE_H_INCLUDED

#include <ProblemSetup.h>
#include <InputOutput/InputOutput.h>
#include <InputOutput/Particles.h>

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
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/XMLUni.hpp>

#include <string>
#include <vector>

namespace Aqua {
namespace InputOutput {

/** \class State State.h InputOutput/State.h
 * @brief Load/Save the XML simulation definition files.
 *
 * In AQUAgpusph the input/output managers are divided in 3 different types:
 *   -# The simulation configuration files manager
 *   -# The report file managers
 *   -# The particles output file managers
 *
 * This class is based in the xerces-c library, to learn more please visit the
 * following web page:
 *
 * http://xerces.apache.org/xerces-c
 *
 * You can find more information about how to create simulation configuration
 * files in the examples provided with the package and in the
 * Aqua::InputOutput::ProblemSetup class documentation.
 *
 * @see Aqua::InputOutput::Report
 * @see Aqua::InputOutput::Particles
 */
class State
{
  public:
	/** @brief Constructor
	 */
	State();

	/// Destructor
	~State();

	/** @brief Save the configuration file.
	 *
	 * The saved XML file can be used to continue the simulation from the last
	 * saved time instant.
	 *
	 * The output XML file will be the first non existing file named
	 * `"AQUAgpusph.save.%d.xml"`, where `"%d"` is an unsigned integer.
	 *
	 * Of course, to can load a saved simulation, the output particles file
	 * should be saved as well.
	 *
	 * @see Aqua::InputOutput::Particles::save()
	 * @param sim_data Simulation data
	 * @param savers Particles savers list
	 */
	void save(ProblemSetup& sim_data, std::vector<Particles*> savers);

	/** @brief Load the simulation XML definition files.
	 *
	 * @param input_file XML file to load
	 * @param sim_data Simulation data
	 */
	void load(std::string input_file, ProblemSetup& sim_data);

  protected:
	/** @brief Parse the XML file
	 * @param filepath file to be parsed.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names
	 * defined in the file.
	 * @return false if all gone right, true otherwise
	 * @note This function is calling himself for each `<Include>` tag found,
	 * conveniently changing \a filepath and \a prefix. If no prefix is
	 * specified in the `<Include>` tag, the same received \a prefix will be
	 * used.
	 */
	void parse(std::string filepath,
	           ProblemSetup& sim_data,
	           std::string prefix = "");

	/** @brief Parse the general settings sections.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphSettings
	 */
	void parseSettings(xercesc::DOMElement* root,
	                   ProblemSetup& sim_data,
	                   std::string prefix = "");

	/** @brief Parse the variables sections.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphVariables
	 */
	void parseVariables(xercesc::DOMElement* root,
	                    ProblemSetup& sim_data,
	                    std::string prefix = "");

	/** @brief Parse the definitions sections.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphDefinitions
	 */
	void parseDefinitions(xercesc::DOMElement* root,
	                      ProblemSetup& sim_data,
	                      std::string prefix = "");

	/** @brief Parse the tools sections.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTool
	 */
	void parseTools(xercesc::DOMElement* root,
	                ProblemSetup& sim_data,
	                std::string prefix = "");

	/** @brief Parse the time control sections.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
	 */
	void parseTiming(xercesc::DOMElement* root,
	                 ProblemSetup& sim_data,
	                 std::string prefix = "");

	/** Look for particles set sections.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphParticlesSet
	 */
	void parseSets(xercesc::DOMElement* root,
	               ProblemSetup& sim_data,
	               std::string prefix = "");

	/** @brief Parse the reports sections.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @param prefix String to be inserted before the variable names.
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTool
	 */
	void parseReports(xercesc::DOMElement* root,
	                  ProblemSetup& sim_data,
	                  std::string prefix = "");

	/** @brief Write the XML file
	 * @param filepath file to be written.
	 * @param sim_data Simulation data
	 * @param savers Particles savers list
	 * @return false if all gone right, true otherwise
	 */
	void write(std::string filepath,
	           ProblemSetup& sim_data,
	           std::vector<Particles*> savers);

	/** @brief Write the settings section.
	 * @param doc XML generated document.
	 * @param root root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphSettings
	 */
	void writeSettings(xercesc::DOMDocument* doc,
	                   xercesc::DOMElement* root,
	                   ProblemSetup& sim_data);

	/** @brief Write the variables section.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphVariables
	 */
	void writeVariables(xercesc::DOMDocument* doc,
	                    xercesc::DOMElement* root,
	                    ProblemSetup& sim_data);

	/** @brief Write the definitions section.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphDefinitions
	 */
	void writeDefinitions(xercesc::DOMDocument* doc,
	                      xercesc::DOMElement* root,
	                      ProblemSetup& sim_data);

	/** @brief Write the tools section.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTool
	 */
	void writeTools(xercesc::DOMDocument* doc,
	                xercesc::DOMElement* root,
	                ProblemSetup& sim_data);

	/** @brief Write the time control section.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
	 */
	void writeTiming(xercesc::DOMDocument* doc,
	                 xercesc::DOMElement* root,
	                 ProblemSetup& sim_data);

	/** @brief Write the particles set sections.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @param savers Particles savers list
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphParticlesSet
	 */
	void writeSets(xercesc::DOMDocument* doc,
	               xercesc::DOMElement* root,
	               ProblemSetup& sim_data,
	               std::vector<Particles*> savers);

	/** @brief Write the reports section.
	 * @param doc XML generated document.
	 * @param root Root XML node.
	 * @param sim_data Simulation data
	 * @return false if all gone right, true otherwise
	 * @see Aqua::InputOutput::ProblemSetup::sphTool
	 */
	void writeReports(xercesc::DOMDocument* doc,
	                  xercesc::DOMElement* root,
	                  ProblemSetup& sim_data);

  private:
	/// Output file
	std::string _output_file;
}; // class InputOutput

}
} // namespaces

#endif // STATE_H_INCLUDED
