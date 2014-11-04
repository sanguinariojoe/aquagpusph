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

#include <InputOutput/InputOutput.h>

namespace Aqua{
namespace InputOutput{

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
 * files in the examples provided with the package, in the user manual (chapter
 * 4), and in the Aqua::InputOutput::ProblemSetup class documentation.
 *
 * @see Aqua::InputOutput::InputOutput
 * @see Aqua::InputOutput::Report
 * @see Aqua::InputOutput::Particles
 */
class State : public InputOutput
{
public:
    /// Constructor
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
     * @see Aqua::InputOutput::Particles::save()
     * @return false if all gone right, true otherwise.
     */
    bool save();

    /** @brief Load the simulation configuration data.
     * @return false if all gone right, true otherwise.
     */
    bool load();

protected:
    /** @brief Parse the XML file
     * @param filepath file to be parsed.
     * @return false if all gone right, true otherwise
     * @note This function is calling himself for each `<Include>` tag found,
     * conveniently changing the input file name @paramname{filepath}.
     */
    bool parse(const char* filepath);

    /** @brief Parse the general settings sections.
     * @param root root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphSettings
     */
    bool parseSettings(xercesc::DOMElement *root);

    /** @brief Parse the variables sections.
     * @param root root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphVariables
     */
    bool parseVariables(xercesc::DOMElement *root);

    /** @brief Parse the tools sections.
     * @param root root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphTool
     */
    bool parseTools(xercesc::DOMElement *root);

    /** @brief Parse the time control sections.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
     */
    bool parseTiming(xercesc::DOMElement *root);

    /** Look for particles set sections.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphParticlesSet
     */
    bool parseSet(xercesc::DOMElement *root);

    /** @brief Write the XML file
     * @param filepath file to be written.
     * @return false if all gone right, true otherwise
     */
    bool write(const char* filepath);

    /** @brief Write the settings section.
     * @param doc XML generated document.
     * @param root root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphSettings
     */
    bool writeSettings(xercesc::DOMDocument* doc, xercesc::DOMElement *root);

    /** @brief Write the variables section.
     * @param doc XML generated document.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphVariables
     */
    bool writeVariables(xercesc::DOMDocument* doc, xercesc::DOMElement *root);

    /** @brief Write the tools section.
     * @param doc XML generated document.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphTool
     */
    bool writeTools(xercesc::DOMDocument* doc, xercesc::DOMElement *root);

    /** @brief Write the time control section.
     * @param doc XML generated document.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphTimingParameters
     */
    bool writeTiming(xercesc::DOMDocument* doc, xercesc::DOMElement *root);

    /** @brief Write the particles set sections.
     * @param doc XML generated document.
     * @param root Root XML node.
     * @return false if all gone right, true otherwise
     * @see Aqua::InputOutput::ProblemSetup::sphParticlesSet
     */
    bool writeSet(xercesc::DOMDocument* doc, xercesc::DOMElement *root);

private:
    /// Output file
    char* _output_file;
};  // class InputOutput

}}  // namespaces

#endif // STATE_H_INCLUDED
