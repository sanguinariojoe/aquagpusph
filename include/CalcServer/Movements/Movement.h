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
 * @brief Motions base class.
 * (See Aqua::CalcServer::Movement::Movement for details)
 */

#ifndef MOVEMENT_H_INCLUDED
#define MOVEMENT_H_INCLUDED

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

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{
/** @namespace Movement
 * @brief Motions.
 */
namespace Movement {

/** @class Movement Movement.h CalcServer/Movements/Movement.h
 * @brief Base class for all the motion models.
 *
 * This class is parsing the input XML file to find the used OpenCL script, but
 * the inherited motion classes may parse more fields with the method _parse().
 */
class Movement : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Movement();

    /// Destructor.
    virtual ~Movement();

    /** @brief Parse the XML definition file.
     *
     * This method will parse the XML file, looking for errors, and extracting
     * the OpenCL kernel to be used.
     *
     * After that the control is passed to _parse() function.
     *
     * @param def XML file with the motion definition.
     * @return false if all gone right, true otherwise.
     * @remarks This method is not virtual, consider overloading _parse()
     * function.
     * @see _parse()
     */
    bool parse(const char* def);

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    virtual bool execute()=0;

protected:
    /** @brief Parse the additional data in the input definition file.
     *
     * Since parse() function is just looking for the OpenCL script to be used,
     * this method should be overloaded to get additional data.
     *
     * @param root Input node of the parser.
     * @return false if all gone right, true otherwise.
     * @see parse()
     */
    virtual bool _parse(xercesc::DOMElement *root)=0;

    /// OpenCL program
    cl_program _program;
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;

private:
    /** @brief Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    virtual bool setupOpenCL();

    /// OpenCL script path
    char* _path;

};  // class Movement

}}} // namespaces

#endif // MOVEMENT_H_INCLUDED
