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

#ifndef MOVEMENT_H_INCLUDED
#define MOVEMENT_H_INCLUDED

// ----------------------------------------------------------------------------
// Include xerces XML parser library
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{
/// @namespace Move Movements space.
namespace Movement {

/** @class Movement Movement.h CalcServer/Movements/Movement.h
 * @brief Abstraction class for all the motion models. As base class, only
 * store some usefull variables as the OpenCL kernel to use.
 */
class Movement : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Movement();

	/** Destructor.
	 */
	virtual ~Movement();

	/** Parse definition file
	 * @param def XML file with the motion definition.
	 * @return false if all gone right, true otherwise.
	 */
	bool parse(const char* def);

	/** Execute the motion.
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool execute()=0;

protected:
	/** Parse the input definition file (motion type specific data).
	 * @param root Input node of the parser.
	 * @return false if all gone right, true otherwise.
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
	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool setupOpenCL();

	/// OpenCL script path
	char* _path;

};  // class Movement

}}} // namespaces

#endif // MOVEMENT_H_INCLUDED
