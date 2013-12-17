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

#ifndef XML_H_INCLUDED
#define XML_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------------
// Include the h5part library
// ----------------------------------------------------------------------------
#ifdef HAVE_H5PART
	#include <H5Part.h>
	#include <hdf5.h>
#endif // HAVE_H5PART

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
// Include the Host-Server fluid transfer layer
// ----------------------------------------------------------------------------
#include <Fluid.h>

// ----------------------------------------------------------------------------
// Include the equation reader
// ----------------------------------------------------------------------------
#include <Tokenizer/Tokenizer.h>

namespace Aqua{ namespace InputOutput{ namespace Input{

/** Reads XML file to set the fluid.
 * @param path Input file path.
 * @param ifluid Index of the fluid.
 * @param i0 Starting index of the particle.
 * @param n number of particles to read.
 * @param refd Reference density.
 * @param h Reference kernel height.
 * @param F Fluid host instance.
 * @return false if all gone right. \n true otherwise.
 */
int loadXML(const char* path, int ifluid, unsigned int i0, unsigned int n, float refd, float h, Fluid *F);

}}} // namespace

#endif // XML_H_INCLUDED
