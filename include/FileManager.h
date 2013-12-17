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

#ifndef FILEMANAGER_H_INCLUDED
#define FILEMANAGER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>
#include <unistd.h>
#include <errno.h>
#include <deque>

#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// h5part library
// ----------------------------------------------------------------------------
#ifdef HAVE_H5PART
	#include <H5Part.h>
	#include <hdf5.h>
#endif

// ----------------------------------------------------------------------------
// xerces XML parser library
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

#include <Singleton.h>
#include <ProblemSetup.h>
#include <Output.h>
#include <AuxiliarMethods.h>

namespace Aqua{
/// @namespace InputOutput Input/Output interfaces.
namespace InputOutput{

/** \class FileManager FileManager.h FileManager.h
 *  Input/Output files manager.
 */
class FileManager : public Aqua::Singleton<Aqua::InputOutput::FileManager>
{
public:
	/** Constructor
	 */
	FileManager();

	/** Destructor
	 */
	~FileManager();

	/** Parse the XML problem definition input file.
	 * The XML file may be structured in several sub-files referenced by
	 * include tags.
	 * @return false if the XML definition has been succesfully loaded, true
	 * otherwise.
	 */
	bool parse();

	/** Creates and initializates the log file.
	 * @return false if the file has been succesfully created, true otherwise.
	 */
	bool startLog();

	/** Open the output files. The output files can be the H5Part ones, the
	 * Tecplot ones, or the VTK ones.
	 * @return false if the files have been succesfully created, true otherwise.
	 */
	bool openFiles();

	/** Close the output files. This function is critical for some formats which
	 * requires an ending data write to be readable.
	 * @return false if the files have been succesfully closed, true otherwise.
	 */
	bool closeFiles();

	/** Set the input file.
	 * @param path Path to input file.
	 */
	void inputFile(const char* path){strcpy(_in_file, path);}

	/** Get input file.
	 * @return Path to input file.
	 */
	const char* inputFile(){return (const char*)_in_file;}

	/** Set the output prefix. The prefix can be used to group all the output
	 * files using the alphanumeric sorting methods of the operating system.
	 * @param prefix Output files prefix.
	 * @remarks This method must be called before the openFiles() one. Othewise
	 * the files will not be renamed.
	 */
	void outputPrefix(const char* prefix){strcpy(_out_prefix, prefix);}

	/** Get the output prefix.
	 * @return Output files prefix.
	 */
	const char* outputPrefix(){return _out_prefix;}

	/** Set the number of printed files.
	 * @param n Number of files printed
	 */
	void nOutput(unsigned int n){_n_out=n;}

	/** Get the number of printed files.
	 * @return Number of output files printed.
	 */
	unsigned int nOutput(){return _n_out;}

	/** Get the log file name.
	 * @return Log file name.
	 */
	const char* logFileName(){return (const char*)_log_file;}

	/** Get the energy file name.
	 * @return Energy file name.
	 */
	const char* enFileName(){return (const char*)_en_file;}

	/** Get the log file decriptor.
	 * @return Log file decriptor.
	 */
	FILE* logFile(){return _log_file_id;}

	/** Get the energy file decriptor.
	 * @return Energy file decriptor.
	 */
	FILE* enFile(){return _en_file_id;}

	/** Get the bounds file decriptor.
	 * @return Bounds file decriptor.
	 */
	FILE* boundsFile(){return _bounds_file_id;}

	#ifdef HAVE_H5PART
	    /** Get the H5Part file decriptor. Of course this file will
	     * only exist if the H5Part output has been activated.
	     * @return H5Part file decriptor.
	     */
	    H5PartFile* h5File(){return _H5Part_file_id;}
	#endif

	/** Switch on/off the H5Part output file.
	 * @param flag true if the H5Part output must be performed, false otherwise.
	 * @remarks This method does not take into account if the H5Part support has
	 * been activated during the program compilation.
	 */
	void isH5Part(bool flag){_must_print_H5=flag;}

	/** Get if H5Part output is enabled/disabled.
	 * @return true if H5Part output must be performed, false otherwise.
	 */
	bool isH5Part(){return _must_print_H5;}

	/** Switch on/off the VTK output.
	 * @param flag true if the VTK output must be performed, false otherwise.
	 * @remarks This method does not take into account if VTK support has
	 * been activated during the program compilation.
	 */
	void isVTK(bool flag){_must_print_VTK=flag;}

	/** Get if the VTK output is enabled/disabled.
	 * @return true if the VTK output must be performed, false otherwise.
	 */
	bool isVTK(){return _must_print_VTK;}

	/** Switch on/off the Tecplot output.
	 * @param flag true if the Tecplot output must be performed, false otherwise.
	 * @warning Tecplot output is outdated and it will be removed in future
	 * releases.
	 */
	void isTecplot(bool flag){_must_print_Tecplot=flag;}

	/** Get if the Tecplot output is enabled/disabled.
	 * @return true if the Tecplot output must be performed, false otherwise.
	 * @warning Tecplot output is outdated and it will be removed in future
	 * releases.
	 */
	bool isTecplot(){return _must_print_Tecplot;}

	/** Mark a field to be excluded from the output files.
	 * @param field Excluded field name.
	 */
	void excludeField(const char* field);

	/** Get if a field must be printed in the output files.
	 * @param field Field name.
	 * @return true if the field must be printed, false otherwise.
	 */
	bool isField(const char* field);

protected:

	/** Starts the log file web page.
	 */
	void initLogFile();

	/** Finish the log file web page.
	 */
	void endLogFile();

	/** Starts the Energy report file.
	 */
	void initEnFile();

	/** Finish the Energy report file.
	 */
	void endEnFile();

	/** Starts the Bounds report file.
	 */
	void initBoundsFile();

	/** Finish the Bounds report file.
	 */
	void endBoundsFile();

private:
	/** Look for general settings sections.
	 * @param root root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseSettings(xercesc::DOMElement *root);
	/** Look for OpenCL settings sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseOpenCL(xercesc::DOMElement *root);
	/** Look for time control sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseTiming(xercesc::DOMElement *root);
	/** Look for SPH settings sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseSPH(xercesc::DOMElement *root);
	/** Look for fluids sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseFluid(xercesc::DOMElement *root);
	/** Look for sensors sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseSensors(xercesc::DOMElement *root);
	/** Look for movements sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseMotions(xercesc::DOMElement *root);
	/** Look for portals sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parsePortals(xercesc::DOMElement *root);
	/** Look for ghost particles sections.
	 * @param root Root XML node.
	 * @return false if all gone right, true otherwise
	 */
	bool parseGhostParticles(xercesc::DOMElement *root);

	/// Name of the input file
	char* _in_file;
	/// Prefix of the output file
	char* _out_prefix;
	/// Number of output files printed
	unsigned int _n_out;
	/// Name of the output log file
	char* _log_file;
	/// Name of the output energy file
	char* _en_file;
	/// Name of the output bounds file
	char* _bounds_file;
	/// Log file handle
	FILE *_log_file_id;
	/// Energy file handle
	FILE *_en_file_id;
	/// Bounds file handle
	FILE *_bounds_file_id;
	#ifdef HAVE_H5PART
        /// H5Part output file handle
	    H5PartFile *_H5Part_file_id;
	#endif
	/// H5Parts files printing flag
	bool _must_print_H5;
	/// Tecplot files printing flag
	bool _must_print_VTK;
	/// Tecplot files printing flag
	bool _must_print_Tecplot;
	/// List of fields excluded from the output files
	std::deque<char*> _excluded_fields;
};  // class FileManager

}}  // namespaces

#endif // FILEMANAGER_H_INCLUDED
