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
 * @brief Base class for all the report file managers.
 * (See Aqua::InputOutput::Report for details)
 */

#ifndef REPORT_H_INCLUDED
#define REPORT_H_INCLUDED

#include <string>
#include <stdexcept>
#include "aquagpusph/sphPrerequisites.hpp"
#include "InputOutput.hpp"

namespace Aqua {
namespace InputOutput {

/** \class Report Report.h InputOutput/Report.h
 * @brief Base class for all the report file managers.
 *
 * In AQUAgpusph the input/output managers are divided in 3 different types:
 *   -# The simulation configuration files manager
 *   -# The report file managers
 *   -# The particles output file managers
 *
 * The reports are just output files that are relatively low hard disk
 * demanding.
 *
 * This class allows the inherited ones to force the output file name (see
 * file(const char* filename)) or to look the first file name which is not
 * already existing in the system (see
 * file(const char* basename, unsigned int startindex)).
 *
 * @see Aqua::InputOutput::InputOutput
 * @see Aqua::InputOutput::State
 * @see Aqua::InputOutput::Particles
 */
class Report : public InputOutput
{
  public:
	/** @brief Save the data.
	 *
	 * @param t Simulation time
	 */
	virtual void save(float t) = 0;

	/** @brief Load the data.
	 *
	 * Since the reports are mainly output files, the load method should be
	 * useless, and therefore this class provide a way to can omit it in the
	 * inherited ones.
	 */
	virtual void load() { throw std::logic_error("Report cannot load files"); }

	/** @brief Get the used output file path.
	 * @return The report file, NULL if it is not a file.
	 */
	std::string file() { return _output_file; }

  protected:
	/// Constructor
	Report();

	/// Destructor
	virtual ~Report();

	/** @brief Set the report file name.
	 * @param filename The file name. Optionally @paramname{filename} = null can
	 * be set in order to clear the stored file name.
	 */
	void file(std::string filename);

	/** @brief Look for the first non-existing file name
	 *
	 * Several scape strings can be applied to look for the appropriate name:
	 *   - `"%d"`/`"{index}"` will be replaced by the first integer which
	 *     results in a non-existing file path
	 *   - `"{mpi_rank}"` will be replaced by the process identifier index
	 *
	 * @param basename The base name of the file
	 * @param start_index First index that will be checked
	 */
	void file(std::string basename, unsigned int start_index);

  private:
	/// Last file printed
	std::string _output_file;

}; // class InputOutput

}
} // namespaces

#endif // REPORT_H_INCLUDED
