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
 * @brief Particles plain text data files loader/saver.
 * (See Aqua::InputOutput::FastASCII for details)
 */

#pragma once

#include "aquagpusph/sphPrerequisites.hpp"
#include "FastASCII.hpp"

namespace Aqua {
namespace InputOutput {

/** @class CSV CSV.h InputOutput/CSV.h
 * @brief CSV particles data files loader/saver.
 *
 * These files are formatted as ASCCI plain text where the particles data are
 * stored by rows, and where the fields are separated by columns.
 * @remarks In the case of integer numbers (signed or unsigned) this class does
 * not care about decimal points, just truncating the value, i.e. 1.5 will be
 * interpreted as 1, and -1.5 will be interpreted as -1.
 * @warning Saving the particles data in plain text format may be heavily hard
 * disk demanding, and therefore it is strongly recommended to consider binary
 * formats like Aqua::InputOutput::VTK.
 */
class CSV : public FastASCII
{
  public:
	/** @brief Constructor
	 * @param sim_data Simulation data
	 * @param iset Particles set index.
	 * @param offset First particle managed by this saver/loader.
	 * @param n Number of particles managed by this saver/loader. If 0,
	 * the number of particles will be obtained from the input file (thus only
	 * valid for loaders)
	 */
	CSV(ProblemSetup& sim_data,
	    unsigned int iset,
	    size_t offset,
	    size_t n = 0,
	    const std::string file_ext = ".csv",
	    const char sep = ',');

	/// Destructor
	~CSV();

	/** @brief Load the data.
	 */
	void load() final;

	/** @brief Print the data to a file
	 * @param sep Fields separator
	 * @param comp_sep Components separator (for vectorial types)
	 * @note This method is public to work with the OpenCL callbacks, but it is
	 * not meant to be called by the users
	 */
	void print_file() final;

  protected:
	/** @brief Write the file header
	 * @param f The file handler
	 */
	void print_header(std::ofstream& f) const final;

	/** @brief Conveniently format a read line.
	 * @param l Line text.
	 */
	void formatLine(std::string& l) final;

	/** @brief Count the number of particles present in the input file.
	 * @param f File to be read.
	 * @return The number of particles found in the file.
	 */
	size_t readNParticles(std::ifstream& f) final;

  private:
	/// The separator
	char _sep;

	/// A flag to mark whether the file header shall be still looked for
	bool _has_header;

}; // class InputOutput

} // InputOutput::
} // Aqua::
