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
 * @brief Particles files manager.
 * (See Aqua::InputOutput::Particles for details)
 */

#ifndef PARTICLES_H_INCLUDED
#define PARTICLES_H_INCLUDED

#include <vector>
#include <map>
#include "aquagpusph/sphPrerequisites.hpp"
#include "aquagpusph/ProblemSetup.hpp"
#include "Logger.hpp"
#include "InputOutput.hpp"

namespace Aqua {
namespace InputOutput {

/** \class Particles Particles.h InputOutput/Particles.h
 * @brief Particles file loader/saver base class.
 *
 * In AQUAgpusph the input/output managers are divided in 3 different types:
 *   -# The simulation configuration files manager
 *   -# The report file managers
 *   -# The particles output file managers
 *
 * The particles files have 2 main objectives:
 *   -# Particles data loading at the start of simulations.
 *   -# Visualization of the simulation results.
 *
 * @see Aqua::InputOutput::InputOutput
 * @see Aqua::InputOutput::Report
 * @see Aqua::InputOutput::State
 */
class Particles : public InputOutput
{
  public:
	/** @brief Constructor
	 * @param sim_data Simulation data
	 * @param offset First particle managed by this saver/loader
	 * @param iset Particles set index
	 * @param n Number of particles managed by this saver/loader. If 0,
	 * the number of particles will be obtained from the input file (thus only
	 * valid for loaders)
	 */
	Particles(ProblemSetup& sim_data,
	          unsigned int iset,
	          unsigned int offset,
	          unsigned int n = 0);

	/// Destructor
	virtual ~Particles();

	/** @brief Get the last printed file path.
	 * @return The last printed file, NULL if a file has not been printed yet.
	 */
	const std::string file() { return _output_file; }

	/** @brief Get the number of particles managed by this instance
	 * @return Number of particles
	 */
	unsigned int n() { return _bounds.y - _bounds.x; }

	/** @brief Get the user event to be waited for before the file saving is
	 * finished
	 * @return The event
	 * @note In general ::waitForSavers() shall be used instead of this
	 * function, which is provided just to work with the OpenCL callbacks
	 */
	inline cl_event getUserEvent() const { return _user_event; }

	/** @brief Wait for the eventual parallel saving threads.
	 *
	 * Some savers may optionally launch parallel threads to save the data, in
	 * an asynchronous way, in order to improve the performance. In such a case,
	 * AQUAgpusph shall wait them to finish before proceeding to destroy the
	 * data
	 */
	virtual inline void waitForSavers()
	{
		if (!_user_event)
			return;
		cl_int err_code;
		err_code = clWaitForEvents(1, &_user_event);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure waiting for a file writer\n");
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		err_code = clReleaseEvent(_user_event);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure releasing the user event\n");
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		_user_event = NULL;
	}

	/** @brief Save the data.
	 *
	 * @param t Simulation time
	 */
	virtual void save(float t);

	/** @brief Print the data to a file
	 *
	 * This function shall be overloaded by the inherited classes. Remember to
	 * call this function anyway to set the user event as completed
	 * @note This method is public to work with the OpenCL callbacks, but it is
	 * not meant to be called by the users
	 */
	virtual void print_file()
	{
		cl_int err_code;
		err_code = clSetUserEventStatus(_user_event, CL_COMPLETE);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure setting the user event as completed.\n");
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		err_code = clReleaseEvent(_user_event);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure releasing the user event.\n");
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}

  protected:
	/** @brief Get the simulation data structure
	 *
	 * @return Simulation data
	 */
	ProblemSetup& simData() { return _sim_data; }

	/** @brief Set the number of particles managed by this instance
	 * @return Number of particles
	 */
	inline void n(unsigned int n)
	{
		_bounds.y = _bounds.x + n;
		for (auto const& mem : _data)
			free(mem.second);
		_data.clear();
	}

	/** @brief Get the particle index bounds of the "set of particles" managed
	 * by this class.
	 * @return The index bounds (first and last particle).
	 */
	inline const uivec2 bounds() const { return _bounds; }

	/** @brief Get the "particles set" index associated with this class
	 * @return The "particles index" index.
	 */
	inline const unsigned int setId() const { return _iset; }

	/** @brief Register some default arrays:
	 *   -# iset
	 *   -# id_sorted
	 *   -# id_unsorted
	 */
	void loadDefault();

	/** @brief Set the file name.
	 * @param filename The new file to save/load. Optionally a null parameter
	 * can be passed in order to clear the stored file name.
	 */
	inline void file(const std::string filename) { _output_file = filename; };

	/** Look for the first non-existing file path
	 * @param basename The base name of the file
	 * @param start_index First index that will be checked.
	 * @param digits Number of digits of the replaced integer number. If the
	 * number of digits of the integer value are greater than this value this
	 * parameter will be ignored, otherwise zeroes will be appended at the left
	 * of the decimal representation of the integer.
	 * @return The next non-existing file index.
	 * @see Aqua::newFilePath()
	 */
	unsigned int file(const std::string basename,
	                  unsigned int start_index,
	                  unsigned int digits = 5);

	/** @brief Get the current simulation time to be written
	 * @return The simulation time
	 */
	inline const float time() const { return _time; }

	/** @brief Download the data from the device and store it
	 * @param fields Fields to download
	 * @return The event that will be marked as completed when the data is
	 * already available
	 * @see ::data()
	 */
	cl_event download(std::vector<std::string> fields);

	/** @brief Get the stored memory objects where the device data has been
	 * downloaded
	 * @return The memory objects
	 * @see ::download()
	 */
	inline std::map<std::string, void*> data() const { return _data; }

  private:
	/// Simulation data
	ProblemSetup _sim_data;

	/// Particles managed bounds
	uivec2 _bounds;

	/// Fluid index
	unsigned int _iset;

	/// Last file printed
	std::string _output_file;

	/// Current simulation time to be written
	float _time;

	/** The user event to be marked as finished when the last call to
	 * InputOutput::save() is dispatched
	 */
	cl_event _user_event;

	/// List of host allocated memories
	std::map<std::string, void*> _data;

}; // class InputOutput

}
} // namespaces

#endif // PARTICLES_H_INCLUDED
