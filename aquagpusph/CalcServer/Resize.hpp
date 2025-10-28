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
 * @brief Resize an array.
 * (See Aqua::CalcServer::Resize for details)
 */

#ifndef RESIZE_H_INCLUDED
#define RESIZE_H_INCLUDED

#include "CalcServer.hpp"
#include "SetScalar.hpp"
#include "Kernel.hpp"

namespace Aqua {
namespace CalcServer {

/** @class Resize Resize.h CalcServer/Resize.h
 * @brief Resize an array to a desired length.
 *
 * Incremental resizing can be optionally asked (enabled by default), in such
 * a way the array is only resized if the new size is larger than the already
 * available one
 */
class Resize final : public Aqua::CalcServer::ScalarExpression
{
  public:
	/** Constructor.
	 * @param name Tool name.
	 * @param var_name Variable to set.
	 * @param length The length expression. It is always evaluated as a
	 * "size_t"
	 * @param shrink Whether the array shall be shrinked or not
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Resize(const std::string name,
	    const std::string var_name,
	    const std::string value,
	    bool shrink = false,
	    bool once = false);

	/** Destructor.
	 */
	~Resize();

	/** Initialize the tool.
	 */
	void setup();

  protected:
	/** @brief Evaluate the expression if possible
	 */
	void _solve();

	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/** Get the input variable
	 */
	void variable();

	/// Input variable name
	std::string _var_name;
	/// Length expression
	std::string _value;
	/// Shall we shrink?
	bool _shrink;

	/// Input variable
	InputOutput::ArrayVariable* _var;
};

}
} // namespace

#endif // RESIZE_H_INCLUDED
