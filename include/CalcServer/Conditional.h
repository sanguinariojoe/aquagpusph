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
 * @brief Check a condition to enable/disable all the tools in its scope
 * (See Aqua::CalcServer::Conditional and Aqua::CalcServer::End for details)
 */

#ifndef IF_H_INCLUDED
#define IF_H_INCLUDED

#include <CalcServer/Tool.h>

namespace Aqua {
namespace CalcServer {

/** @class Conditional Conditional.h CalcServer/Conditional.h
 * @brief Base class for conditional tools like While or If
 *
 * If the result of evaluating the condition expression is equal to 0, the
 * result will be considered false, and therefore all the tools until the next
 * End tool will be disabled
 */
class Conditional : public Aqua::CalcServer::Tool
{
  public:
	/** @brief Constructor.
	 * @param name Tool name.
	 * @param condition Condition to evaluate. If the result is 0, false will be
	 * considered and all the subsequent tools will be disabled until an End
	 * tool is reached.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Conditional(const std::string name,
	            const std::string condition,
	            bool once = false);

	/// Destructor.
	~Conditional();

	/** @brief Initialize the tool.
	 */
	void setup();

	/** Get the next tool to be executed in the pipeline.
	 *
	 * Depending on the condition value, that tool is the next one or the tool
	 * after the End tool which closes the scope
	 * @return Next tool to be executed. NULL if no more tools shall be executed
	 * in the pipeline
	 */
	virtual Tool* next_tool();

	/** Open a new scope
	 * @return 1
	 * @note The scope shall be close at some point by an EndIf tool
	 */
	const int scope_modifier() { return 1; }

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/// Condition expression to evaluate
	std::string _condition;
	/// The next tool in the pipeline when the condition is not fulfilled
	Tool* _ending_tool;

  protected:
	/// Condition result
	bool _result;
};

/** @class While Conditional.h CalcServer/Conditional.h
 * @brief Execute all the tools in its scope until the condition becomes
 * unfulfilled
 * @see Aqua::CalcServer::Conditional
 */
class While : public Aqua::CalcServer::Conditional
{
  public:
	/** @brief Constructor.
	 * @param name Tool name.
	 * @param condition Condition to evaluate. If the result is 0, false will be
	 * considered and all the subsequent tools will be disabled until an End
	 * tool is reached.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	While(const std::string name,
	      const std::string condition,
	      bool once = false);

	/// Destructor.
	~While();

	/** @brief Initialize the tool.
	 */
	void setup();
};

/** @class If Conditional.h CalcServer/Conditional.h
 * @brief Execute all the tools in its scope if the condition is fulfilled
 * @see Aqua::CalcServer::Conditional
 */
class If : public Aqua::CalcServer::Conditional
{
  public:
	/** @brief Constructor.
	 * @param name Tool name.
	 * @param condition Condition to evaluate. If the result is 0, false will be
	 * considered and all the subsequent tools will be disabled until an End
	 * tool is reached.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	If(const std::string name, const std::string condition, bool once = false);

	/// Destructor.
	~If();

	/** @brief Initialize the tool.
	 */
	void setup();

	/** Get the next tool to be executed in the pipeline.
	 *
	 * Such tool will be the next one if the condition is fulfilled, or the tool
	 * after the Aqua::CalcServer::End tool which closes the scope otherwise.
	 *
	 * Since the closing Aqua::CalcServer::End tool is invariably giving back
	 * the control to this tool, the condition will be marked as unfulfilled
	 * when this tool is called, avoiding the loop execution
	 * @return Next tool to be executed. NULL if no more tools shall be executed
	 * in the pipeline
	 */
	Tool* next_tool();

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);
};

/** @class End Conditional.h CalcServer/Conditional.h
 * @brief Close the scope open by a previous conditional tool, like While or If
 *
 * This tool is always giving back the control to the conditional tool, which is
 * responsible to redirect the pipeline
 */
class End : public Aqua::CalcServer::Tool
{
  public:
	/** @brief Constructor.
	 * @param name Tool name.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	End(const std::string name, bool once = false);

	/// Destructor.
	~End();

	/** @brief Initialize the tool.
	 */
	void setup();

	/** Close an already open scope
	 * @return -1
	 */
	const int scope_modifier() { return -1; }
};

}
} // namespace

#endif // IF_H_INCLUDED
