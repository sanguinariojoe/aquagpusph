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
 * @brief Math expression evaluator.
 * (See Aqua::Tokenizer for details)
 */

#ifndef TOKENIZER_H_INCLUDED
#define TOKENIZER_H_INCLUDED

#include <map>
#include <string>
#include <muParser.h>
#include <mutex>
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {

/** \class Tokenizer Tokenizer_muparser.hpp Tokenizer/Tokenizer_muparser.hpp
 * @brief Tool to evaluate math expressions.
 *
 * This tool is based in libmatheval, to learn more please visit the following
 * web page:
 *
 * http://www.gnu.org/software/libmatheval
 */
class Tokenizer_muparser
{
  public:
	/// Constructor
	Tokenizer_muparser();

	/// Destructor
	~Tokenizer_muparser();

	/** @brief Register a variable.
	 *
	 * In case that the variable already exist, it will be modified.
	 *
	 * The registered variables can be used later in the expression to be
	 * evaluated.
	 *
	 * @param name Name of the variable.
	 * @param value Value of the variable.
	 * @return true if the variable already exists, false otherwise.
	 */
	template<typename T>
	bool registerVariable(const std::string name, T value)
	{
		const std::lock_guard<std::mutex> lock(this->mutex);
		bool overwritten = false;
		if (isVariable(name))
			overwritten = true;
		p.DefineConst(name, (mu::value_type)value);
		return overwritten;
	}

	/** @brief Get the list of variables used on an mathematical expression.
	 * @param eq Math expression to solve.
	 * @return List of variable names
	 */
	std::vector<std::string> exprVariables(const std::string eq);

	/** @brief Solve a math expression.
	 * @param eq Math expression to solve.
	 * @return Expression value, 0.0 if the evaluation failed (it will be
	 * reported by terminal).
	 */
	template<typename T=float>
	T solve(const std::string eq)
	{
		const std::lock_guard<std::mutex> lock(this->mutex);
		T res;

		// First try straight number conversions
		size_t sz;
		try {
			res = narrow_cast<T>((int64_t)std::stoll(eq, &sz));
			if (sz == eq.size()) {
				return res;
			}
		} catch (...) {
		}
		try {
			std::string::size_type sz;
			res = narrow_cast<T>(std::stod(eq, &sz));
			if (sz == eq.size()) {
				return res;
			}
		} catch (...) {
		}

		// No way, let's evaluate it as an expression
		mu::value_type res_org;
		p.SetExpr(eq);
		try {
			res_org = p.Eval();
			res = narrow_cast<T>(res_org);
		} catch (mu::Parser::exception_type& e) {
			std::ostringstream msg;
			msg << "Error evaluating \"" << e.GetExpr() << "\"" << std::endl;
			msg << "input = \"" << eq << "\"" << std::endl;
			LOG(L_WARNING, msg.str());
			msg.str("");
			msg << "\t" << e.GetMsg() << std::endl;
			LOG0(L_DEBUG, msg.str());
			msg.str("");
			msg << "\tToken " << e.GetToken() << " in position " << e.GetPos()
				<< std::endl;
			LOG0(L_DEBUG, msg.str());
			throw;
		} catch(std::out_of_range&) {
			LOG(L_ERROR, std::string("Error evaluating \"") + eq + "\":\n");
			LOG0(L_DEBUG, "\n");
			LOG0(L_DEBUG, std::string("The result ") +
			     std::to_string(res_org) + " overflows the type \"" +
			     typeid(res).name() + "\"\n");
			LOG0(L_DEBUG, "\n");
			throw;
		}

		return res;
	}

  protected:
	/** @brief Checks if a variable has been registered.
	 * @param name Name of the variable
	 * @return true if already exist a variable with the given name, false
	 * otherwise.
	 */
	bool isVariable(const std::string name);

	/** @brief Returns a variable value.
	 * @param name Name of the variable
	 * @return The Value of the variable, or 0.0 if the variable cannot be
	 * found.
	 */
	template<typename T=float>
	T variable(const std::string name)
	{
		if (!isVariable(name))
			return 0.f;
		mu::valmap_type cmap = p.GetConst();
		return narrow_cast<T>(cmap[name.c_str()]);
	}

  private:
	/// Mathematical expressions variables inspector
	mu::Parser q;
	/// Mathematical expressions parser
	mu::Parser p;

	/// Shared mutex to avoid race conditions
	static std::mutex mutex;
}; // class Tokenizer

} // Aqua::

#endif // TOKENIZER_H_INCLUDED
