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
#include <mutex>
#include "aquagpusph/ext/exprtk.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {

typedef long double tokenizer_t;

/** \class Tokenizer Tokenizer.h Tokenizer/Tokenizer.h
 * @brief Tool to evaluate math expressions.
 *
 * This tool is based in libmatheval, to learn more please visit the following
 * web page:
 *
 * http://www.gnu.org/software/libmatheval
 */
class Tokenizer_exprtk
{
  public:
	/// Constructor
	Tokenizer_exprtk();

	/// Destructor
	~Tokenizer_exprtk();

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
	void registerVariable(const std::string name, T value)
	{
		tokenizer_t val = narrow_cast<tokenizer_t>(value);
		const std::lock_guard<std::mutex> lock(this->mutex);
		removeVariable(name);
		vars.add_variable(name, val);
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
		exprtk::parser<tokenizer_t> parser;
		exprtk::expression<tokenizer_t> expr;

		expr.register_symbol_table(vars);
		// expr.register_symbol_table(ivars);

		if (!parser.compile(eq, expr))
		{
			LOG(L_ERROR, std::string("Error parsing \"") + eq + "\":\n");
			LOG0(L_DEBUG, "\n");
			LOG0(L_DEBUG, parser.error() + "\n");
			LOG0(L_DEBUG, "\n");
			throw std::runtime_error("Invalid expression");
		}

		T res = 0;
		try {
			res = narrow_cast<T>(expr.value());
			return narrow_cast<T>(res);
		} catch(std::out_of_range) {
			LOG(L_ERROR, std::string("Error parsing \"") + eq + "\":\n");
			LOG0(L_DEBUG, "\n");
			LOG0(L_DEBUG, std::string("The result ") + std::to_string(res) +
				"overflows the type float\n");
			LOG0(L_DEBUG, "\n");
			throw;
		}
		return res;
	}

  protected:
	/** @brief Remove a variable if it has been registered.
	 * @param name Name of the variable
	 */
	inline void removeVariable(const std::string name)
	{
		if (vars.is_variable(name))
			vars.remove_variable(name);
	}

  private:
	/// List of registered variables
	/// @{
	exprtk::symbol_table<tokenizer_t> vars;
	/// @}

	/// Shared mutex to avoid race conditions
	static std::mutex mutex;
}; // class Tokenizer

} // Aqua::

#endif // TOKENIZER_H_INCLUDED
