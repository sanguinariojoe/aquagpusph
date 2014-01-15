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

#ifndef TOKENIZER_H_INCLUDED
#define TOKENIZER_H_INCLUDED

#include <map>
#include <string>

/** \class Tokenizer Tokenizer.h Tokenizer/Tokenizer.h
 * @brief Tool to evaluate math expressions.
 */
class Tokenizer
{
public:
	/** Constructor
	 */
	Tokenizer();

	/** Destructor
	 */
	~Tokenizer();

	/** Register a variable to be used later, or modify its value if it has
	 * been already generated.
	 * @param name Name of the variable.
	 * @param value Value of the variable.
	 * @return true if the variable already exists, false otherwise.
	 */
	bool registerVariable(const char* name, float value);

	/** Unregister a variable.
	 * @param name Name of the variable.
	 * @return true if the variable has been unregistered, false if the
	 * variable cannot be unregistered (for instance because it does not
     * exist)
	 */
	bool unregisterVariable(const char* name);

    /** Clear all the registered variables.
     */
    void clearVariables();

	/** Returns if a variable has been registered.
	 * @param name Name of the desired variable
	 * @return true if exist a variable with the given name, false otherwise.
	 */
	bool isVariable(const char* name);

	/** Returns a variable value.
	 * @param name Name of the desired variable
	 * @return The Value of the variable, or 0.f if the variable cannot be found.
	 */
	float variable(const char* name);

	/** Solve the math expresion.
	 * @param eq Math expression to solve.
	 * @return Expression value, 0.f if failed (it will be reported).
	 */
	float solve(const char* eq);

protected:
    /** Register the default variables
     */
    virtual void defaultVariables();

private:
	/// Registered variables
	std::map<std::string, float> _variables;
};   // class Tokenizer

#endif // TOKENIZER_H_INCLUDED
