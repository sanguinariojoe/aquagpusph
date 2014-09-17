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

namespace Aqua{

/** \class Tokenizer Tokenizer.h Tokenizer/Tokenizer.h
 * @brief Tool to evaluate math expressions.
 *
 * This tool is based in libmatheval, to learn more please visit the following
 * web page:
 *
 * http://www.gnu.org/software/libmatheval
 */
class Tokenizer
{
public:
    /// Constructor
    Tokenizer();

    /// Destructor
    ~Tokenizer();

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
    bool registerVariable(const char* name, float value);

    /** @brief Unregister a variable.
     * @param name Name of the variable.
     * @return true if the variable has been unregistered, false if the
     * variable cannot be unregistered (for instance because it does not
     * exist)
     */
    bool unregisterVariable(const char* name);

    /** @brief Clear/unregister all the registered variables.
     */
    void clearVariables();

    /** @brief Checks if a variable has been registered.
     * @param name Name of the variable
     * @return true if already exist a variable with the given name, false
     * otherwise.
     */
    bool isVariable(const char* name);

    /** @brief Returns a variable value.
     * @param name Name of the variable
     * @return The Value of the variable, or 0.0 if the variable cannot be
     * found.
     */
    float variable(const char* name);

    /** @brief Solve a math expression.
     * @param eq Math expression to solve.
     * @param error true if the expression evaluation has failed, false
     * otherwise. NULL if not errors feedback is required.
     * @return Expression value, 0.0 if the evaluation failed (it will be
     * reported by terminal).
     */
    float solve(const char* eq, bool *error=NULL);

protected:
    /** @brief Register the default variables.
     *
     * After the libmatheval implementation, variables like \f$ e \f$ or
     * \f$ \pi \f$ are defined out of the box, but the function is retained just
     * in case.
     */
    virtual void defaultVariables();

private:
    /// Registered variables
    std::map<std::string, float> _variables;
};   // class Tokenizer

}   // namespaces

#endif // TOKENIZER_H_INCLUDED
