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

#ifndef M_PI
	#define M_PI 3.14159265359
#endif
#ifndef M_E
	#define M_E 2.71828182845
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>

using namespace std;

/** \class Tokenizer Tokenizer.h Tokenizer/Tokenizer.h
 * @brief Tool to evaluate expressions.
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

	/** Unregister a variable, so it will be undefined.
	 * @param name Name of the variable.
	 * @return true if the variable has been unregistered, false if the
	 * variable cannot be unregistered (for instance because it does not
     * exist)
	 */
	bool unregisterVariable(const char* name);

	/** Returns if a variable has been registered.
	 * @param name Name of the desired variable
	 * @param whole if true whole word must match with a variable name,
	 * otherwise partial matching will be considered valid as well.
	 * @param readed Output readed characters.
	 * @return true if exist a variable with the given name, false otherwise.
	 */
	bool isVariable(const char* name, bool whole=true, unsigned int *readed=0);

	/** Returns a varaible value.
	 * @param name Name of the desired variable
	 * @param whole if true whole word must match with a variable name,
	 * otherwise partial matching will be considered valid as well.
	 * @param readed Output readed characters.
	 * @return Value of a variable, or 0.f if the variable cannot be found.
	 */
	float variable(const char* name, bool whole=true, unsigned int *readed=0);

	/** Solve the expresion.
	 * @param eq Expression to solve.
	 * @param value Output computed value.
	 * @return True if all gone right, false otherwise.
	 */
	bool solve(const char* eq, float *value);

protected:
	/** Solve expresion.
	 * @param eq Expression to solve.
	 * @param value Output value.
	 * @return True if all gone right, false otherwise
	 * false otherwise.
	 */
	bool solver(const char* eq, float *value);

	/** Format the sentence, removing the spaces and including parentheses.
	 * @param eq Expression to format.
	 * @return Formated sentence, NULL if errors are detected.
	 */
	const char* format(const char* eq);

    /** Removes the spaces into a sentence.
     * @param eq Expression to format.
	 * @return Formated sentence, NULL if errors are detected.
     */
    const char* removeSpaces(const char* eq);

    /** Insert parentheses in order to get right operators solving sort.
     * @param eq Expression to format.
	 * @return Formated sentence, NULL if errors are detected.
     */
    const char* sortSentence(const char* eq);

    /** Accumulate the terms surrounded by parentheses.
     * @param eq Expression to format.
     * @param n Operator position.
	 * @return Formated sentence, NULL if errors are detected.
     */
    const char* joinTerms(const char* eq, int n);

	/** Find a variable.
	 * @param name Name of the variable.
	 * @param name Name of the desired variable
	 * @param whole if true whole word must match with a variable name,
	 * otherwise partial matching will be considered valid as well.
	 * @param readed Output readed characters.
	 * @return Position of the variable, -1 if the variable cannot be found.
	 */
	int findVariable(const char* name, bool whole=true, unsigned int *readed=0);

	/** Analyze if the expression is a number
	 * @param word Word to analyze if it is a number.
	 * @param whole if true whole word must match with a variable name,
	 * otherwise partial matching will be considered valid as well.
	 * @param readed Output readed characters.
	 * @return true if the expression is a number, false otherwise.
	 */
	bool isFloat(const char* word, bool whole=true, unsigned int *readed=0);

	/** Operate.
	 * @param leftOp Left operator.
	 * @param op Operation to apply.
	 * @param value Output result value.
	 * @return True if all gone right, false otherwise.
	 * @remarks First character must be a valid operator:
	 *   - +
	 *   - -
	 *   - *
	 *   - /
	 *   - ^
	 */
	bool operate(float leftOp, const char* op, float *value);

	/** Solve function.
	 * @param func function to solve.
	 * @param value Output value.
	 * @param readed Output readed characters.
	 * @return true if all gone right, false otherwise.
	 */
	bool solveFunction(const char* func, float *value, unsigned int *readed=0);

    /** Returns the number of arguments between parentheses.
     * @param func function to solve (with the parentheses).
     * @return Number of arguments passed.
     */
    unsigned int getNumberOfArguments(const char* func);

    /** Returns a selected argument of a list.
     * @param func function to solve (with the parentheses).
     * @param id index of the argument (0 is the first argument).
     * @return Argument selected, NULL if the argument cannot be found.
     */
    const char* getArgument(const char* func, unsigned int id);

    /** Solve a expression between parentheses.
     * @param func function to solve (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveParentheses(const char* func, float *value, unsigned int *readed=0);

    /** Returns the arguments between the parentheses.
     * @param eq Expression that start with parentheses
     * @param readed Output readed characters.
     * @return Argument between the parentheses, NULL if it is not a valid
     * expression.
     */
    const char* getParenthesesArg(const char* eq, unsigned int *readed=0);

    /** Solve negative conversion.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveNegative(const char* func, float *value, unsigned int *readed=0);

    /** Converts degrees to radians.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveRad(const char* func, float *value, unsigned int *readed=0);

    /** Converts radians to degrees.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveDeg(const char* func, float *value, unsigned int *readed=0);

    /** Solve sin function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveSin(const char* func, float *value, unsigned int *readed=0);

    /** Solve cos function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveCos(const char* func, float *value, unsigned int *readed=0);

    /** Solve tan function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveTan(const char* func, float *value, unsigned int *readed=0);

    /** Solve asinf function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveAsin(const char* func, float *value, unsigned int *readed=0);

    /** Solve acosf function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveAcos(const char* func, float *value, unsigned int *readed=0);

    /** Solve atanf function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveAtan(const char* func, float *value, unsigned int *readed=0);

    /** Solve atan2f function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveAtan2(const char* func, float *value, unsigned int *readed=0);

    /** Solve exp function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveExp(const char* func, float *value, unsigned int *readed=0);

    /** Solve log function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveLog(const char* func, float *value, unsigned int *readed=0);

    /** Solve log10 function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveLog10(const char* func, float *value, unsigned int *readed=0);

    /** Solve pow function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solvePow(const char* func, float *value, unsigned int *readed=0);

    /** Solve sqrt function.
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveSqrt(const char* func, float *value, unsigned int *readed=0);

    /** Solve clamp function (adjust a value to their limits).
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveClamp(const char* func, float *value, unsigned int *readed=0);

    /** Solve min function (Maximum value of twice values).
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveMin(const char* func, float *value, unsigned int *readed=0);

    /** Solve max function (Maximum value of twice values).
     * @param func function argument (with the parentheses).
     * @param value Output value.
     * @param readed Output readed characters
     * @return True if all gone right, false otherwise
     */
    bool solveMax(const char* func, float *value, unsigned int *readed=0);

private:
	/// Registered variable names
	std::vector<char*> regNames;
	/// Registered variable values (float value type is considered)
	std::vector<float> regValues;

	/** Get the minimum of two values.
	 * @param a First value.
	 * @param b Second value.
	 * @return Minimum value.
	 */
	float min(float a, float b){return (a>b)?b:a;}

	/** Get the maximum of two values.
	 * @param a First value.
	 * @param b Second value.
	 * @return Maximum value.
	 */
	float max(float a, float b){return (a<b)?b:a;}

	/** Clamps a value between the bounds.
	 * @param x Value to adjust into the bounds.
	 * @param a Minimum value.
	 * @param b Maximum value.
	 * @return Clamped value.
	 */
	float clamp(float x, float a, float b){return x < a ? a : (x > b ? b : x);}
};   // class Tokenizer

#endif // TOKENIZER_H_INCLUDED
