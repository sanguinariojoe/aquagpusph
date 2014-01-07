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

#ifndef C1Interpolation_H_INCLUDED
#define C1Interpolation_H_INCLUDED

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>
#include <unistd.h>
#include <errno.h>
#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class Poly C1Interpolation.h CalcServer/Movements/C1Interpolation.h
 * @brief Polynomial description of a curve such that points and derivatives
 * can be imposed.
 */
class Poly
{
public:
    /** Constructor.
     * @param p Points which the curve must pass through.
     * @param dp Derivatives imposed.
     */
    Poly(std::deque<vec2> p, std::deque<vec2> dp);

    /** Destructor.
     */
    ~Poly();

    /** Get the imposed points.
     * @return imposed points array.
     */
    std::deque<vec2> points(){return mP;}

    /** Get the imposed derivatives.
     * @return imposed derivatives array.
     */
    std::deque<vec2> derivatives(){return mD;}

    /** Get the resultant coefficients.
     * @return Coefficients.
     */
    std::deque<float> coefficients(){return mK;}

    /** Evaluates the curve in the desired x coordinate.
     * @return Curve y coordinate.
     */
    float evaluate(float x);

    /** Evaluates the curve derivative in the desired x coordinate.
     * @return Curve derivative \f$\frac{dy}{dx}\f$.
     */
    float derivate(float x);

private:
    /** Compute the coefficients
     */
    void compute();

    /// Starting coordinate
    float mX0;
    /// Points imposed
    std::deque<vec2> mP;
    /// Derivatives imposed
    std::deque<vec2> mD;
    /// Resultant coefficients
    std::deque<float> mK;
};

/** @class C1Interpolation C1Interpolation.h CalcServer/Movements/C1Interpolation.h
 * @brief Data file linear interpolator. With a provided data file this class can perform data
 * interpolation using the time field (assuming that it is placed in the first column).
 */
class C1Interpolation
{
public:
	/** Constructor.
	 * @param dataFile Data file path.
	 * @note Data file can be omissed at construction, but ensure yourself
	 * to provide it later.
	 */
	C1Interpolation(const char *dataFile=NULL);

	/** Destructor.
	 */
	~C1Interpolation();

	/** Update data.
	 * @param t Desired time instant.
	 * @return Data array. The first component is the time.
	 */
	std::deque<float> update(float t);

	/** Get number of data fields.
	 * @return Number of data fields.
	 */
	unsigned int nFields(){return mData.size();}

	/** Get data fields.
	 * @return Data array. The first component is the time.
	 */
	std::deque<float> data(){return mData;}

	/** Get data fields derivative.
	 * @return Data derivative. The first component is the time.
	 */
	std::deque<float> derivative();

	/** Set the data file
	 * @param dataFile Data file path.
	 * @return true if file was opened ok. \n false otherwise.
	 * @note Seek will move to last time selected with update,
	 * if any t=0s will selected.
	 */
	bool open(const char *dataFile);

private:
	/** Reads a line of file.
	 * @return Data array. If bad formated line or EOF reached,
	 * clear array will sent.
	 */
	std::deque<float> readLine();

	/// Data file
	FILE *mDataFile;
	/// Last requested time
	float _time;
	/// Last computed curves
	std::deque<Poly*> mPoly;
	/// Valid time for the computed curve
	float mPolyTime;
	/// Interpolated output data
	std::deque<float> mData;
};

}}} // namespace

#endif // C1Interpolation_H_INCLUDED
