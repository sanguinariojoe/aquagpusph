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
     * @param p Points imposed.
     * @param dp Derivatives imposed.
     */
    Poly(std::deque<vec2> p, std::deque<vec2> dp);

    /** Destructor.
     */
    ~Poly();

    /** Get the imposed points.
     * @return Imposed points array.
     */
    std::deque<vec2> points(){return _p;}

    /** Get the imposed derivatives.
     * @return Imposed derivatives array.
     */
    std::deque<vec2> derivatives(){return _dp;}

    /** Get the resultant coefficients.
     * @return Coefficients.
     */
    std::deque<float> coefficients(){return _k;}

    /** Evaluates the curve in the desired x coordinate.
     * @return Curve y coordinate.
     */
    float evaluate(float x);

    /** Evaluates the curve derivative in the desired x coordinate.
     * @return Curve derivative \f$ \frac{dy}{dx} \f$.
     */
    float derivate(float x);

private:
    /** Compute the coefficients
     */
    void compute();

    /// Starting coordinate
    float _x0;
    /// Points imposed
    std::deque<vec2> _p;
    /// Derivatives imposed
    std::deque<vec2> _dp;
    /// Resultant coefficients
    std::deque<float> _k;
};

/** @class C1Interpolation C1Interpolation.h CalcServer/Movements/C1Interpolation.h
 * @brief Data file fields linear interpolator. It will read the data file
 * and perform data interpolation using the time field (assuming that it is
 * placed in the first column).
 */
class C1Interpolation
{
public:
	/** Constructor.
	 * @param data_file Data file path.
	 * @note Data file can be omissed at the construction, but ensure yourself
	 * to provide it later.
	 */
	C1Interpolation(const char *data_file=NULL);

	/** Destructor.
	 */
	~C1Interpolation();

	/** Update data for the new time instant.
	 * @param t Desired time instant.
	 * @return Data array. The first component is the time.
	 */
	std::deque<float> update(float t);

	/** Get the number of data fields.
	 * @return Number of data fields.
	 */
	unsigned int nFields(){return _data.size();}

	/** Get the data fields.
	 * @return Data array. The first component is the time.
	 */
	std::deque<float> data(){return _data;}

	/** Get the data fields derivative with respect to the time.
	 * @return Data fields derivative. The first component is the time (non
     * derivated).
	 */
	std::deque<float> derivative();

	/** Set the data file
	 * @param data_file Data file path.
	 * @return true if file was opened ok, false otherwise.
	 * @note Seek point will be moved to the last time selected in the last
	 * update calling, or \f$ t = 0 \f$ s if update has not been called yet.
	 */
	bool open(const char *data_file);

private:
	/** Reads a line of the file.
	 * @return Data array. If a bad formated line or EOF is reached, clear
	 * data array will be sent.
	 */
	std::deque<float> readLine();

	/// Data file
	FILE *_data_file;
	/// Last requested time
	float _time;
	/// Computed curves for each field
	std::deque<Poly*> _poly;
	/// Maximum time where the curve still becomes valid
	float _poly_time;
	/// Interpolated output data
	std::deque<float> _data;
};

}}} // namespace

#endif // C1Interpolation_H_INCLUDED
