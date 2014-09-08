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
 * @brief Continuous C1 data interpolation.
 * (See Aqua::CalcServer::Movement::C1Interpolation for details)
 */

#ifndef C1Interpolation_H_INCLUDED
#define C1Interpolation_H_INCLUDED

#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class Poly C1Interpolation.h CalcServer/Movements/C1Interpolation.h
 * @brief Polynomial description of a curve.
 */
class Poly
{
public:
    /** @brief Constructor.
     * @param p List of points \f$ y(x) \f$ that the curve must cross.
     * @param dp List of derivatives \f$ \frac{d y}{d x}(x) \f$ imposed to the
     * polynomial function.
     */
    Poly(std::deque<vec2> p, std::deque<vec2> dp);

    /// Destructor.
    ~Poly();

    /** @brief Get the imposed points \f$ y(x) \f$.
     * @return Imposed points list.
     */
    std::deque<vec2> points(){return _p;}

    /** @brief Get the imposed derivatives \f$ \frac{d y}{d x}(x) \f$.
     * @return Imposed derivatives list.
     */
    std::deque<vec2> derivatives(){return _dp;}

    /** @brief Get the polynomial function coefficients.
     * @return Coefficients.
     */
    std::deque<float> coefficients(){return _k;}

    /** @brief Evaluates the curve in the desired \f$ x \f$ coordinate.
     * @param x X coordinate to evaluate the function.
     * @return Curve y coordinate.
     */
    float evaluate(float x);

    /** Evaluates the curve derivative in the desired \f$ x \f$ coordinate.
     * @param x X coordinate to evaluate the function derivative.
     * @return Curve derivative \f$ \frac{d y}{d x}(x) \f$.
     */
    float derivate(float x);

private:
    /// Compute the coefficients
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
 * @brief Data file fields continuous C1 interpolation.
 *
 * The data read from the file will be interpolated with a polynomial function
 * such that:
 *    - The curve will cross the points at \f$ t_{n-1} \f$ and \f$ t_{n} \f$
 *    - The curve will share the derivative value at \f$ t_{n-1} \f$ with the
 *      polynomial function of the time interval \f$ t_{n-2}, t_{n-1} \f$
 *
 * where \f$ t_{n-1} < t < t_{n} \f$
 *
 * @see Aqua::CalcServer::Movement::Poly
 * @see Aqua::CalcServer::Movement::LinearInterpolation
 */
class C1Interpolation
{
public:
    /** @brief Constructor.
     * @param data_file Data file path.
     * @note Data file can be omitted at the construction, providing it later.
     */
    C1Interpolation(const char *data_file=NULL);

    /// Destructor.
    ~C1Interpolation();

    /** @brief Update data for a new time instant.
     * @param t Desired time instant.
     * @return Data array. The first component is the time.
     */
    std::deque<float> update(float t);

    /** @brief Get the number of data fields.
     * @return Number of data fields.
     */
    unsigned int nFields(){return _data.size();}

    /** @brief Get the data fields.
     * @return Data array. The first component is the time.
     */
    std::deque<float> data(){return _data;}

    /** @brief Get the data fields derivative with respect to the time.
     * @return Data fields derivative. The first component is the time.
     */
    std::deque<float> derivative();

    /** @brief Set the data file
     * @param data_file Data file path.
     * @return true if file was successfully opened, false otherwise.
     * @note Seek point will be moved to the last time selected in the last
     * update calling, or \f$ t = 0 \f$ s if update has not been called yet.
     */
    bool open(const char *data_file);

private:
    /** @brief Reads a line of the file.
     * @return Data array. If a bad formatted line or EOF is reached, empty
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
