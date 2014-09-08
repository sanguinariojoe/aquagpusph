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
 * @brief Linear data interpolation.
 * (See Aqua::CalcServer::Movement::LinearInterpolation for details)
 */

#ifndef LINEARINTERPOLATION_H_INCLUDED
#define LINEARINTERPOLATION_H_INCLUDED

#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class LinearInterpolation LinearInterpolation.h CalcServer/Movements/LinearInterpolation.h
 * @brief Data file fields linear interpolation.
 *
 * The data read from the file will be interpolated with a polynomial function
 * such that \f$ f \left( t \right) =
 *   \left(
 *     f \left( t_{n} \right) - f \left( t_{n-1} \right)
 *   \right)
 *   \frac{
 *     t - t_{n-1}
 *   }{
 *     t_{n} - t_{n-1}
 *   } \f$ ,
 * where \f$ t_{n-1} < t < t_{n} \f$
 *
 * @warning It is strongly recommended to use
 * Aqua::CalcServer::Movement::C1Interpolation instead of this data
 * interpolation tool.
 * @see Aqua::CalcServer::Movement::C1Interpolation
 */
class LinearInterpolation
{
public:
    /** @brief Constructor.
     * @param data_file Data file path.
     * @note Data file can be omitted at the construction, providing it later.
     */
    LinearInterpolation(const char *data_file=NULL);

    /// Destructor.
    ~LinearInterpolation();

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
    /// Last requested time time
    float _time;
    /// Data array for time _time
    std::deque<float> _data;
    /// Previous time into the file.
    std::deque<float> _prev_data;
    /// Next time into the file.
    std::deque<float> _next_data;
};

}}} // namespace

#endif // LINEARINTERPOLATION_H_INCLUDED
