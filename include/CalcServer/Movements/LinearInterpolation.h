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

#ifndef LINEARINTERPOLATION_H_INCLUDED
#define LINEARINTERPOLATION_H_INCLUDED

#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class LinearInterpolation LinearInterpolation.h CalcServer/Movements/LinearInterpolation.h
 * @brief Data file fields linear interpolator. It will read the data file
 * and perform data interpolation (C0 such that the function is contiguous)
 * using the time field (assuming that it is placed in the first column).
 * @warning It is strongly recommended to use C1Interpolation intead of this
 * data interpolator.
 */
class LinearInterpolation
{
public:
	/** Constructor.
	 * @param data_file Data file path.
	 * @note Data file can be omissed at construction, but ensure yourself
	 * to provide it later.
	 */
	LinearInterpolation(const char *data_file=NULL);

	/** Destructor.
	 */
	~LinearInterpolation();

	/** Update data.
	 * @param t Time.
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
