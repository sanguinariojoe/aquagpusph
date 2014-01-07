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

/** @class LinearInterpolation LinearInterpolation.h CalcServer/Movements/LinearInterpolation.h
 * @brief Data file linear interpolator. With a provided data file this class can perform data
 * interpolation using the time field (assuming that it is placed in the first column).
 */
class LinearInterpolation
{
public:
	/** Constructor.
	 * @param dataFile Data file path.
	 * @note Data file can be omissed at construction, but ensure yourself
	 * to provide it later.
	 */
	LinearInterpolation(const char *dataFile=NULL);

	/** Destructor.
	 */
	~LinearInterpolation();

	/** Update data.
	 * @param t Time.
	 * @return Data array. First component may be t.
	 */
	std::deque<float> update(float t);

	/** Get number of data fields.
	 * @return Number of data fields.
	 */
	unsigned int nFields(){return mData.size();}

	/** Get data fields.
	 * @return Active data. First component may be t.
	 */
	std::deque<float> data(){return mData;}

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
	/// Active time
	float _time;
	/// Active data array
	std::deque<float> mData;
	/// Previous time into the file.
	std::deque<float> mPrevData;
	/// Next time into the file.
	std::deque<float> mNextData;
};

}}} // namespace

#endif // LINEARINTERPOLATION_H_INCLUDED
