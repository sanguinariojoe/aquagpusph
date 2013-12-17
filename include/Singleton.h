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

#ifndef SINGLETON_H_INCLUDED
#define SINGLETON_H_INCLUDED

// ----------------------------------------------------------------------------
// Include standar header
// ----------------------------------------------------------------------------
#include <stdlib.h>

/// @namespace Aqua Main AQUAgpusph namespace.
namespace Aqua{

/** \class Singleton Singleton.h Singleton.h
 * \brief Simple but versatile singleton instance (Meyers singleton).
 */
template<typename T> class Singleton
{
public:
	/** Returns the instance of the class.
	 * @return singleton instance.
	 */
	static T* singleton(){return singletonPtr;}

protected:
	/** Constructor
	 */
	Singleton(){
	    singletonPtr = (T*)this;
	}
	/** Destructor
	 */
	~Singleton(){singletonPtr = NULL;}

private:
	/// Static singleton pointer store
	static T *singletonPtr;
};

template <typename T> T* Singleton <T>::singletonPtr = 0;

}   // namespaces

#endif // SINGLETON_H_INCLUDED
