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
 * @brief Singleton instance.
 * (See Aqua::Singleton for details)
 */

#ifndef SINGLETON_H_INCLUDED
#define SINGLETON_H_INCLUDED

namespace Aqua{

/** \class Singleton Singleton.h Singleton.h
 * \brief Simple but versatile singleton instance (Meyers singleton).
 *
 * The singletons are instances that can be accessed everywhere even though it
 * has not been passed as argument.
 *
 * To create a singleton instance just define the class as an inherited of this
 * one:
 *
 * @code{.cpp}
    #include <Singleton.h>
    class MyClass : public Aqua::Singleton<MyClass>
    {
        ...
    };
   @endcode
 *
 * Then you can access your class everywhere you want (please, don't forget to
 * include your class header):
 *
 * @code{.cpp}
    MyClass *c = MyClass::singleton();
   @endcode
 *
 * @warning This implementation will not support multithreading.
 */
template<typename T> class Singleton
{
public:
    /** @brief Returns the instance of the class.
     * @return singleton instance.
     */
    static T* singleton(){return _singletonPtr;}

protected:
    /// Constructor
    Singleton(){
        _singletonPtr = (T*)this;
    }
    /// Destructor
    virtual ~Singleton(){_singletonPtr = NULL;}

private:
    /// Static singleton pointer store
    static T *_singletonPtr;
};

/// Initialization of the singleton as a NULL pointer
template <typename T> T* Singleton <T>::_singletonPtr = NULL;

}   // namespaces

#endif // SINGLETON_H_INCLUDED
