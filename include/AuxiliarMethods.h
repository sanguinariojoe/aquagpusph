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
 * @brief Set of auxiliar functions.
 */

#ifndef AUXILIARMETHODS_H_INCLUDED
#define AUXILIARMETHODS_H_INCLUDED

#include <sphPrerequisites.h>
#include <string>

namespace Aqua{

/** @brief Check wkether a key press event has been triggered.
 * 
 * @return false if no keys have been pressed, true otherwise.
 */
const bool isKeyPressed();

/** @brief Check if a string ends with an specific suffix
 * 
 * @param str String to be checked
 * @param suffix Suffix to be looked for
 * @return true if #str ends with #suffix, false otherwise.
 */
const bool hasSuffix(const std::string &str, const std::string &suffix);

/** @brief Replace all substring occurrences by another substring
 * 
 * @param str String to be modified
 * @param search Substring to be replaced
 * @param replace Replacement substring
 */
void replaceAll(std::string &str,
                const std::string &search,
                const std::string &replace);

/** @brief Replace all substring occurrences by another substring
 * 
 * @param str String to be modified
 * @param search Substring to be replaced
 * @param replace Replacement substring
 * @return Modified string
 */
std::string replaceAllCopy(std::string str,
                           const std::string &search,
                           const std::string &replace);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix.
 * 
 * @param s String to become trimmed
 */
void ltrim(std::string &s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 * 
 * @param s String to become trimmed
 */
void rtrim(std::string &s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix and suffix.
 * 
 * @param s String to become trimmed
 */
void trim(std::string &s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 * 
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string ltrimCopy(std::string s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 * 
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string rtrimCopy(std::string s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix and suffix.
 * 
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string trimCopy(std::string s);

/** @brief Transform a xxd exported file into a C++ string.
 *
 * xdd can be used to embed external files into the program source code, by
 * means of an include statement. However, the data is exported as a non-null
 * terminated char array and its length.
 * This method is appending the null character, and building a C++ string
 * 
 * @param arr C-Like chars array
 * @param len C-Like chars array length
 * @return C++ string
 */
std::string xxd2string(const unsigned char* arr, const unsigned int &len);

/** @brief Convert a string to lower case
 */
void toLower(std::string &str);

/** @brief Convert a string to lower case
 * @return Modified string
 */
std::string toLowerCopy(std::string str);

/** @brief Check if a number is power of 2.
 *
 * @param x Value to test.
 * @return true if it is a power of two, false otherwise.
 */
inline const bool isPowerOf2(const unsigned int &x) {
    return x && ((x & (x - 1)) == 0);
}

/** @brief Next number which is power of 2.
 *
 * Compute a power of value greater or equal than @paramname{x}.
 *
 * @param x Number to round up to a power of two.
 * @return Next value which is power of two (it can be the same @paramname{x}
 * input value).
 */
const unsigned int nextPowerOf2(const unsigned int &x);

/** @brief Rounded up value which is divisible by @paramname{divisor}.
 *
 * Compute a value which is divisible by @paramname{divisor} and greater
 * or equal than @paramname{x}.
 *
 * @param x Number to round up.
 * @param divisor Divisor.
 * @return Rounded up number.
 */
const unsigned int roundUp(const unsigned int &x, const unsigned int &divisor);

/** @brief Round an float value to an integer one.
 * @param n Number to round.
 * @return The closest integer to @paramname{n}.
 */
const int round(const float &n);

/** @brief Gets the folder path which contains the file @paramname{file_path}.
 *
 * @param file_path The file path.
 * @return The folder.
 */
const std::string getFolderFromFilePath(const std::string &file_path);

/** @brief Gets the file name of the path @paramname{file_path}.
 *
 * @param file_path The file path.
 * @return The file name.
 */
const std::string getFileNameFromFilePath(const std::string &file_path);

/** @brief Gets the file extension.
 *
 * Get the file extension from the full file path @paramname{file_path}.
 *
 * @param file_path The file path.
 * @return Extension of the file.
 */
const std::string getExtensionFromFilePath(const std::string &file_path);

/** @brief Check if the file @paramname{file_name} exist on the system.
 *
 * @param file_name The file path.
 * @return false if the file can not be found in the system, true otherwise.
 */
const bool isFile(const std::string &file_name);

/** @brief Check if the path @paramname{path} is a relative or an absolute one.
 *
 * @param path The path.
 * @return true if it is a relative path, false otherwise.
 */
const bool isRelativePath(const std::string &path);

/** @brief Compute the maximum local work size allowed by a device.
 *
 * @param queue Command queue.
 * @return The local work size, 0 if it is not possible to find a valid value.
 */
const size_t getLocalWorkSize(const cl_command_queue &queue);

/** @brief Compute the global work size needed to compute @paramname{n} threads.
 *
 * @param n Amount of data to operate in kernel (aiming threads to launch).
 * @param local_work_size The local work size which will be applied.
 * @return The required global work size.
 *
 * @see roundUp()
 */
inline const size_t getGlobalWorkSize(const cl_uint &n,
                                      const size_t &local_work_size)
{
    return roundUp(n, local_work_size);
}

/** @brief Gets the minimum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return Minimum value.
 */
template <typename T> inline const T min(const T &a, const T &b)
{
    return (a > b) ? b : a;
}

/** @brief Gets the maximum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return Maximum value.
 */
template <typename T> inline const T max(const T &a, const T &b)
{
    return (a < b) ? b : a;
}

/** @brief Clamps a value between the bounds.
 *
 * @param x Value to adjust into the bounds.
 * @param a Minimum value.
 * @param b Maximum value.
 * @return Clamped value.
 */
inline const float clamp(const float &x, const float &a, const float &b)
{
    return x < a ? a : (x > b ? b : x);
}

/** @brief Return a null vector.
 *
 * @return zeroes vector.
 */
inline const vec Vzero() { return VEC_ZERO; }

/** @brief Return the x direction unit vector.
 *
 * @return x direction unit vector.
 */
inline const vec Vx() { return VEC_X; }

/** @brief Return the y direction unit vector.
 *
 * @return y direction unit vector.
 */
inline const vec Vy() { return VEC_Y; }

#ifdef HAVE_3D
/** @brief Return the z direction unit vector.
 *
 * @remarks Only available in the 3D version.
 * @return z direction unit vector.
 */
inline const vec Vz() { return VEC_Z; }
#endif

/** @brief Multiply a vector by a scalar.
 *
 * @param n Scalar to operate.
 * @param v Vector to operate.
 * @return @paramname{n} \f$ \cdot \f$ @paramname{v} Resulting vector.
 */
const vec mult(const float &n, const vec &v);

/** @brief Adding operation.
 *
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} + @paramname{b}.
 */
const vec add(const vec &a, const vec &b);

/** @brief Subtracting operation.
 *
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} - @paramname{b}.
 */
const vec sub(const vec &a, const vec &b);

/** @brief Inner product.
 *
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} \f$ \cdot \f$ @paramname{b} scalar product value.
 */
const float dot(const vec &a, const vec &b);

/** @brief Compute the vector length.
 *
 * @param v Input vector.
 * @return \f$ \vert \f$ @paramname{v} \f$ \vert \f$ vector length.
 */
const float length(const vec &v);

/** @brief Compute a normalized vector copy (such that length() will return 1.0.
 *
 * @param v Vector to normalize.
 * @return Normalized copy of the vector.
 *
 * @see length()
 */
const vec normalize(const vec &v);

#ifdef HAVE_3D
/** @brief Cross product.
 *
 * @remarks Only available in the 3D version.
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} \f$ \times \f$ @paramname{b} crossed product vector.
 */
const vec cross(const vec &a, const vec &b);
#endif

/** @brief Get the number of digits of an integer decimal text representation.
 *
 * @param number Number from which the number of digits should be computed.
 * @return Number of digits.
 */
const unsigned int numberOfDigits(const unsigned int number);

}   // namespace

#endif // AUXILIARMETHODS_H_INCLUDED
