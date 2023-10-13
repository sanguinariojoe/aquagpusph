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
#include <vector>

namespace Aqua {

/// Returns if a key press event has been registered.
/**
 * @return 0 if no keys have been pressed, 1 otherwise.
 */
int
isKeyPressed();

/** @brief Check if a string starts with an specific prefix
 *
 * @param str String to be checked
 * @param prefix Prefix to be looked for
 * @return true if #str starts with #prefix, false otherwise.
 */
bool
hasPrefix(const std::string& str, const std::string& prefix);

/** @brief Alias for hasPrefix()
 */
inline bool
startswith(const std::string& str, const std::string& prefix)
{
    return hasPrefix(str, prefix);
}

/** @brief Check if a string ends with an specific suffix
 *
 * @param str String to be checked
 * @param suffix Suffix to be looked for
 * @return true if #str ends with #suffix, false otherwise.
 */
bool
hasSuffix(const std::string& str, const std::string& suffix);

/** @brief Alias for hasSuffix()
 */
inline bool
endswith(const std::string& str, const std::string& prefix)
{
    return hasSuffix(str, prefix);
}

/** @brief Replace all substring occurrences by another substring
 *
 * @param str String to be modified
 * @param search Substring to be replaced
 * @param replace Replacement substring
 */
void
replaceAll(std::string& str,
           const std::string& search,
           const std::string& replace);

/** @brief Replace all substring occurrences by another substring
 *
 * @param str String to be modified
 * @param search Substring to be replaced
 * @param replace Replacement substring
 * @return Modified string
 */
std::string
replaceAllCopy(std::string str, std::string search, std::string replace);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix.
 *
 * @param s String to become trimmed
 */
void
ltrim(std::string& s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 *
 * @param s String to become trimmed
 */
void
rtrim(std::string& s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix and suffix.
 *
 * @param s String to become trimmed
 */
void
trim(std::string& s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 *
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string
ltrimCopy(std::string s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string suffix.
 *
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string
rtrimCopy(std::string s);

/** @brief Remove all the blank spaces (including line breaks, tabulators...)
 * string prefix and suffix.
 *
 * @param s String to become trimmed
 * @return Trimmed string
 */
std::string
trimCopy(std::string s);

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
std::string
xxd2string(unsigned char* arr, unsigned int len);

/** @brief Convert a string to lower case
 */
void
toLower(std::string& str);

/** @brief Convert a string to lower case
 * @return Modified string
 */
std::string
toLowerCopy(std::string str);

/** @brief Set several constants into a string
 *
 * Following scape strings can be used, that will be replaced as follows:
 *   - `"{mpi_rank}"` will be replaced by the process identifier index
 *   - `"{version}"` will be replaced by the AQUAgpusph package version
 *
 * @param str The string where constants shall be set
 */
void
setStrConstants(std::string& str);

/** @brief Set several constants into a string
 *
 * Following scape strings can be used, that will be replaced as follows:
 *   - `"{mpi_rank}"` will be replaced by the process identifier index
 *   - `"{version}"` will be replaced by the AQUAgpusph package version
 *
 * @param str The string where constants shall be set
 * @return The modified string
 */
std::string
setStrConstantsCopy(std::string str);

/** @brief Split a string by a character
 * @param str String to be split
 * @param chr Splitting cahracter
 * @return List of substrings
 */
std::vector<std::string>
split(std::string str, char chr);

/** @brief Split a list of split_formulae
 *
 * Formulae can be separated either by semicolon, ';', or by simple comma, ','.
 * This function already takes care of parentheses, but it is not checking for
 * syntax errors
 *
 * @param str String within formulae
 * @return List of formulae string
 */
std::vector<std::string>
split_formulae(std::string str);

/** @brief Look for a file path which is not already taken
 *
 * Several scape strings can be used, that will be replaced as follows:
 *   - `"%d"`/`"{index}"` will be replaced by the first integer which
 *     results in a non-existing file path
 *   - `"{mpi_rank}"` will be replaced by the process identifier index
 *   - `"{version}"` will be replaced by the AQUAgpusph package version
 *
 * In case a file path cannot be obtained, a std::invalid_argument exception
 * will be raised.
 *
 * @param basename The base name of the file
 * @param i First index that will be checked. Also, the resulting file index is
 * stored in this variable.
 * @return The available file path
 * @note In case overwriting files is allowed, you can call this method with a
 * try, catching std::invalid_argument exception. If exception
 * std::invalid_argument is eventually raised, you can call
 * setStrConstantsCopy()
 */
std::string
newFilePath(const std::string& basename,
            unsigned int& i,
            unsigned int digits = 5);

/// Next number which is power of 2.
/** Compute a value which, being power of two, is greater or equal than
 * @paramname{x}.
 *
 * @param x Number to round up to a power of two.
 * @return Next value which is power of two (it can be the same input value
 * @paramname{x}).
 */
unsigned int
nextPowerOf2(unsigned int x);

/// Check if a number is power of 2.
/** Compute if a value is power of 2.
 *
 * @param x Value to test.
 * @return 1 if it is a power of two, 0 otherwise.
 */
unsigned int
isPowerOf2(unsigned int x);

/// Rounded up value which is divisible by @paramname{divisor}.
/** Compute a value, which being divisible by @paramname{divisor}, is greater
 * or equal than @paramname{x}.
 *
 * @param x Number to round up.
 * @param divisor Divisor.
 * @return Rounded up number.
 */
unsigned int
roundUp(unsigned int x, unsigned int divisor);

/** @brief Round an float value to an integer one.
 * @param n Number to round.
 * @return The closest integer to @paramname{n}.
 */
int
round(float n);

/// Gets the folder path which contains the file @paramname{file_path}.
/**
 * @param file_path The file path.
 * @return The folder.
 */
const std::string
getFolderFromFilePath(const std::string file_path);

/// Gets the file name of the path @paramname{file_path}.
/**
 * @param file_path The file path.
 * @return The file name.
 */
const std::string
getFileNameFromFilePath(const std::string file_path);

/// Gets the file extension.
/** Get the file extension from the full file path @paramname{file_path}.
 *
 * @param file_path The file path.
 * @return Extension of the file.
 */
const std::string
getExtensionFromFilePath(const std::string file_path);

/** @brief Check if the file @paramname{file_name} exist on the system.
 *
 * @param file_name The file path.
 * @return false if the file can not be found in the system, true otherwise.
 */
bool
isFile(const std::string file_name);

/** @brief Check if the path @paramname{path} is a relative or an absolute one.
 *
 * @param path The path.
 * @return true if it is a relative path, false otherwise.
 */
bool
isRelativePath(const std::string path);

/// Compute the maximum local work size allowed by a device.
/**
 * @param n Amount of data to operate in kernel (aiming threads to launch).
 * @param queue Command queue.
 * @return The local work size, 0 if it is not possible to find a valid value.
 */
size_t
getLocalWorkSize(cl_uint n, cl_command_queue queue);

/// Compute the global work size needed to compute @paramname{n} threads.
/**
 * @param n Amount of data to operate in kernel (aiming threads to launch).
 * @param local_work_size The local work size which will be applied.
 * @return The required global work size.
 *
 * @see roundUp()
 */
size_t
getGlobalWorkSize(cl_uint n, size_t local_work_size);

/// Gets the minimum of two values.
/**
 * @param a First value.
 * @param b Second value.
 * @return Minimum value.
 */
template<typename T>
inline T
min(T a, T b)
{
	return (a > b) ? b : a;
}

/// Gets the maximum of two values.
/**
 * @param a First value.
 * @param b Second value.
 * @return Maximum value.
 */
template<typename T>
inline T
max(T a, T b)
{
	return (a < b) ? b : a;
}

/// Clamps a value between the bounds.
/**
 * @param x Value to adjust into the bounds.
 * @param a Minimum value.
 * @param b Maximum value.
 * @return Clamped value.
 */
inline float
clamp(float x, float a, float b)
{
	return x < a ? a : (x > b ? b : x);
}

/// Return a null vector.
/**
 * @return zeroes vector.
 */
vec
Vzero();

/// Return the x direction unit vector.
/**
 * @return x direction unit vector.
 */
vec
Vx();

/// Return the y direction unit vector.
/**
 * @return y direction unit vector.
 */
vec
Vy();

#ifdef HAVE_3D
/// Return the z direction unit vector.
/**
 * @remarks Only available in the 3D version.
 * @return z direction unit vector.
 */
vec
Vz();
#endif

/// Multiply a vector by a scalar.
/**
 * @param n Scalar to operate.
 * @param v Vector to operate.
 * @return @paramname{n} \f$ \cdot \f$ @paramname{v} Resulting vector.
 */
vec
mult(float n, vec v);

/// Adding operation.
/**
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} + @paramname{b}.
 */
vec
add(vec a, vec b);

/// Subtracting operation.
/**
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} - @paramname{b}.
 */
vec
sub(vec a, vec b);

/// Inner product.
/**
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} \f$ \cdot \f$ @paramname{b} scalar product value.
 */
float
dot(vec a, vec b);

/// Compute the vector length.
/**
 * @param v Input vector.
 * @return \f$ \vert \f$ @paramname{v} \f$ \vert \f$ vector length.
 */
float
length(vec v);

/// Compute a normalized vector copy (such that length() will return 1.0.
/**
 * @param v Vector to normalize.
 * @return Normalized copy of the vector.
 *
 * @see length()
 */
vec
normalize(vec v);

#ifdef HAVE_3D
/// Cross product.
/**
 * @remarks Only available in the 3D version.
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return @paramname{a} \f$ \times \f$ @paramname{b} crossed product vector.
 */
vec
cross(vec a, vec b);
#endif

/// Get the number of digits of an integer decimal text representation.
/**
 * @param number Number from which the number of digits should be computed.
 * @return Number of digits.
 */
unsigned int
numberOfDigits(unsigned int number);

} // namespace

#endif // AUXILIARMETHODS_H_INCLUDED
