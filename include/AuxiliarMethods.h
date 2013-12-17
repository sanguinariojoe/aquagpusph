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

#ifndef AUXILIARMETHODS_H_INCLUDED
#define AUXILIARMETHODS_H_INCLUDED

#include <sphPrerequisites.h>

namespace Aqua{

/** Returns if a key press event has been registered.
 * @return 0 if no keys have been pressed, 1 otherwise.
 */
int isKeyPressed();

/** Compute a value which, being power of two, is greather or equal than n.
 * @param n Number which may be not power of two.
 * @return next value which is power of two.
 */
unsigned int nextPowerOf2( unsigned int n );

/** Compute if a value is power of 2.
 * @param x Value to test.
 * @return 1 if it is a power of two, 0 otherwise.
 */
unsigned int isPowerOf2(unsigned int x);

/** Compute a value, which being divisible by divisor, is greather or equal
 * than n.
 * @param n Number to round up.
 * @param divisor Divisor.
 * @return Rounded up number.
 */
unsigned int roundUp(unsigned int n, unsigned int divisor);

/** Load an OpenCL kernel from a file.
 * @param kernel The output kernel identifier.
 * @param program The output program identifier.
 * @param context Context where the program must be loaded.
 * @param device Device that must use the kernel to compute.
 * @param path Path of the kernel file.
 * @param entry_point Method into the kernel which will be called.
 * @param flags Compilation flags.
 * @param header Source code which will be inserted at the start of the readed
 * source code.
 * @return The work group size, 0 if errors happened.
 * @note Folder of the source file, and -cl-fast-relaxed-math flags will be
 * appended.
 */
size_t loadKernelFromFile(cl_kernel* kernel, cl_program* program,
                          cl_context context, cl_device_id device,
	                      const char* path, const char* entry_point,
	                      const char* flags, const char* header=NULL);

/** Gets the folder path which contains the file file_path.
 * @param file_path The file path.
 * @return The folder path.
 */
const char* getFolderFromFilePath(const char* file_path);

/** Gets the file extension.
 * @param file_path The file path.
 * @return Extension of the file.
 */
const char* getExtensionFromFilePath(const char* file_path);

/** Method that returns if a file exist on the system
 * @param file_name The file path.
 * @return 0 if the file can not be found in the system, 1 otherwise.
 */
int isFile(const char* file_name);

/** Load a file, and returns it as characters array.
 * @param source_code Readed file content.
 * @param file_name The file path.
 * @return Length of the source code array.
 * @note If SourceCode is a NULL pointer, only the length of the
 * source code will be returned, but any array will be built.
 * @warning Be sure that SourceCode is an array with, at least,
 * "length" + 1 bytes allocated space.
 */
size_t readFile(char* source_code, const char* file_name);

/** Send an argument to an OpenCL kernel.
 * @param kernel Kernel that must receive the argument.
 * @param index Index of the argument into the kernel.
 * @param size Memory size of the argument.
 * @param ptr Pointer to the argument.
 * @return 0 if the argument is succesfully sent, 1 otherwise.
 */
int sendArgument(cl_kernel kernel, int index, size_t size, void* ptr);

/** Compute the maximum local work size allowed by the device.
 * @param n Amount of data to operate in kernel (desired threats to launch).
 * @param queue Command queue.
 * @return The local work size, 0 if it is not possible to find a valid value.
 */
size_t getLocalWorkSize(cl_uint n, cl_command_queue queue);

/** Compute the global work size needed.
 * @param n Amount of data to operate in kernel (desired threats to launch).
 * @param LocalWorkSize The local work size which will be applied.
 * @return The global work size.
 */
size_t getGlobalWorkSize(cl_uint n, size_t local_work_size);

/** Gets the minimum of two values.
 * @param a First value.
 * @param b Second value.
 * @return Minimum value.
 */
template <typename T> inline T min(T a, T b){return (a>b)?b:a;}

/** Gets the maximum of two values.
 * @param a First value.
 * @param b Second value.
 * @return Maximum value.
 */
template <typename T> inline T max(T a, T b){return (a<b)?b:a;}

/** Clamps a value between the bounds.
 * @param x Value to adjust into the bounds.
 * @param a Minimum value.
 * @param b Maximum value.
 * @return Clamped value.
 */
inline float clamp(float x, float a, float b){return x < a ? a : (x > b ? b : x);}

/** Return a zeroes vector.
 * @return zeroes vector.
 */
vec Vzero();

/** Return the x direction unit vector.
 * @return x direction unit vector.
 */
vec Vx();

/** Return the y direction unit vector.
 * @return y direction unit vector.
 */
vec Vy();

#ifdef HAVE_3D
/** Return the z direction unit vector.
 * @return z direction unit vector.
 */
vec Vz();
#endif

/** Multiply a vector by a scalar.
 * @param n Value to multiply.
 * @param v Vector to operate.
 * @return product value.
 */
vec mult(float n, vec v);

/** Adding operation.
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return a+b.
 */
vec add(vec a, vec b);

/** Substracting operation.
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return a-b.
 */
vec sub(vec a, vec b);

/** Mean product.
 * @param a Vector to operate.
 * @param b Vector to operate.
 * @return Scalar product value.
 */
float dot(vec a, vec b);

/** Compute the vector length.
 * @param v Input vector.
 * @return Vector length.
 */
float length(vec v);

/** Compute a normalized vector copy.
 * @param v Vector to normalize.
 * @return Normalized copy.
 */
vec normalize(vec v);

#ifdef HAVE_3D
	/** Cross product.
	 * @param a Vector to operate.
	 * @param b Vector to operate.
	 * @return Cross product value.
	 */
	vec cross(vec a, vec b);
#endif

}   // namespace

#endif // AUXILIARMETHODS_H_INCLUDED
