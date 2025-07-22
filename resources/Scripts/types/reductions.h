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
 * @brief Some handy reductions for vectorial types
 */

inline float
reduce_max_vec2(const vec2 in)
{
	return max(in.x, in.y);
}

#define reduce_max_float2 reduce_max_vec2

inline float
reduce_max_vec4(const vec4 in)
{
	return max(reduce_max_vec2(in.xy), reduce_max_vec2(in.zw));
}

#define reduce_max_float4 reduce_max_vec4

inline float
reduce_max_vec8(const vec8 in)
{
	return max(reduce_max_vec4(in.s0123), reduce_max_vec4(in.s4567));
}

#define reduce_max_float8 reduce_max_vec8

inline float
reduce_max_vec16(const vec16 in)
{
	return max(reduce_max_vec8(in.s01234567), reduce_max_vec8(in.s89ABCDEF));
}

#define reduce_max_float16 reduce_max_vec16

inline float
reduce_min_vec2(const vec2 in)
{
	return min(in.x, in.y);
}

#define reduce_min_float2 reduce_min_vec2

inline float
reduce_min_vec4(const vec4 in)
{
	return min(reduce_min_vec2(in.xy), reduce_min_vec2(in.zw));
}

#define reduce_min_float4 reduce_min_vec4

inline float
reduce_min_vec8(const vec8 in)
{
	return min(reduce_min_vec4(in.s0123), reduce_min_vec4(in.s4567));
}

#define reduce_min_float8 reduce_min_vec8

inline float
reduce_min_vec16(const vec16 in)
{
	return min(reduce_min_vec8(in.s01234567), reduce_min_vec8(in.s89ABCDEF));
}

#define reduce_min_float16 reduce_min_vec16

inline float
reduce_sum_vec2(const vec2 in)
{
	return in.x + in.y;
}

#define reduce_sum_float2 reduce_sum_vec2

inline float
reduce_sum_vec4(const vec4 in)
{
	return reduce_sum_vec2(in.xy) + reduce_sum_vec2(in.zw);
}

#define reduce_sum_float4 reduce_sum_vec4

inline float
reduce_sum_vec8(const vec8 in)
{
	return reduce_sum_vec4(in.s0123) + reduce_sum_vec4(in.s4567);
}

#define reduce_sum_float8 reduce_sum_vec8

inline float
reduce_sum_vec16(const vec16 in)
{
	return reduce_sum_vec8(in.s01234567) + reduce_sum_vec8(in.s89ABCDEF);
}

#define reduce_sum_float16 reduce_sum_vec16
