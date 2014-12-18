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
 * @brief Ghost particles boundary element helper functions.
 * (See GhostParticles.cl for details)
 */

#ifdef HAVE_3D
    /** @brief Compute the point projected over the wall.
     * 
     * It is useful to mirror the particle or to check if the particle is in
     * the wall bounds
     * @param pos Particle position.
     * @param p1 corner of the considered wall.
     * @param n normal of the considered wall.
     * @return Reflection point.
     * @note Consider using the macro #wallProjection instead of this method
     * directly.
     */
    vec _wallProjection(vec pos, vec p1, vec n){
        vec p = pos - p1;
        vec d = dot(p, n) * n;
        return pos - d;
    }

    /** @def wallProjection
     * @brief Compute the point projected over the wall.
     * 
     * It is useful to mirror the particle or to check if the particle is in
     * the wall bounds.
     */
    #define wallProjection(pos) _wallProjection(pos, p1_w, n_w)

    /** @brief Check if a point is into the wall bounds.
     * @param pos Particle position.
     * @param p1 1st vertex of the wall.
     * @param p2 2nd vertex of the wall.
     * @param p3 3rd vertex of the wall.
     * @param p4 4th vertex of the wall.
     * @return true if point is on wall bounds, false otherwise.
     * @warning The method is assuming that the point is already on the wall
     * plane.
     * @note Consider using the macro #isOnWallBounds instead of this method
     * directly.
     */
    bool _isOnWallBounds(vec pos, vec p1, vec p2, vec p3, vec p4){
        // Test if the point is in the 1->2->3 triangle
        float u = dot((pos - p2), (p1 - p2)) / dot((p1 - p2), (p1 - p2));
        float v = dot((pos - p2), (p3 - p2)) / dot((p3 - p2), (p3 - p2));
        if((u + v >= 0.f) && (u + v <= 1.f))
            return true;
        // Before test the other triangle (1->3->4) we must assert that it is
        // not a triangular wall.
        const float lv =  dot((p1 - p4), (p1 - p4));
        if(!lv)
            return false;
        // Test if the point is in the 1->3->4 triangle
        u = dot((pos - p4), (p3 - p4)) / dot((p1 - p4), (p3 - p4));
        v = dot((pos - p4), (p1 - p4)) / lv;
        if((u + v >= 0.f) && (u + v <= 1.f))
            return true;
        // Tes failed
        return false;
    }

    /** @def isOnWallBounds
     * @brief Check if a point is into the wall bounds.
     */
    #define isOnWallBounds(pos) _isOnWallBounds(pos, p1_w, p2_w, p3_w, p4_w)

    /** @brief Compute a wall point velocity.
     * @param pos Wall point.
     * @param p1 1st vertex of the wall.
     * @param p2 2nd vertex of the wall.
     * @param p3 3rd vertex of the wall.
     * @param p4 4th vertex of the wall.
     * @param n normal of the considered wall.
     * @param v1 1st vertex velocity.
     * @param v2 2nd vertex velocity.
     * @param v3 3rd vertex velocity.
     * @param v4 4th vertex velocity.
     * @return Velocity at wall point.
     * @remarks The method is assuming that the point is already on the wall
     * plane.
     * @note Phong approach will be applied. To do it we launch a ray from p1
     * to pos, getting the point on the opposite edge and interpolating the
     * value into it. Then a linear interpolation is performed aver the
     * initially launched ray.
     * @note Consider using the macro #wallVelocity instead of this method
     * directly.
     */
    vec _wallVelocity(vec pos, vec p1, vec p2, vec p3, vec p4, vec n,
                      vec v1, vec v2, vec v3, vec v4){
        // Find the opposite edge where the point will be located
        // ======================================================
        vec d = pos - p1;
        float ld1 = fast_length(d);
        vec d /= ld1;
        float ld2 = dot(p3 - pos, d);
        vec p = pos + ld2 * d;
        // We may try to find the opposite point on the edge 3->2, due to
        // triangular faces will accept it ever.
        vec td = p2 - p3;
        vec vv = v2;
        if(dot(td, p - p3) < 0.f)
            // We failed because the point is on the edge 4->3.
            td = p4 - p3;
            vv = v4;
        float ltd = length(td);
        // Compute the point on the opposite edge
        // ======================================
        vec n_e = cross(td / ltd, n);
        float ld = dot(n_e, p3 - pos) / dot(n_e, d);
        p = pos + ld * d;
        // Interpolate the value on the edge point
        // =======================================
        float lt = fast_length(p - p3);
        float f = lt / ltd;
        vv = f * vv + (1.f - f) * v3;
        // Interpolate the value on the launched ray
        // =========================================
        f = ld1 / fast_length(p - p1);
        return f * vv + (1.f - f) * v1;
    }

    /** @def wallVelocity(pos)
     * @brief Compute a wall point velocity.
     */
    #define wallVelocity(pos) _wallVelocity(pos, p1_w, p2_w, p3_w, p4_w, n_w, v1_w, v2_w, v3_w, v4_w)

#else
    /** @brief Compute the point projected over the wall.
     * 
     * It is useful to mirror the particle or to check if the particle is in
     * the wall bounds
     * @param pos Particle position.
     * @param p1 corner of the considered wall.
     * @param n normal of the considered wall.
     * @return Reflection point.
     * @note Consider using the macro #wallProjection instead of this method
     * directly.
     */
    vec _wallProjection(vec pos, vec p1, vec n){
        const vec p = pos - p1;
        const vec d = dot(p, n) * n;
        return pos - d;
    }

    /** @def wallProjection
     * @brief Compute the point projected over the wall.
     * 
     * It is useful to mirror the particle or to check if the particle is in
     * the wall bounds.
     */
    #define wallProjection(pos) _wallProjection(pos, p1_w, n_w)

    /** @brief Check if a point is into the wall bounds.
     * @param pos Particle position.
     * @param p1 1st vertex of the wall.
     * @param p2 2nd vertex of the wall.
     * @return true if point is on wall bounds, false otherwise.
     * @warning The method is assuming that the point is already on the wall
     * plane.
     * @note Consider using the macro #isOnWallBounds instead of this method
     * directly.
     */
    bool _isOnWallBounds(vec pos, vec p1, vec p2){
        vec p = pos - p1;
        vec d = p2 - p1;
        float l2 = dot(p, d);
        if((l2 < 0.f) || (l2 > dot(d,d)))
            return false;
        return true;
    }

    /** @def isOnWallBounds
     * @brief Check if a point is into the wall bounds.
     */
    #define isOnWallBounds(pos) _isOnWallBounds(pos, p1_w, p2_w)

    /** @brief Compute a wall point velocity.
     * @param pos Wall point.
     * @param p1 1st vertex of the wall.
     * @param p2 2nd vertex of the wall.
     * @param v1 1st vertex velocity.
     * @param v2 2nd vertex velocity.
     * @return Velocity at wall point.
     * @remarks The method is assuming that the point is already on the wall
     * plane.
     * @note Consider using the macro #wallVelocity instead of this method
     * directly.
     */
    vec _wallVelocity(vec pos, vec p1, vec p2, vec v1, vec v2){
        vec p = pos - p1;
        vec d = p2 - p1;
        float f = length(p) / length(d);
        return f * v2 + (1.f - f) * v1;
    }

    /** @def wallVelocity(pos)
     * @brief Compute a wall point velocity.
     */
    #define wallVelocity(pos) _wallVelocity(pos, p1_w, p2_w, v1_w, v2_w)

#endif
