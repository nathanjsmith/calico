// Copyright (c) 2014, Nathan Smith <nathanjsmith@gmail.com>
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef __CALICO__ACCELERATOR__BRUTEFORCE__HPP__
#define __CALICO__ACCELERATOR__BRUTEFORCE__HPP__

#include <calico/math.hpp>

#include <limits>

namespace calico {
namespace accelerator {

/**
    Given a mesh and a containment test, this will iterate across all surfaces
    in the mesh and test first the plane in which the surface lies for
    intersection with the provided ray and then as to whether the intersection
    point is contained inside the bounds of the surface.
*/
template <typename Float, 
          typename Mesh, 
          typename Containment, 
          typename limits=std::numeric_limits<Float> >
class BruteForce {
public:
    BruteForce(const Mesh &mesh) : _mesh(mesh) {}

    /**
        Find the intersection of the provided ray and a face in the assigned
        Mesh object by performing an intersection test against all faces in
        Mesh.

        If an intersection is found, then t, face and hit_{x,y,z} will all be
        set according to the intsection. Ray misses are denoted by setting face
        to Mesh::ray_miss_id_c. If face is not equal to Mesh::ray_miss_id_c,
        then an intersection was found. If face is equal to Mesh::ray_miss_id_c, 
        then t, hit_x, hit_y, and hit_z are undefined on output.

        @param start_x      X-component of the ray's start position (input)
        @param start_y      Y-component of the ray's start position (input)
        @param start_z      Z-component of the ray's start position (input)
        @param direction_x  X-component of the ray's unit direction (input)
        @param direction_y  Y-component of the ray's unit direction (input)
        @param direction_z  Z-component of the ray's unit direction (input)
        @param face         Face for which intersection testing will be skipped 
                            (input)

        @param face         Id of the face on which the closest intersection
                            was discovered (output)
        @param t            Distance along the path that the intersection 
                            is found (output)
        @param hit_x        X-component of the ray intersection (output)
        @param hit_y        Y-component of the ray intersection (output)
        @param hit_z        Z-component of the ray intersection (output)
    */
    void find_intersection(const Float start_x, 
                           const Float start_y, 
                           const Float start_z,
                           const Float direction_x, 
                           const Float direction_y, 
                           const Float direction_z,
                           typename Mesh::FaceId &face, Float &t, 
                           Float &hit_x, Float &hit_y, Float &hit_z)
    {
        typename Mesh::FaceId ignore_face = face;
        face = Mesh::ray_miss_id_c;
        t    = limits::infinity();

        for (typename Mesh::FaceId f = 0u; f < _mesh.size(); ++f) {

            if (f == ignore_face) {
                continue;
            }

            // If the ray lies in the plane, this should return Inf as the
            // intersection distance.  That will get filtered out in the
            // following tests.
            Float tmp_t = math::ray_plane_intersection(_mesh.normal_x(f),
                                                       _mesh.normal_y(f),
                                                       _mesh.normal_z(f),
                                                       _mesh.d(f),
                                                       start_x,
                                                       start_y,
                                                       start_z,
                                                       direction_x,
                                                       direction_y,
                                                       direction_z);

            // Only keep hits in front of the start point.  Only calculate the
            // intersection point if the hit distance is better than our
            // previous best result.
            //
            // TODO: Try considering all hits and see if it does better for
            //       branch prediction
            if (Float(0.) < tmp_t && tmp_t < t) {
                Float tmp_hit_x = start_x + (direction_x * tmp_t);
                Float tmp_hit_y = start_y + (direction_y * tmp_t);
                Float tmp_hit_z = start_z + (direction_z * tmp_t);

                if (Containment::contains(_mesh, f, tmp_hit_x, tmp_hit_y, tmp_hit_z)) {
                    hit_x = tmp_hit_x;
                    hit_y = tmp_hit_y;
                    hit_z = tmp_hit_z;
                    face  = f;
                    t     = tmp_t;
                }
            }
        }// end accelerator loop
    }

private:
    const Mesh &_mesh;
};
//=============================================================================

}// end namespace accelerator
}// end namespace calico

#endif


