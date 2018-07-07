// Copyright (c) 2018, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef __CALICO__ACCELERATOR__WALD_BVH__HPP__
#define __CALICO__ACCELERATOR__WALD_BVH__HPP__

#include <calico/math.hpp>

#include <limits>

namespace calico {
namespace accelerator {

/*
template <typename Float,
          typename Mesh,
          std::size_t bin_count=16,
          typename limits=std::numeric_limits<Float> >
std::size_t compute_bin_id(
*/

template <typename Float>
Float min_v(Float x1, Float x2, Float x3) {
    return std::min(x1, std::min(x2, x3));
}

template <typename Float>
Float max_v(Float x1, Float x2, Float x3) {
    return std::max(x1, std::max(x2, x3));
}

template <typename Float, typename Mesh>
void centroid(Mesh const &mesh, std::size_t index, Float &x, Float &y, Float &z) {
    constexpr Float third{Float(1.)/Float(3.)};
    x = (mesh.x(index, 0) + mesh.x(index, 1) + mesh.x(index, 2)) * third;
    y = (mesh.y(index, 0) + mesh.y(index, 1) + mesh.y(index, 2)) * third;
    z = (mesh.z(index, 0) + mesh.z(index, 1) + mesh.z(index, 2)) * third;
}


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
class WaldBvh {
public:
    typedef Float       FloatType;
    typedef Mesh        MeshType;
    typedef Containment ContainmentType;

    WaldBvh(const Mesh &mesh) : _mesh(mesh) {
        // Compute bounding boxes for each face of the mesh
        _min_x = new Float[mesh.size()];
        _min_y = new Float[mesh.size()];
        _min_z = new Float[mesh.size()];

        _max_x = new Float[mesh.size()];
        _max_y = new Float[mesh.size()];
        _max_z = new Float[mesh.size()];

        _centroid_x = new Float[mesh.size()];
        _centroid_y = new Float[mesh.size()];
        _centroid_z = new Float[mesh.size()];

        for (std::size_t i{0u}; i < mesh.size(); ++i) {
            _min_x[i] = min_v(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _max_x[i] = max_v(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _min_y[i] = min_v(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _max_y[i] = max_v(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _min_z[i] = min_v(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            _max_z[i] = max_v(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            centroid(_mesh, i, _centroid_x[i], _centroid_y[i], _centroid_z[i]);
        }

        // TODO: build a hierarchy of bounding boxes
    }

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
    	typename Mesh::FaceId const ignore_face{face};
    }

private:
    const Mesh &_mesh;

    Float *_min_x;
    Float *_min_y;
    Float *_min_z;
    Float *_max_x;
    Float *_max_y;
    Float *_max_z;

    Float *_centroid_x;
    Float *_centroid_y;
    Float *_centroid_z;
};

}// end namespace calico::accelerator
}// end namespace calico

#endif
