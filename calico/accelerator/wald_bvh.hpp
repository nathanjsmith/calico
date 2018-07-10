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

namespace bvh_wald_util {
    template <typename T>
    T* build_array(std::size_t count, T fill_value) {
        T* a{new T[count]};
        std::fill(a, a+count, fill_value);
        return a;
    }
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
        // From section 3.2 of the paper, we need:
        // 1. bounding box and centroid for each triangle (tb and c respectively)
        // 2. Bounding box for _all_ the triangles (vb, or voxel-bounds)
        // 3. Bounding box for _all_ the centroids (cb)
        _triangle_min_x = bvh_wald_util::build_array<Float>(mesh.size(), limits::max());
        _triangle_min_y = bvh_wald_util::build_array<Float>(mesh.size(), limits::max());
        _triangle_min_z = bvh_wald_util::build_array<Float>(mesh.size(), limits::max());

        _triangle_max_x = bvh_wald_util::build_array<Float>(mesh.size(), limits::lowest());
        _triangle_max_y = bvh_wald_util::build_array<Float>(mesh.size(), limits::lowest());
        _triangle_max_z = bvh_wald_util::build_array<Float>(mesh.size(), limits::lowest());

        _centroid[0] = bvh_wald_util::build_array<Float>(mesh.size(), 0.);
        _centroid[1] = bvh_wald_util::build_array<Float>(mesh.size(), 0.);
        _centroid[2] = bvh_wald_util::build_array<Float>(mesh.size(), 0.);

        _centroid_min_x = limits::max();
        _centroid_min_y = limits::max();
        _centroid_min_z = limits::max();
        _centroid_max_x = limits::lowest();
        _centroid_max_y = limits::lowest();
        _centroid_max_z = limits::lowest();

        for (std::size_t i{0u}; i < mesh.size(); ++i) {
            // Compute tb_i
            _triangle_min_x[i] = min_v(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _triangle_max_x[i] = max_v(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _triangle_min_y[i] = min_v(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _triangle_max_y[i] = max_v(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _triangle_min_z[i] = min_v(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            _triangle_max_z[i] = max_v(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            centroid(_mesh, i, _centroid[0][i], _centroid[1][i], _centroid[2][i]);

            // Compute triangle's contribution to vb
            _min_x = std::min(_triangle_min_x[i], _min_x);
            _min_y = std::min(_triangle_min_y[i], _min_y);
            _min_z = std::min(_triangle_min_z[i], _min_z);

            _max_x = std::max(_triangle_max_x[i], _max_x);
            _max_y = std::max(_triangle_max_y[i], _max_y);
            _max_z = std::max(_triangle_max_z[i], _max_z);

            // Compute triangle's contribution to cb
            _centroid_min_x = std::min(_triangle_min_x[i], _centroid_min_x);
            _centroid_min_y = std::min(_triangle_min_y[i], _centroid_min_y);
            _centroid_min_z = std::min(_triangle_min_z[i], _centroid_min_z);

            _centroid_max_x = std::max(_triangle_max_x[i], _centroid_max_x);
            _centroid_max_y = std::max(_triangle_max_y[i], _centroid_max_y);
            _centroid_max_z = std::max(_triangle_max_z[i], _centroid_max_z);
        }

        // We need a continuous array of storage representing the indices of
        // all triangles.
        triangle_indices = new std::size_t[mesh.size()];
        std::iota(triangle_indices, triangle_indices+mesh.size(), std::size_t(0u));

        // allocate space for the maximum number of nodes (2N-1).

        // First, bounding boxes for each node, initialized to be empty (min > max at extremes)
        _node_box_min_x   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::max());
        _node_box_min_y   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::max());
        _node_box_min_z   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::max());
        _node_box_max_x   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::lowest());
        _node_box_max_y   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::lowest());
        _node_box_max_z   = bvh_wald_util::build_array<Float>(mesh.size()*2u, limits::lowest());

        // Next, the child nodes. If the child node is 0, then this is not a node, but a leaf.
        _node_child_left  = bvh_wald_util::build_array<std::size_t>(mesh.size()*2u, std::size_t(0));
        _node_child_right = bvh_wald_util::build_array<std::size_t>(mesh.size()*2u, std::size_t(0));

        // Finally, how many triangles are in this leaf? The element _start_triangle_index[n] 
        // has meaning only if _node_child_left[n] != 0
        _start_triangle_index = bvh_wald_util::build_array<std::size_t>(mesh.size()*2u, std::size_t(0));
        _triangle_count = bvh_wald_util::build_array<std::size_t>(mesh.size()*2u, std::size_t(0));

        // TODO: partition triangles across bins
    }

    ~WaldBvh() {
        delete[] _triangle_min_x;
        delete[] _triangle_min_y;
        delete[] _triangle_min_z;

        delete[] _triangle_max_x;
        delete[] _triangle_max_y;
        delete[] _triangle_max_z;

        delete[] _centroid[0];
        delete[] _centroid[1];
        delete[] _centroid[2];
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

    Float *_triangle_min_x; // O(N)
    Float *_triangle_min_y; // O(N)
    Float *_triangle_min_z; // O(N)
    Float *_triangle_max_x; // O(N)
    Float *_triangle_max_y; // O(N)
    Float *_triangle_max_z; // O(N)

    Float *_centroid[3];    // O(N)

    Float _min_x;
    Float _min_y;
    Float _min_z;
    Float _max_x;
    Float _max_y;
    Float _max_z;

    Float _centroid_min_x;
    Float _centroid_min_y;
    Float _centroid_min_z;
    Float _centroid_max_x;
    Float _centroid_max_y;
    Float _centroid_max_z;

    std::size_t *triangle_indices;

    Float *_node_box_min_x;
    Float *_node_box_min_y;
    Float *_node_box_min_z;
    Float *_node_box_max_x;
    Float *_node_box_max_y;
    Float *_node_box_max_z;

    std::size_t *_node_child_left;
    std::size_t *_node_child_right;

    std::size_t *_start_triangle_index;
    std::size_t *_triangle_count;
};

}// end namespace calico::accelerator
}// end namespace calico

#endif
