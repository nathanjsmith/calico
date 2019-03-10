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

#ifndef __CALICO__ACCELERATOR__SIMPLE_BVH__HPP__
#define __CALICO__ACCELERATOR__SIMPLE_BVH__HPP__

#include <calico/math.hpp>
#include <calico/accelerator/bounding_box.hpp>

namespace calico {
namespace accelerator {

/// Find the minimum component of the 3-tuple x1, x2, x3
template <typename Float, typename FloatMeta>
Float min_v(Float x1, Float x2, Float x3) {
    return FloatMeta::min(x1, FloatMeta::min(x2, x3));
}

/// Find the maximum component of the 3-tuple x1, x2, x3
template <typename Float, typename FloatMeta>
Float max_v(Float x1, Float x2, Float x3) {
    return FloatMeta::max(x1, FloatMeta::max(x2, x3));
}

/// Compute the centroid of element index from the Mesh mesh and store the xyz
/// 3-tuple into x, y and z
template <typename Float, typename Mesh>
void centroid(Mesh const &mesh, std::size_t index, Float &x, Float &y, Float &z) {
    constexpr Float third{Float(1.)/Float(3.)};
    x = (mesh.x(index, 0) + mesh.x(index, 1) + mesh.x(index, 2)) * third;
    y = (mesh.y(index, 0) + mesh.y(index, 1) + mesh.y(index, 2)) * third;
    z = (mesh.z(index, 0) + mesh.z(index, 1) + mesh.z(index, 2)) * third;
}


/**
    Given a mesh and a containment test, this will build a Bounding Volume
    Hierarchy (BVH) using axis-aligned bounding boxes sorted into a hierarchy.
    When a ray is provided to the accelerator, it will traverse the tree
    according to the boxes that encapsulate part of the tree to determine which
    triangles a ray *may* strike. Once the subset of triangles that are
    potential hit candidates have been determined, it will then test those
    triangles for intersection.
*/
template <typename Float, 
          typename Mesh, 
          typename Containment, 
          typename FloatMeta=calico::math::StdTypeInterface<Float> >
class SimpleBvh {
public:
    typedef Float       FloatType;
    typedef Mesh        MeshType;
    typedef Containment ContainmentType;

    SimpleBvh(Mesh const &mesh) : _mesh(mesh) {
        // Compute:
        // 1. bounding box and centroid for each triangle (tb and c respectively)
        // 2. Bounding box for _all_ the triangles (vb, or voxel-bounds)
        // 3. Bounding box for _all_ the centroids (cb)
        _triangle_min_x.resize(mesh.size(), FloatMeta::greatest());
        _triangle_min_y.resize(mesh.size(), FloatMeta::greatest());
        _triangle_min_z.resize(mesh.size(), FloatMeta::greatest());

        _triangle_max_x.resize(mesh.size(), FloatMeta::lowest());
        _triangle_max_y.resize(mesh.size(), FloatMeta::lowest());
        _triangle_max_z.resize(mesh.size(), FloatMeta::lowest());

        _centroid[0].resize(mesh.size(), Float(0.));
        _centroid[1].resize(mesh.size(), Float(0.));
        _centroid[2].resize(mesh.size(), Float(0.));

        std::size_t const mesh_size = mesh.size();
        for (std::size_t i{0u}; i < mesh_size; ++i) {
            // Compute tb_i
            _triangle_min_x[i] = min_v<Float, FloatMeta>(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _triangle_max_x[i] = max_v<Float, FloatMeta>(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _triangle_min_y[i] = min_v<Float, FloatMeta>(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _triangle_max_y[i] = max_v<Float, FloatMeta>(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _triangle_min_z[i] = min_v<Float, FloatMeta>(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            _triangle_max_z[i] = max_v<Float, FloatMeta>(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            centroid(_mesh, i, _centroid[0][i], _centroid[1][i], _centroid[2][i]);
        }

        // Allocate space for the maximum number of nodes (2N-1).
        _nodes.resize(mesh.size()*2+1);

        _triangle_indices.resize(mesh.size());
        std::iota(std::begin(_triangle_indices), std::end(_triangle_indices), std::size_t(0));

        // Node 0 holds the whole tree. So its bounding box is the bounding box
        // of the whole mesh
        mesh.bounding_box(_nodes[0u].min[0], _nodes[0u].min[1], _nodes[0u].min[2],
                          _nodes[0u].max[0], _nodes[0u].max[1], _nodes[0u].max[2]);
        _nodes[0u].leaf  = true;
        _nodes[0u].start = 0u;
        _nodes[0u].stop  = mesh.size();
        _next_available_node = 2u;

        recursive_subdivide(0u);
    }

    ~SimpleBvh() = default;
    SimpleBvh(SimpleBvh const &) = default;
    SimpleBvh(SimpleBvh &&) = default;

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
                            is found. Will be set to FloatMeta::greatest() if no 
                            intersection is found (output)
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
        Float min_t = 0., max_t = FloatMeta::greatest();
        typename Mesh::FaceId const ignore_face{face};
        face = Mesh::ray_miss_id_c;
        t = FloatMeta::greatest();

        // Pre-compute these to accelerate bounding box calculations
        Float inverse_direction_x{Float(1) / direction_x};
        Float inverse_direction_y{Float(1) / direction_y};
        Float inverse_direction_z{Float(1) / direction_z};
        recursive_find_intersection(start_x, start_y, start_z,
                          direction_x, direction_y, direction_z, 
                          inverse_direction_x, inverse_direction_y, inverse_direction_z, 
                          0u, ignore_face, face, min_t, max_t, hit_x, hit_y, hit_z);
        if (face != Mesh::ray_miss_id_c) {
            t = min_t;
        }
    }

private:


    struct Node {
        std::size_t start; ///< If leaf, index of first triangle in the node
                           ///< Else, index of left child; right child at start+1
        std::size_t stop;  ///< If leaf, terminating index of children (like
                           ///< std::vector::end(). supports loops of i =
                           ///< start; i < stop; ++i); Else, index of right child.
        bool leaf;         ///< True, Node is a leaf; else node in BVH
        Float min[3];      ///< AABB min of elements in the node
        Float max[3];      ///< AABB max of elements in the node

        std::size_t count() const {return stop - start;}
    };

    /**
      Given a leaf node, determine whether it should be subdivided. If so,
      convert the leaf into a node, create a left and right child, and divide
      the elements between the left and right child, reordering them in the
      _triangle_indices array in the process.

      @param node_index     Index of the leaf-node to consider for a split.
    */
    void recursive_subdivide(std::size_t node_index)
    {
        if (_nodes[node_index].count() <= std::size_t(4)) {
            // Don't subdivide.
            return;
        }

        // Make the next level by arbitrarily splitting the node in the middle
        // of its widest axis
        std::uint8_t widest_axis{0};
        Float const magnitudes[3] = {_nodes[node_index].max[0] - _nodes[node_index].min[0],
                                     _nodes[node_index].max[1] - _nodes[node_index].min[1],
                                     _nodes[node_index].max[2] - _nodes[node_index].min[2]};
        if (magnitudes[0] > magnitudes[1] && magnitudes[0] > magnitudes[2]) {
            widest_axis = std::uint8_t(0);
        }
        else if (magnitudes[1] > magnitudes[2]) {
            widest_axis = std::uint8_t(1);
        }
        else {
            widest_axis = std::uint8_t(2);
        }

        // Midpoint split point
        Float split_point = _nodes[node_index].min[widest_axis] + magnitudes[widest_axis] * 0.5;

        // Copy the triangle indices out of the set of elements we're going to
        // overwrite
        // TODO: I'm sure there's something really simple we can do to avoid
        //       making this copy, but for my first pass I'm going to just do
        //       this for the simplicity and ensure I get the expected outcome.
        std::vector<std::size_t> old_triangle_indices(_nodes[node_index].count());
        std::copy(_triangle_indices.begin() + _nodes[node_index].start, // start index
                  _triangle_indices.begin() + _nodes[node_index].count(), // stop index
                  old_triangle_indices.begin());

        // Initialize two new leaf nodes that we'll made children after
        // promoting the provided leaf-node into a tree node.
        std::size_t left_node = _next_available_node++;
        std::size_t right_node = _next_available_node++;
        _nodes[left_node].leaf = true;
        _nodes[left_node].min[0] = FloatMeta::greatest();
        _nodes[left_node].min[1] = FloatMeta::greatest();
        _nodes[left_node].min[2] = FloatMeta::greatest();
        _nodes[left_node].max[0] = FloatMeta::lowest();
        _nodes[left_node].max[1] = FloatMeta::lowest();
        _nodes[left_node].max[2] = FloatMeta::lowest();

        _nodes[right_node].leaf = true;
        _nodes[right_node].min[0] = FloatMeta::greatest();
        _nodes[right_node].min[1] = FloatMeta::greatest();
        _nodes[right_node].min[2] = FloatMeta::greatest();
        _nodes[right_node].max[0] = FloatMeta::lowest();
        _nodes[right_node].max[1] = FloatMeta::lowest();
        _nodes[right_node].max[2] = FloatMeta::lowest();

        // Sort the triangles to the left/right of the split point
        std::size_t const count = old_triangle_indices.size();
        std::size_t left_index{_nodes[node_index].start};
        std::size_t right_index{_nodes[node_index].count()-1};
        std::size_t left_count{0};
        std::size_t right_count{0};
        std::size_t stored_node{0};
        std::size_t triangle_index{0u};
        for (std::size_t i{0u}; i < count; ++i) {
            triangle_index = old_triangle_indices[i];
            if (_centroid[widest_axis][triangle_index] < split_point) {
                _triangle_indices[left_index++] = triangle_index;
                ++left_count;
                stored_node = left_node;
            }
            else {
                _triangle_indices[right_index--] = triangle_index;
                ++right_count;
                stored_node = right_node;
            }

            // Update the node's bounding box using the triangle's actual
            // bounds
            _nodes[stored_node].min[0] = FloatMeta::min(_nodes[stored_node].min[0],
                                                        _triangle_min_x[triangle_index]);
            _nodes[stored_node].min[1] = FloatMeta::min(_nodes[stored_node].min[1],
                                                        _triangle_min_y[triangle_index]);
            _nodes[stored_node].min[2] = FloatMeta::min(_nodes[stored_node].min[2],
                                                        _triangle_min_z[triangle_index]);

            _nodes[stored_node].max[0] = FloatMeta::max(_nodes[stored_node].max[0],
                                                        _triangle_max_x[triangle_index]);
            _nodes[stored_node].max[1] = FloatMeta::max(_nodes[stored_node].max[1],
                                                        _triangle_max_y[triangle_index]);
            _nodes[stored_node].max[2] = FloatMeta::max(_nodes[stored_node].max[2],
                                                        _triangle_max_z[triangle_index]);
        }

        if (left_count == 0u || right_count == 0u) {
            // There's no point in performing this split. All of the triangles
            // ended up on the same side. Just leave the current node a leaf.
            // Restore the _next_available_node to its previous position.
            _next_available_node -= std::size_t(2);
        }
        else {
            // Update the start/stop indices
            _nodes[left_node].start = _nodes[node_index].start;
            _nodes[left_node].stop  = _nodes[node_index].start + left_count;

            _nodes[right_node].start = _nodes[node_index].start + left_count;
            _nodes[right_node].stop  = _nodes[node_index].stop;

            // Update the node we're subdividing to indicate it isn't a leaf,
            // but a node
            _nodes[node_index].leaf = false;
            _nodes[node_index].start = left_node;
            _nodes[node_index].stop  = right_node;

            // Subdivide the left and right halves
            recursive_subdivide(left_index);
            recursive_subdivide(right_index);
        }
    }
    
    /// Internal search of nodes for intersections
    void recursive_find_intersection(
                           const Float start_x, 
                           const Float start_y, 
                           const Float start_z,
                           const Float direction_x, 
                           const Float direction_y, 
                           const Float direction_z,
                           const Float inverse_direction_x, 
                           const Float inverse_direction_y, 
                           const Float inverse_direction_z,
                           const std::size_t node_id,
                           typename Mesh::FaceId const ignore_face,
                           typename Mesh::FaceId &face, 
                           Float &min_t, Float &max_t,
                           Float &hit_x, Float &hit_y, Float &hit_z)
    {
        // Should we explore this node?
        Float local_min_t = 0., local_max_t = FloatMeta::greatest();
        if (!calico::accelerator::intersects<Float, FloatMeta>(
                start_x, start_y, start_z,
                direction_x, direction_y, direction_z,
                inverse_direction_x, inverse_direction_y, inverse_direction_z,
                _nodes[node_id].min[0], _nodes[node_id].min[1], _nodes[node_id].min[2], 
                _nodes[node_id].max[0], _nodes[node_id].max[1], _nodes[node_id].max[2], 
                local_min_t, local_max_t))
        {
            return;
        }
        if (local_min_t > min_t) {
            // It's a candidate, but the closest hit possible is further away
            // than the closest hit found thus far. Ignore this box.
            return;
        }

        // Yes, this node is a potential hit candidate.  Is it a leaf node?
        if (_nodes[node_id].leaf) {
            // Leaf node. Search the node's triangles for an intersection
            for (std::size_t i{_nodes[node_id].start}; i < _nodes[node_id].stop; ++i) {
                std::size_t const face_index{_triangle_indices[i]};
                typename Mesh::FaceId const f{_mesh.index_to_face_id(face_index)};
                if (f == ignore_face) {
                    continue;
                }

                // If the ray lies in the plane, this should return Inf as the
                // intersection distance.  That will get filtered out in the
                // following tests.
                Float tmp_t = math::ray_plane_intersection(_mesh.normal_x(face_index),
                                                           _mesh.normal_y(face_index),
                                                           _mesh.normal_z(face_index),
                                                           _mesh.d(face_index),
                                                           start_x,
                                                           start_y,
                                                           start_z,
                                                           direction_x,
                                                           direction_y,
                                                           direction_z);

                // Only keep hits in front of the start point.  Only calculate the
                // intersection point if the hit distance is better than our
                // previous best result.
                if (Float(0.) < tmp_t && tmp_t <= min_t) {
                    // This is closer than the closest we've found yet. Make
                    // sure the intersection is inside of the triangle.
                    Float tmp_hit_x = start_x + (direction_x * tmp_t);
                    Float tmp_hit_y = start_y + (direction_y * tmp_t);
                    Float tmp_hit_z = start_z + (direction_z * tmp_t);

                    if (Containment::contains(_mesh, face_index, tmp_hit_x, tmp_hit_y, tmp_hit_z))
                    {
                        hit_x = tmp_hit_x;
                        hit_y = tmp_hit_y;
                        hit_z = tmp_hit_z;
                        face  = f;
                        min_t = tmp_t;
                    }
                }
            }// end loop over triangles in this node
        }// end if node is leaf
        else {
            // This node is the parent of two more bounding boxes.
            recursive_find_intersection(start_x, start_y, start_z,
                                        direction_x, direction_y, direction_z,
                                        inverse_direction_x, inverse_direction_y, inverse_direction_z,
                                        _nodes[node_id].start, ignore_face, face, 
                                        local_min_t, local_max_t,
                                        hit_x, hit_y, hit_z);

            recursive_find_intersection(start_x, start_y, start_z,
                                        direction_x, direction_y, direction_z,
                                        inverse_direction_x, inverse_direction_y, inverse_direction_z,
                                        _nodes[node_id].stop, ignore_face, face, 
                                        local_min_t, local_max_t,
                                        hit_x, hit_y, hit_z);
        }
    }

    /// Reference to the Mesh we're accelerating search
    Mesh const &_mesh;

    /// All of the BVH nodes (and leafs) in the BVH tree
    std::vector<Node> _nodes;

    // TODO: Would the bounding boxes be more efficiently accessed as a
    //       contiguous array of 6 values? I.e. is AoS better than SoA for
    //       bounding box storage/access? I think it might be so that it
    //       occupies one cache line, not six.

    std::vector<Float> _triangle_min_x; ///< Component of per-triangle bounds (tb)
    std::vector<Float> _triangle_min_y; ///< Component of per-triangle bounds (tb)
    std::vector<Float> _triangle_min_z; ///< Component of per-triangle bounds (tb)
    std::vector<Float> _triangle_max_x; ///< Component of per-triangle bounds (tb)
    std::vector<Float> _triangle_max_y; ///< Component of per-triangle bounds (tb)
    std::vector<Float> _triangle_max_z; ///< Component of per-triangle bounds (tb)

    std::vector<Float> _centroid[3];    // Triangle centroid (c)

    /// The nodes point to these indices. The indices point to triangles in the Mesh
    std::vector<std::size_t> _triangle_indices; 

    std::size_t _next_available_node; ///< Next available node index
};

}// end namespace calico::accelerator
}// end namespace calico

#endif
