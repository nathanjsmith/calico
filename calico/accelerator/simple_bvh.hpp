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

#include <iterator>
#include <algorithm>
#include <cassert>

namespace calico {
namespace accelerator {

template <typename Float>
bool close(Float const a, Float const b) {
  Float const delta = std::abs(a - b);
  Float const largest = std::max(std::abs(a), std::abs(b));
  return (delta / largest) < 1e-6;
}

template <typename Float>
void assert_close(Float const a, Float const b) {
  if (!close(a,b)) {
    std::cerr << a << " is not close to " << b << std::endl;
    assert(false);
  }
}

class RaiiScopeMarker
{
public:
  RaiiScopeMarker(std::ostream &stream, std::string start, std::string stop) : _stream(stream), _stop(stop) {
    _stream << start;
  }
  ~RaiiScopeMarker() {
    _stream << _stop;
  }
private:
  std::ostream &_stream;
  std::string   _stop;
};


template <typename Node>
void export_bvh_to_obj(std::vector<Node> const &nodes, std::size_t index, std::size_t parent, std::ostream &stream)
{
    if (nodes[index].leaf) {
        stream << "g Leaf " << index << ", child of " << parent << "\n";
    }
    else {
        stream << "g Branch " << index << ", child of " << parent << "\n";
    }
    stream << "v " << nodes[index].min[0] << " " << nodes[index].min[1] << " " << nodes[index].min[2] << "\n"; // -8 (-1,-1,-1)
    stream << "v " << nodes[index].min[0] << " " << nodes[index].min[1] << " " << nodes[index].max[2] << "\n"; // -7 (-1,-1, 1)
    stream << "v " << nodes[index].min[0] << " " << nodes[index].max[1] << " " << nodes[index].min[2] << "\n"; // -6 (-1, 1,-1)
    stream << "v " << nodes[index].min[0] << " " << nodes[index].max[1] << " " << nodes[index].max[2] << "\n"; // -5 (-1, 1, 1)
    stream << "v " << nodes[index].max[0] << " " << nodes[index].min[1] << " " << nodes[index].min[2] << "\n"; // -4 ( 1,-1,-1)
    stream << "v " << nodes[index].max[0] << " " << nodes[index].min[1] << " " << nodes[index].max[2] << "\n"; // -3 ( 1,-1, 1)
    stream << "v " << nodes[index].max[0] << " " << nodes[index].max[1] << " " << nodes[index].min[2] << "\n"; // -2 ( 1, 1,-1)
    stream << "v " << nodes[index].max[0] << " " << nodes[index].max[1] << " " << nodes[index].max[2] << "\n"; // -1 ( 1, 1, 1)
    stream << "f -8 -7 -5 -6\n"; // left
    stream << "f -4 -2 -1 -3\n"; // right
    stream << "f -8 -4 -3 -7\n"; // front
    stream << "f -6 -5 -1 -2\n"; // back
    stream << "f -7 -3 -1 -5\n"; // top
    stream << "f -8 -6 -2 -4\n"; // bottom

    if (nodes[index].leaf) {
        return;
    }

    // Recurse
    export_bvh_to_obj(nodes, nodes[index].start, index, stream); // left
    export_bvh_to_obj(nodes, nodes[index].stop, index, stream);  // right
}

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
        _centroid[0].resize(mesh.size(), Float(0.));
        _centroid[1].resize(mesh.size(), Float(0.));
        _centroid[2].resize(mesh.size(), Float(0.));
        _triangle_bounds[0].resize(mesh.size(), Limits{FloatMeta::greatest(), FloatMeta::lowest()});
        _triangle_bounds[1].resize(mesh.size(), Limits{FloatMeta::greatest(), FloatMeta::lowest()});
        _triangle_bounds[2].resize(mesh.size(), Limits{FloatMeta::greatest(), FloatMeta::lowest()});

        std::cout << "Constructing bounding box around " << mesh.size() << " triangles." << std::endl;

        std::size_t const mesh_size = mesh.size();
        for (std::size_t i{0u}; i < mesh_size; ++i) {
            // Compute tb_i
            _triangle_bounds[0][i].min = min_v<Float, FloatMeta>(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));
            _triangle_bounds[0][i].max = max_v<Float, FloatMeta>(_mesh.x(i, 0), _mesh.x(i, 1), _mesh.x(i, 2));

            _triangle_bounds[1][i].min = min_v<Float, FloatMeta>(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));
            _triangle_bounds[1][i].max = max_v<Float, FloatMeta>(_mesh.y(i, 0), _mesh.y(i, 1), _mesh.y(i, 2));

            _triangle_bounds[2][i].min = min_v<Float, FloatMeta>(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));
            _triangle_bounds[2][i].max = max_v<Float, FloatMeta>(_mesh.z(i, 0), _mesh.z(i, 1), _mesh.z(i, 2));

            centroid(_mesh, i, _centroid[0][i], _centroid[1][i], _centroid[2][i]);
        }

        // Allocate space for the maximum number of nodes (2N-1).
        _nodes.resize(mesh.size()*2+1);

        _triangle_indices.resize(mesh.size());
        std::iota(std::begin(_triangle_indices), std::end(_triangle_indices), std::size_t(0));

        // Node 0 holds the whole tree. So its bounding box is the bounding box
        // of the whole mesh
        std::cout << "Before querying bounding box, node 0's bounds are:\n";
        std::cout << "(" << _nodes[0].min[0] << ", " << _nodes[0].min[1] << ", " << _nodes[0].min[2] << ") - ("
                         << _nodes[0].max[0] << ", " << _nodes[0].max[1] << ", " << _nodes[0].max[2] << ")" 
                         << std::endl;

        mesh.bounding_box(_nodes[0u].min[0], _nodes[0u].min[1], _nodes[0u].min[2],
                          _nodes[0u].max[0], _nodes[0u].max[1], _nodes[0u].max[2]);
        _nodes[0u].leaf  = true;
        _nodes[0u].start = 0u;
        _nodes[0u].stop  = mesh.size();
        _next_available_node = 1u;

        std::cout << "Root node is now:\n(" << _nodes[0].min[0] << ", " << _nodes[0].min[1] << ", " << _nodes[0].min[2] << ") - ("
                                      << _nodes[0].max[0] << ", " << _nodes[0].max[1] << ", " << _nodes[0].max[2] << ")" 
                                      << std::endl;

        recursive_subdivide(/*node_index=*/0u, /*level=*/1u, /*max_levels=*/3000u);

        std::cout << "Traversing computed BVH and writing out boxes at each level:" << std::endl;
        export_bvh_to_obj(_nodes, 0, 0, std::cout);
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
                           Float &hit_x, Float &hit_y, Float &hit_z) const
    {
        Float min_t = FloatMeta::greatest(), max_t = FloatMeta::greatest();
        typename Mesh::FaceId const ignore_face{face};
        face = Mesh::ray_miss_id_c;
        t = FloatMeta::greatest();

        //// std::cout << "find_intersection(start=(" << start_x << ", "
        ////                                          << start_y << ", "
        ////                                          << start_z << "),\n"
        ////              "                  direction=(" << direction_x << ", "
        ////                                              << direction_y << ", "
        ////                                              << direction_z << "),\n"
        ////              "                  face=" << face << ",\n"
        ////              "                  t=" << t << ",\n"
        ////              "                  hit=(" << hit_x << ", " << hit_y << ", " << hit_z << ")\n";

        // Pre-compute these to accelerate bounding box calculations
        Float inverse_direction_x{Float(1) / direction_x};
        Float inverse_direction_y{Float(1) / direction_y};
        Float inverse_direction_z{Float(1) / direction_z};
        recursive_find_intersection(start_x, start_y, start_z,
                          direction_x, direction_y, direction_z, 
                          inverse_direction_x, inverse_direction_y, inverse_direction_z, 
                          0u, ignore_face, face, min_t, max_t, hit_x, hit_y, hit_z);
        //// std::cout << "After recursive_find_intersection:" << std::endl;
        //// std::cout << "  min_t = " << min_t << std::endl;
        //// std::cout << "  max_t = " << max_t << std::endl;
        //// std::cout << "  face  = " << face << std::endl;
        //// std::cout << "  hit_x = " << hit_x << std::endl;
        //// std::cout << "  hit_y = " << hit_y << std::endl;
        //// std::cout << "  hit_z = " << hit_z << std::endl;

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

    struct Limits {
        Float min;
        Float max;
    };

    /**
      Given a leaf node, determine whether it should be subdivided. If so,
      convert the leaf into a node, create a left and right child, and divide
      the elements between the left and right child, reordering them in the
      _triangle_indices array in the process.

      @param node_index     Index of the leaf-node to consider for a split.
    */
    void recursive_subdivide(std::size_t node_index, std::size_t level, std::size_t max_levels)
    {
        //// std::cout << "Considering split of node " << node_index << " {" << std::endl;
        if (level == max_levels) {
            //// std::cout << "  Hit max levels. Stopping subdivision" << std::endl;
            //// std::cout << "}" << std::endl;
            // Don't subdivide.
            return;
        }
        if (_nodes[node_index].count() <= std::size_t(4)) {
            //// std::cout << "  Node is too small to subdivide. It has only " << _nodes[node_index].count() << " elements." << std::endl;
            //// std::cout << "}" << std::endl;
            // Don't subdivide.
            return;
        }

        // Perform a debugging/sanity check
        {
          Float min[3] = {FloatMeta::greatest(), FloatMeta::greatest(), FloatMeta::greatest()};
          Float max[3] = {FloatMeta::lowest(), FloatMeta::lowest(), FloatMeta::lowest()};
          std::size_t d_index{_nodes[node_index].start};
          std::size_t const d_right_index{d_index + _nodes[node_index].count()-1};
          while (d_index <= d_right_index) {
            std::size_t i = _triangle_indices[d_index++];
            min[0] = std::min(_triangle_bounds[0][i].min, min[0]);
            min[1] = std::min(_triangle_bounds[1][i].min, min[1]);
            min[2] = std::min(_triangle_bounds[2][i].min, min[2]);

            max[0] = std::max(_triangle_bounds[0][i].max, max[0]);
            max[1] = std::max(_triangle_bounds[1][i].max, max[1]);
            max[2] = std::max(_triangle_bounds[2][i].max, max[2]);
          }

          // std::cerr << "Current box holds triangle_indices [" << _nodes[node_index].start << ", " << _nodes[node_index].count

          assert_close(min[0], _nodes[node_index].min[0]);
          assert_close(min[1], _nodes[node_index].min[1]);
          assert_close(min[2], _nodes[node_index].min[2]);
          assert_close(max[0], _nodes[node_index].max[0]);
          assert_close(max[1], _nodes[node_index].max[1]);
          assert_close(max[2], _nodes[node_index].max[2]);
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

        // Wald BVH says we should create bins around the bounding box of the
        // widest axis. Use a K of 16
        constexpr std::size_t SORT_BINS(16);
        constexpr Float INV_SORT_BINS(Float(1.) / SORT_BINS);
        Float const bin_width         = magnitudes[widest_axis] / Float(SORT_BINS);
        Float const inverse_bin_width = Float(1) / bin_width;
        Float const axis_min          = _nodes[node_index].min[widest_axis];
        std::vector<std::size_t> bin_count(SORT_BINS, std::size_t(0));

        // Sort triangles into bins according to their midpoints
        const std::size_t triangle_count(_nodes[node_index].count());
        Float const k_1 = (SORT_BINS * (Float(1.) - 1e-6)) / magnitudes[widest_axis];
        for (std::size_t local_triangle_index{0}; local_triangle_index < triangle_count; ++local_triangle_index) {
            std::size_t const triangle_index = _triangle_indices[_nodes[node_index].start + local_triangle_index];
            std::size_t const bin_index = k_1 * (_centroid[widest_axis][triangle_index] - axis_min);

            if (bin_index >= SORT_BINS) {
              std::cout << "Expected a value [0, " << SORT_BINS-1 << "], but computed " << bin_index << std::endl;
              std::cout << "  BVH Level:         " << level << std::endl;
              std::cout << "  Axis min:          " << axis_min << std::endl;
              std::cout << "  Triangle centroid: " << _centroid[widest_axis][triangle_index] << std::endl;
              std::cout << "  Triangle ID:       " << triangle_index << std::endl;
              std::cout << "  Triangle AABB:     " << _triangle_bounds[widest_axis][triangle_index].min << ", " << _triangle_bounds[widest_axis][triangle_index].max << std::endl;
              std::cout << "  Node AABB:         " << _nodes[node_index].min[widest_axis] << ", " << _nodes[node_index].max[widest_axis] << std::endl;
              std::cout << "  Widest axis:       " << static_cast<int>(widest_axis) << std::endl;
              std::cout << "  Magnitudes:        " << magnitudes[0] << ", " << magnitudes[1] << ", " << magnitudes[2] << std::endl;
              std::cout << "  k_1:               " << k_1 << std::endl;
            }
            assert(bin_index < SORT_BINS);
            bin_count[bin_index] += 1;
        }

        // Iterate across bins to see where the least expensive split point is
        // located.
        std::size_t champ_bin{0};
        std::size_t left_triangles{0};
        Float       champ_cost{FloatMeta::greatest()};
        for (std::size_t bin_index{1}; bin_index < SORT_BINS; ++bin_index) {
            left_triangles += bin_count[bin_index];
            Float sah_cost = (INV_SORT_BINS * bin_index)               * left_triangles +                   // Cost for left-half
                             (INV_SORT_BINS * (SORT_BINS - bin_index)) * (triangle_count - left_triangles); // Cost for right-half
            if (sah_cost < champ_cost) {
                champ_bin = bin_index;
                champ_cost = sah_cost;
            }
        }

        // SAH split point
        Float const split_point = (champ_bin + 1) * bin_width + axis_min;
        //// std::cout << "  Found split point along axis-" << int(widest_axis) 
        ////           << " (" << _nodes[node_index].min[widest_axis] << " to " << _nodes[node_index].max[widest_axis] 
        ////           << ") at " << split_point << " in bin(" << champ_bin << "/" << SORT_BINS << ")" << std::endl;

        // Copy the triangle indices out of the set of elements we're going to
        // overwrite
        // TODO: I'm sure there's something really simple we can do to avoid
        //       making this copy, but for my first pass I'm going to just do
        //       this for the simplicity and ensure I get the expected outcome.
        std::vector<std::size_t> old_triangle_indices;// (_nodes[node_index].count());
        old_triangle_indices.reserve(_nodes[node_index].count());
        std::copy(_triangle_indices.begin() + _nodes[node_index].start, // start index
                  _triangle_indices.begin() + _nodes[node_index].start + _nodes[node_index].count(), // stop index
                  std::back_inserter(old_triangle_indices));
        // std::copy(_triangle_indices.begin() + _nodes[node_index].start, // start index
        //           _triangle_indices.begin() + _nodes[node_index].start + _nodes[node_index].count(), // stop index
        //           old_triangle_indices.begin());
        //// std::cout << "  Sought to copy " << _nodes[node_index].count() << " triangle indices" << std::endl;
        //// std::cout << "    old_triangle_indices holds " << old_triangle_indices.size() << " indices" << std::endl;
        assert(old_triangle_indices.size() == _nodes[node_index].count());

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
        std::size_t right_index{left_index + _nodes[node_index].count()-1};
        std::size_t left_count{0};
        std::size_t right_count{0};
        std::size_t stored_node{0};
        std::size_t triangle_index{0u};
        for (std::size_t i{0u}; i < count; ++i) {
            triangle_index = old_triangle_indices[i];
            if (_centroid[widest_axis][triangle_index] < split_point) {
                //// std::cout << "  Triangle " << triangle_index << " centroid at " << _centroid[widest_axis][triangle_index] << " < " << split_point << ", so sorted to the left at index (" << left_index << ")" << std::endl;
                _triangle_indices[left_index++] = triangle_index;
                ++left_count;
                stored_node = left_node;
            }
            else {
                //// std::cout << "  Triangle " << triangle_index << " centroid at " << _centroid[widest_axis][triangle_index] << " >= " << split_point << ", so sorted to the right at index (" << right_index << ")" << std::endl;
                _triangle_indices[right_index--] = triangle_index;
                ++right_count;
                stored_node = right_node;
            }

            // Update the node's bounding box using the triangle's actual
            // bounds
            // TODO: It may be more efficient to use a second loop and compute
            //       the x-dimension bounds, then the y-dimension bounds, then
            //       the z-dimension bounds. It is likely more cache (stream)
            //       friendly.
            _nodes[stored_node].min[0] = FloatMeta::min(_nodes[stored_node].min[0],
                                                        _triangle_bounds[0][triangle_index].min);
            _nodes[stored_node].min[1] = FloatMeta::min(_nodes[stored_node].min[1],
                                                        _triangle_bounds[1][triangle_index].min);
            _nodes[stored_node].min[2] = FloatMeta::min(_nodes[stored_node].min[2],
                                                        _triangle_bounds[2][triangle_index].min);

            _nodes[stored_node].max[0] = FloatMeta::max(_nodes[stored_node].max[0],
                                                        _triangle_bounds[0][triangle_index].max);
            _nodes[stored_node].max[1] = FloatMeta::max(_nodes[stored_node].max[1],
                                                        _triangle_bounds[1][triangle_index].max);
            _nodes[stored_node].max[2] = FloatMeta::max(_nodes[stored_node].max[2],
                                                        _triangle_bounds[2][triangle_index].max);
        }

        if (left_count == 0u || right_count == 0u) {
            // There's no point in performing this split. All of the triangles
            // ended up on the same side. Just leave the current node a leaf.
            // Restore the _next_available_node to its previous position
            // (subtract 2 to get rid of the left _and_ right nodes).
            _next_available_node -= std::size_t(2);
            std::cout << "  All triangles sorted to one side (left=" << left_count << ", right=" << right_count << "). Reverting to leaf." << std::endl;
        }
        else {
            //// std::cout << "  Success! We divided the triangles into (left=" << left_count << ", right=" << right_count << ") bins!" << std::endl;
            //// std::cout << "    The bounding box looks like:" << std::endl;
            //// std::cout << "    Left:  (" << _nodes[left_node].min[0] << ", " << _nodes[left_node].min[1] << ", " << _nodes[left_node].min[2] << ")" << std::endl;
            //// std::cout << "           (" << _nodes[left_node].max[0] << ", " << _nodes[left_node].max[1] << ", " << _nodes[left_node].max[2] << ")" << std::endl;
            //// std::cout << "    Right: (" << _nodes[right_node].min[0] << ", " << _nodes[right_node].min[1] << ", " << _nodes[right_node].min[2] << ")" << std::endl;
            //// std::cout << "           (" << _nodes[right_node].max[0] << ", " << _nodes[right_node].max[1] << ", " << _nodes[right_node].max[2] << ")" << std::endl;

            // Update the start/stop indices
            _nodes[left_node].start = _nodes[node_index].start;
            _nodes[left_node].stop  = _nodes[node_index].start + left_count;

            _nodes[right_node].start = _nodes[node_index].start + left_count;
            _nodes[right_node].stop  = _nodes[node_index].stop;

            //// std::cout << "  Left box:" << std::endl;
            //// std::cout << "    Start: " << _nodes[left_node].start << std::endl;
            //// std::cout << "    Stop:  " << _nodes[left_node].stop << std::endl;
            //// std::cout << "    Triangles:" << std::endl;
            //// for (std::size_t idx{_nodes[left_node].start}; idx < _nodes[left_node].stop; ++idx) {
            ////   std::cout << "      (" << idx << ") " << _triangle_indices[idx] << "\n";
            //// }
            //// std::cout << "  Right box:" << std::endl;
            //// std::cout << "    Start: " << _nodes[right_node].start << std::endl;
            //// std::cout << "    Stop:  " << _nodes[right_node].stop << std::endl;
            //// for (std::size_t idx{_nodes[right_node].start}; idx < _nodes[right_node].stop; ++idx) {
            ////   std::cout << "      (" << idx << ") " << _triangle_indices[idx] << "\n";
            //// }
            //// std::cout << "---------------" << std::endl;

            // Update the node we're subdividing to indicate it isn't a leaf,
            // but a node
            _nodes[node_index].leaf = false;
            _nodes[node_index].start = left_node;
            _nodes[node_index].stop  = right_node;

            // Subdivide the left and right halves
            recursive_subdivide(left_node, level+1, max_levels);
            recursive_subdivide(right_node, level+1, max_levels);
        }
        //// std::cout << "}" << std::endl;
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
                           Float &hit_x, Float &hit_y, Float &hit_z) const
    {
        //// RaiiScopeMarker scope(std::cout, "{\n", "}\n");
        //// std::cout << "recursive_find_intersection(start=(" << start_x << ", "
        ////                                                    << start_y << ", "
        ////                                                    << start_z << "),\n"
        ////              "                            direction=(" << direction_x << ", "
        ////                                                        << direction_y << ", "
        ////                                                        << direction_z << "),\n"
        ////              "                            node_id=" << node_id << ",\n"
        ////              "                            ignore_face=" << ignore_face << ",\n"
        ////              "                            face=" << face << ",\n"
        ////              "                            min_t=" << min_t << ",\n"
        ////              "                            max_t=" << max_t << ",\n"
        ////              "                            hit=(" << hit_x << ", " << hit_y << ", " << hit_z << ")\n";

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
            //// std::cout << "The ray doesn't strike this node's bounding box" << std::endl;
            return;
        }
        if (local_min_t > min_t) {
            //// std::cout << "Ray strikes the node's bounding box, but min_t is " << local_min_t << " > " << min_t << std::endl;
            // It's a candidate, but the closest hit possible is further away
            // than the closest hit found thus far. Ignore this box.
            return;
        }
        //// std::cout << "local_min_t is within bounds" << std::endl;

        // Yes, this node is a potential hit candidate.  Is it a leaf node?
        if (_nodes[node_id].leaf) {
            // Leaf node. Search the node's triangles for an intersection
            for (std::size_t i{_nodes[node_id].start}; i < _nodes[node_id].stop; ++i) {
                std::size_t const face_index{_triangle_indices[i]};
                typename Mesh::FaceId const f{_mesh.index_to_face_id(face_index)};
                // std::cout << "Considering triangle " << f << std::endl;
                if (f == ignore_face) {
                    // std::cout << "  Ignoring this triangle as instructed" << std::endl;
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
                //// std::cout << "Intersection with plane at " << tmp_t << std::endl;

                // Only keep hits in front of the start point.  Only calculate the
                // intersection point if the hit distance is better than our
                // previous best result.
                if (Float(0.) < tmp_t && tmp_t <= min_t) {
                    //// std::cout << "Hit at " << tmp_t << " appears valid 0 < " << tmp_t << " <= " << min_t << std::endl;

                    // This is closer than the closest we've found yet. Make
                    // sure the intersection is inside of the triangle.
                    Float tmp_hit_x = start_x + (direction_x * tmp_t);
                    Float tmp_hit_y = start_y + (direction_y * tmp_t);
                    Float tmp_hit_z = start_z + (direction_z * tmp_t);

                    if (Containment::contains(_mesh, face_index, tmp_hit_x, tmp_hit_y, tmp_hit_z))
                    {
                        //// std::cout << "Valid intersection!" << std::endl;
                        //// std::cout << "  hit (" << hit_x << ", " << hit_y << ", " << hit_z << ")\n";
                        //// std::cout << "  face " << face << "\n";
                        //// std::cout << "  min_t " << min_t << "\n";

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
            //// std::cout << "Node is not a leaf. Searching children" << std::endl;
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

    std::vector<Limits> _triangle_bounds[3];  ///< Min/max of the triangle along each dimension (tb)
    std::vector<Float>  _centroid[3];         ///< Triangle centroid (c)

    /// The nodes point to these indices. The indices point to triangles in the Mesh
    std::vector<std::size_t> _triangle_indices; 

    std::size_t _next_available_node; ///< Next available node index
};

}// end namespace calico::accelerator
}// end namespace calico

#endif
