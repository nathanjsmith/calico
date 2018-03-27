// Copyright (c) 2017, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef __CALICO__ACCELERATOR__BOUNDING_BOX__HPP__
#define __CALICO__ACCELERATOR__BOUNDING_BOX__HPP__

#include <calico/math.hpp>

#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <memory>

#include <stdlib.h> // for aligned_alloc

namespace calico {
namespace accelerator {
template <typename Float>
constexpr Float ize_scalar();

template <>
constexpr float ize_scalar<float>() {return 1.00000024f;}

template <>
constexpr double ize_scalar<double>() {return 1.0000000000000004;}

/**
    MaxMult bounding box intersection algorithm (mostly) as described in:

    Ize, Thiago. "Robust BVH ray traversal." Journal of Computer Graphics
                 Techniques (JCGT) 2.2 (2013): 12-27.

    A key difference is that this does not require that the bounding box values
    be stored in an array indexable via 0/1 (min/max). Instead, it uses ternary
    expressions to choose between the min/max bounds. It should be nearly as
    performant.

    Also, this uses the template function ize_scalar, which is specialized on
    the floating-point type selected. If you use a type other than float or
    double, you must also provide a specialization of

        template<>
        YourFloat calico::accelerator::ize_scalar<YourFloat>() {return 1.;}

    for your type.  An arbitrary precision type can have the specialization
    just return 1. See the source paper for details on choosing a MaxMult
    scalar.

    @param ray_start_x      Start x value for the ray (ray origin.x)
    @param ray_start_y      Start y value for the ray (ray origin.y)
    @param ray_start_z      Start z value for the ray (ray origin.z)

    @param ray_direction_x  Direction x value for the ray ; direction must be unit length
    @param ray_direction_y  Direction y value for the ray ; direction must be unit length
    @param ray_direction_z  Direction z value for the ray ; direction must be unit length

    @param min_x            Bounding box minimum X value
    @param min_y            Bounding box minimum Y value
    @param min_z            Bounding box minimum Z value
    @param max_x            Bounding box maximum X value
    @param max_y            Bounding box maximum Y value
    @param max_z            Bounding box maximum Z value

    Inputs/Outputs:

    @param min_t            Minimum distance along the vector that should be considered
    @param max_t            Maximum distance along the vector that should be considered

    @return                 true if the ray intersects the bounding box, false
                            if it does not.
*/
template <typename Float, typename interface=calico::math::StdTypeInterface<Float>>
bool intersects(
        const Float ray_start_x,             const Float ray_start_y,             const Float ray_start_z,
        const Float ray_direction_x,         const Float ray_direction_y,         const Float ray_direction_z,
        const Float inverse_ray_direction_x, const Float inverse_ray_direction_y, const Float inverse_ray_direction_z,
        const Float min_x,                   const Float min_y,                   const Float min_z, 
        const Float max_x,                   const Float max_y,                   const Float max_z,
        Float &min_t, Float &max_t)
{
    min_t = 0.;

    // This implementation requires that min(a,b) return b if a is NaN. We can
    // accomplish that using a ternary operator, but std::min or some other min
    // may not behave that way. See the Ize paper for details. In IEEE 754,
    // comparing a number with NaN will yield false (unless the comparison is
    // !=, which none of the below are), so the right-hand value will always be
    // returned in these min/max calculations.
    struct nan_safe {
        static inline Float min(const Float &left, const Float &right) {
            return (left < right ? left : right);
        }
    
        static inline float max(const float &left, const float &right) {
            return (left > right ? left : right);
        }
    };

    const Float min_x_t = ((ray_direction_x < Float(0) ? max_x : min_x) - ray_start_x) * inverse_ray_direction_x;
    const Float max_x_t = ((ray_direction_x < Float(0) ? min_x : max_x) - ray_start_x) * inverse_ray_direction_x;

    const Float min_y_t = ((ray_direction_y < Float(0) ? max_y : min_y) - ray_start_y) * inverse_ray_direction_y;
    const Float max_y_t = ((ray_direction_y < Float(0) ? min_y : max_y) - ray_start_y) * inverse_ray_direction_y;

    const Float min_z_t = ((ray_direction_z < Float(0) ? max_z : min_z) - ray_start_z) * inverse_ray_direction_z;
    const Float max_z_t = ((ray_direction_z < Float(0) ? min_z : max_z) - ray_start_z) * inverse_ray_direction_z;

    min_t = nan_safe::max(min_z_t, nan_safe::max(min_y_t, nan_safe::max(min_x_t, min_t)));
    max_t = nan_safe::min(max_z_t, nan_safe::min(max_y_t, nan_safe::min(max_x_t, max_t))) * ize_scalar<Float>();

    return min_t <= max_t;
}


/**
    For a pointer that was allocated using malloc, calloc, aligned_alloc, etc.,
    free that memory and assign nullptr. The pointer is passed by reference so
    that it can be assigned a value.
*/
template <typename T>
void release(T &p) {
    free(p);
    p = nullptr;
}

/**
    Collection of bounding boxes and mesh-face centroids stored using SoA
    layout.  Memory is automatically allocated on construction and released on
    destruction. See compute_mesh_bounds_and_centroids() for constructing a new
    bounding-box collection from a Mesh.
*/
template <typename Float>
class BoundingBoxes {
public:
    BoundingBoxes() : min_x(nullptr), min_y(nullptr), min_z(nullptr),
                      max_x(nullptr), max_y(nullptr), max_z(nullptr),
                      centroid_x(nullptr), centroid_y(nullptr), centroid_z(nullptr),
                      count(0u) {}

    BoundingBoxes(const std::size_t count_) 
                    : min_x(nullptr), min_y(nullptr), min_z(nullptr),
                      max_x(nullptr), max_y(nullptr), max_z(nullptr),
                      centroid_x(nullptr), centroid_y(nullptr), centroid_z(nullptr),
                      count(count_) 
    {
        // Use the aligned_alloc function from C11 (/not/ C++11, just C11)
        min_x = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        min_y = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        min_z = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));

        max_x = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        max_y = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        max_z = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));

        centroid_x = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        centroid_y = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
        centroid_z = static_cast<Float*>(aligned_alloc(16, sizeof(Float)*count));
    }

    ~BoundingBoxes() {
        release(min_x);
        release(min_y);
        release(min_z);
        release(max_x);
        release(max_y);
        release(max_z);
        release(centroid_x);
        release(centroid_y);
        release(centroid_z);
        count = 0u;
    }

    std::size_t size() const {return count;}

    Float      *min_x;
    Float      *min_y;
    Float      *min_z;
    Float      *max_x;
    Float      *max_y;
    Float      *max_z;
    Float      *centroid_x;
    Float      *centroid_y;
    Float      *centroid_z;
    std::size_t count;
};


/**
    Compute the bounding boxes of all triangles in the mesh. Also compute the
    centroid of each mesh face.

    @return Unique pointer to a new BoundingBoxes object
*/
template <typename MeshAdapter, typename FloatInterface=calico::math::StdTypeInterface<typename MeshAdapter:: FloatType>>
std::unique_ptr<BoundingBoxes<typename MeshAdapter::FloatType>>
    compute_mesh_bounds_and_centroids(const MeshAdapter &mesh) 
{
    constexpr typename MeshAdapter::FloatType one_third{1./3.};
    auto bounds = 
        std::unique_ptr<BoundingBoxes<typename MeshAdapter::FloatType>>(
                new BoundingBoxes<typename MeshAdapter::FloatType>(mesh.size()));
    const typename MeshAdapter::FaceIdIterator last{mesh.end_face_id()};
    for (typename MeshAdapter::FaceIdIterator i{mesh.begin_face_id()}; i != last; ++i) {
        bounds->min_x[i] = FloatInterface::min(mesh.x(i, 0), FloatInterface::min(mesh.x(i, 1), mesh.x(i, 2)));
        bounds->min_y[i] = FloatInterface::min(mesh.y(i, 0), FloatInterface::min(mesh.y(i, 1), mesh.y(i, 2)));
        bounds->min_z[i] = FloatInterface::min(mesh.z(i, 0), FloatInterface::min(mesh.z(i, 1), mesh.z(i, 2)));

        bounds->max_x[i] = FloatInterface::max(mesh.x(i, 0), FloatInterface::max(mesh.x(i, 1), mesh.x(i, 2)));
        bounds->max_y[i] = FloatInterface::max(mesh.y(i, 0), FloatInterface::max(mesh.y(i, 1), mesh.y(i, 2)));
        bounds->max_z[i] = FloatInterface::max(mesh.z(i, 0), FloatInterface::max(mesh.z(i, 1), mesh.z(i, 2)));

        bounds->centroid_x[i] = (mesh.x(i, 0) + mesh.x(i, 1) + mesh.x(i, 2)) * one_third;
        bounds->centroid_y[i] = (mesh.y(i, 0) + mesh.y(i, 1) + mesh.y(i, 2)) * one_third;
        bounds->centroid_z[i] = (mesh.z(i, 0) + mesh.z(i, 1) + mesh.z(i, 2)) * one_third;
    }

    return std::move(bounds);
}


}// end namespace calico::accelerator
}// end namespace calico

#endif
