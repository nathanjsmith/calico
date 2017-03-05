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

namespace calico {
namespace accelerator {

/**
  interface to the routines in std:: that we need for running the routines in
  accelerator. It provides:

    max_infinity   -- A routine that returns +infinity
    min_infinity   -- A routine that returns -infinity
    max            -- A routine that returns the larger of a and b
    min            -- A routine that returns the smaller of a and b
*/
template <typename Float>
struct StdTypeInterface {
  inline constexpr static Float min_infinity() {
    return -std::numeric_limits<Float>::infinity();
  }
  
  inline constexpr static Float max_infinity() {
    return std::numeric_limits<Float>::infinity();
  }
  
  inline static Float min(Float a, Float b) {return std::min(a,b);}
  
  inline static Float max(Float a, Float b) {return std::max(a,b);}
};
//=============================================================================


/**
    Find the intersection of a ray and an axis-aligned bounding box. This
    returns the part of the ray that lies within the bounding box in the
    near/far values.
*/
template <typename Float, typename interface=StdTypeInterface<Float>>
bool intersects(const Float ray_start_x,     const Float ray_start_y,     const Float ray_start_z,
                const Float ray_direction_x, const Float ray_direction_y, const Float ray_direction_z,
                const Float min_x,           const Float min_y,           const Float min_z, 
                const Float max_x,           const Float max_y,           const Float max_z)
{
    // This implementation was informed by the post at 
    // https://tavianator.com/fast-branchless-raybounding-box-intersections/,
    // but is an independent implementation of that algorithm.
    Float min_t = interface::min_infinity();
    Float max_t = interface::max_infinity();

    if (ray_direction_x != Float(0)) {
        const Float mn_x = (min_x - ray_start_x) / ray_direction_x;
        const Float mx_x = (max_x - ray_start_x) / ray_direction_x;

        min_t = interface::max(min_t, interface::min(mn_x, mx_x));
        max_t = interface::min(max_t, interface::max(mn_x, mx_x));
    }

    if (ray_direction_y != Float(0)) {
        const Float mn_y = (min_y - ray_start_y) / ray_direction_y;
        const Float mx_y = (max_y - ray_start_y) / ray_direction_y;

        min_t = interface::max(min_t, interface::min(mn_y, mx_y));
        max_t = interface::min(max_t, interface::max(mn_y, mx_y));
    }

    if (ray_direction_z != Float(0)) {
        const Float mn_z = (min_z - ray_start_z) / ray_direction_z;
        const Float mx_z = (max_z - ray_start_z) / ray_direction_z;

        min_t = interface::max(min_t, interface::min(mn_z, mx_z));
        max_t = interface::min(max_t, interface::max(mn_z, mx_z));
    }

    return (max_t > min_t && max_t > Float(0.)) || max_t == Float(0.);
}


/**
    Find the intersection of a ray and an axis-aligned bounding box. This
    returns the part of the ray that lies within the bounding box in the
    near/far values.
*/
template <typename Float, typename interface=StdTypeInterface<Float>>
bool unsafe_intersects(
        const Float ray_start_x,             const Float ray_start_y,             const Float ray_start_z,
        const Float inverse_ray_direction_x, const Float inverse_ray_direction_y, const Float inverse_ray_direction_z,
        const Float min_x,                   const Float min_y,                   const Float min_z, 
        const Float max_x,                   const Float max_y,                   const Float max_z)
{
    // This implementation was informed by the post at 
    // https://tavianator.com/fast-branchless-raybounding-box-intersections/,
    // but is an independent implementation of that algorithm.
    Float min_t = interface::min_infinity();
    Float max_t = interface::max_infinity();

    const Float mn_x = (min_x - ray_start_x) * inverse_ray_direction_x;
    const Float mx_x = (max_x - ray_start_x) * inverse_ray_direction_x;

    min_t = interface::max(min_t, interface::min(mn_x, mx_x));
    max_t = interface::min(max_t, interface::max(mn_x, mx_x));

    const Float mn_y = (min_y - ray_start_y) * inverse_ray_direction_y;
    const Float mx_y = (max_y - ray_start_y) * inverse_ray_direction_y;

    min_t = interface::max(min_t, interface::min(mn_y, mx_y));
    max_t = interface::min(max_t, interface::max(mn_y, mx_y));

    const Float mn_z = (min_z - ray_start_z) * inverse_ray_direction_z;
    const Float mx_z = (max_z - ray_start_z) * inverse_ray_direction_z;

    min_t = interface::max(min_t, interface::min(mn_z, mx_z));
    max_t = interface::min(max_t, interface::max(mn_z, mx_z));

    return max_t > interface::max(min_t, 0.0) || max_t == Float(0.);
}

template <typename Float>
class BoundingBoxes {
public:
    BoundingBoxes() : min_x(nullptr), min_y(nullptr), min_z(nullptr),
                      max_x(nullptr), max_y(nullptr), max_z(nullptr) {}

    Float *min_x;
    Float *min_y;
    Float *min_z;
    Float *max_x;
    Float *max_y;
    Float *max_z;
};


}// end namespace calico::accelerator
}// end namespace calico

#endif
