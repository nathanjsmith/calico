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
  
  inline static Float min(Float a, Float b) {std::min(a, b);}
  
  inline static Float max(Float a, Float b) {std::max(a, b);}
};
//=============================================================================


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
template <typename Float, typename interface=StdTypeInterface<Float>>
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


template <typename Float>
class BoundingBoxes {
public:
    BoundingBoxes() : min_x(nullptr), min_y(nullptr), min_z(nullptr),
                      max_x(nullptr), max_y(nullptr), max_z(nullptr),
                      count(0u) {}

    Float      *min_x;
    Float      *min_y;
    Float      *min_z;
    Float      *max_x;
    Float      *max_y;
    Float      *max_z;
    std::size_t count;
};


}// end namespace calico::accelerator
}// end namespace calico

#endif
