// Copyright (c) 2014, Nathan Smith <nathanjsmith@gmail.com>
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#ifndef __CALICO__TRACER__HPP__
#define __CALICO__TRACER__HPP__

#include <calico/math.hpp>

#include <iostream>

namespace calico {


template <typename Float, typename MeshAdapter, typename Accelerator>
class Tracer
{
public:
    Tracer(Accelerator &accelerator) 
        : _accelerator(accelerator)
    { }

    /// Trace rays against the mesh (via the accelerator) and find their
    /// intersection distance
    /**
        Rays are provided to the trace_rays routine in an array-like interface
        InputArray.  This interface is required to provide seven accessor methods:

          - InputArray::RayId size()                    // number of rays to trace
          - InputArray::RayId start_face(size_t index)  // face from which the ray 
                                                        // is being launched, or
                                                        // Mesh::ray_miss_id_c
          - Float             start_x(size_t index)     // start point of ray
          - Float             start_y(size_t index)
          - Float             start_z(size_t index)
          - Float             direction_x(size_t index) // (unit) direction of ray
          - Float             direction_y(size_t index)
          - Float             direction_z(size_t index)

        The Float routines must return a floating-point value of any precision
        so long as it is consistent. The RayId routines must return a value of
        InputArray::RayId.

        If InputArray::start_face(size_t index) returns a value other than
        MeshAdapter::ray_miss_id_c, then the face with the returned Id will
        *not* be considered for intersection tests. This allows the user to
        easily launch a ray from the surface of a face (e.g. in transmission
        or reflection calculations) without concern for whether the ray is
        numerically beyond the face. If all faces are hit candidates, then
        InputArray::start_face(size_t) should return MeshAdapter::ray_miss_id_c.

        The ResultArray must provide just a single method named set_result that
        has the signature:

            set_hit(MeshAdapter::FaceId face, Real t, 
                    Real hit_x, Real hit_y, Real hit_z);

        where face is the face struck by the ray (or MeshAdapter::ray_miss_id_c
        if it did not hit any face), t is the distance along the ray that the
        intersection occured (e.g. hit = start + direction * t), and
        hit_{x,y,z} are the x, y and z components of the intersection position.
    */
    template <typename InputArray, typename ResultArray>
    void trace_rays(const InputArray &input, ResultArray &result)
    {
        const typename InputArray::RayId ray_count = input.size();

        // Outerloop on rays so that we can use an accelerator structure inside
        Float t(0.), hit_x(0.), hit_y(0.), hit_z(0.);
        for (std::size_t i = 0u; i < ray_count; ++i) 
        {
            typename MeshAdapter::FaceId face(input.start_face(i));

            _accelerator.find_intersection(input.start_x(i),
                                           input.start_y(i),
                                           input.start_z(i),
                                           input.direction_x(i),
                                           input.direction_y(i),
                                           input.direction_z(i),
                                           face, t,
                                           hit_x, hit_y, hit_z);
            result.set_hit(i, face, t, hit_x, hit_y, hit_z);
        }// end loop over rays
    }

private:
    /// Adapter interface to accelerator being traced.  The accelerator may
    /// _not_ change during tracing.
    Accelerator   &_accelerator;
};
//=============================================================================


template <typename Float, typename MeshAdapter, typename Accelerator>
Tracer<Float, MeshAdapter, Accelerator>
  make_tracer(const MeshAdapter &adapter, Accelerator &accelerator) {
    return Tracer<Float, MeshAdapter, Accelerator>(accelerator);
  }


}// end namespace calico

#endif
