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

#ifndef __CALICO__SOA_RESULT__HPP__
#define __CALICO__SOA_RESULT__HPP__

namespace calico {
namespace result {

template <typename Real, typename FaceId,
          typename RealArray, typename FaceIdArray>
class SoaResult
{
public:
    SoaResult(FaceIdArray &face_array, RealArray &t_array, 
              RealArray &x_hit_array, RealArray &y_hit_array, 
              RealArray &z_hit_array) :
        _face_array(face_array), _t_array(t_array), 
        _x_hit_array(x_hit_array), _y_hit_array(y_hit_array),
        _z_hit_array(z_hit_array)
    {
    }

    void set_hit(size_t index, const FaceId &face_id, 
                 const Real &t, const Real &x_hit, 
                 const Real &y_hit, const Real &z_hit)
    {
      _face_array[index] = face_id;
      _t_array[index] = t;
      _x_hit_array[index] = x_hit;
      _y_hit_array[index] = y_hit;
      _z_hit_array[index] = z_hit;
    }

private:
    FaceIdArray &_face_array;
    RealArray   &_t_array;
    RealArray   &_x_hit_array;
    RealArray   &_y_hit_array;
    RealArray   &_z_hit_array;
};
//=============================================================================


template <typename Real, typename FaceId, 
          typename RealArray, typename FaceIdArray>
SoaResult<Real, FaceId, RealArray, FaceIdArray>
  make_soa_result(FaceIdArray &face_array,
                  RealArray &t_array, 
                  RealArray &x_hit_array, 
                  RealArray &y_hit_array, 
                  RealArray &z_hit_array) {
  return SoaResult<Real, FaceId, RealArray, FaceIdArray>(face_array, 
                                                         t_array, 
                                                         x_hit_array, 
                                                         y_hit_array, 
                                                         z_hit_array);
}
//=============================================================================


}// end namespace calico::result
}// end namespace calico

#endif
