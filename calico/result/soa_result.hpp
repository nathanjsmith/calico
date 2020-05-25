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

#ifndef CALICO_SOA_RESULT_HPP
#define CALICO_SOA_RESULT_HPP

namespace calico {
namespace result {

template <typename Real, typename FaceId,
          typename RealContainer, typename FaceIdContainer>
class SoaResult
{
public:
    SoaResult(FaceIdContainer &face_array, RealContainer &t_array, 
              RealContainer &x_hit_array, RealContainer &y_hit_array, 
              RealContainer &z_hit_array) :
        _face_array(face_array), _t_array(t_array), 
        _x_hit_array(x_hit_array), _y_hit_array(y_hit_array),
        _z_hit_array(z_hit_array)
    {
    }

    void set_hit(std::size_t index, const FaceId &face_id, 
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
    FaceIdContainer &_face_array;
    RealContainer   &_t_array;
    RealContainer   &_x_hit_array;
    RealContainer   &_y_hit_array;
    RealContainer   &_z_hit_array;
};
//=============================================================================


template <typename Real, typename FaceId, 
          typename RealContainer, typename FaceIdContainer>
SoaResult<Real, FaceId, RealContainer, FaceIdContainer>
  make_soa_result(FaceIdContainer &face_array,
                  RealContainer &t_array, 
                  RealContainer &x_hit_array, 
                  RealContainer &y_hit_array, 
                  RealContainer &z_hit_array) {
  return SoaResult<Real, FaceId, RealContainer, FaceIdContainer>(face_array, 
                                                         t_array, 
                                                         x_hit_array, 
                                                         y_hit_array, 
                                                         z_hit_array);
}
//=============================================================================


}// end namespace calico::result
}// end namespace calico

#endif
