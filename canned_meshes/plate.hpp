// Copyright (c) 2014, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef __CALICO__CANNED_MESHES__PLATE__HPP__
#define __CALICO__CANNED_MESHES__PLATE__HPP__

namespace calico {
namespace canned_meshes {

/**
    This Plate follows the interface requirements of the Mesh adapter.  It
    defines a flat plate located on the XY plane with normals pointing along
    +Z.
*/
template <typename Float>
class Plate {
public:
    typedef int FaceId;
    typedef int VertexId;

    Plate() {
        _normal_x.resize(2);
        _normal_y.resize(2);
        _normal_z.resize(2);
        _x.resize(2*3);
        _y.resize(2*3);
        _z.resize(2*3);
        _area.resize(3);

        // Face 1
        _x[0] = 1;
        _y[0] = 1;
        _z[0] = 0;

        _x[1] = 0;
        _y[1] = 1;
        _z[1] = 0;

        _x[2] = 1;
        _y[2] = 0;
        _z[2] = 0;

        _normal_x[0] = 0;
        _normal_y[0] = 0;
        _normal_z[0] = 1;


        // Face 2
        _x[3] = 1;
        _y[3] = 0;
        _z[3] = 0;

        _x[4] = 0;
        _y[4] = 0;
        _z[4] = 0;

        _x[5] = 1;
        _y[5] = 0;
        _z[5] = 0;

        _normal_x[1] = 0;
        _normal_y[1] = 0;
        _normal_z[1] = 1;

        for (FaceId i = 0u; i < 3u; ++i) {
            Float xx, yy, zz;
            math::cross(x(i, 1) - x(i, 0),
                        y(i, 1) - y(i, 0),
                        z(i, 1) - z(i, 0),
                        x(i, 2) - x(i, 0),
                        y(i, 2) - y(i, 0),
                        z(i, 2) - z(i, 0),
                        xx, yy, zz);
            _area[i] = math::length(xx, yy, zz) * Float(0.5);

            _min_x = std::min(x(i, 0), _min_x);
            _min_y = std::min(y(i, 0), _min_y);
            _min_z = std::min(z(i, 0), _min_z);

            _min_x = std::min(x(i, 1), _min_x);
            _min_y = std::min(y(i, 1), _min_y);
            _min_z = std::min(z(i, 1), _min_z);

            _min_x = std::min(x(i, 2), _min_x);
            _min_y = std::min(y(i, 2), _min_y);
            _min_z = std::min(z(i, 2), _min_z);


            _max_x = std::max(x(i, 0), _max_x);
            _max_y = std::max(y(i, 0), _max_y);
            _max_z = std::max(z(i, 0), _max_z);

            _max_x = std::max(x(i, 1), _max_x);
            _max_y = std::max(y(i, 1), _max_y);
            _max_z = std::max(z(i, 1), _max_z);

            _max_x = std::max(x(i, 2), _max_x);
            _max_y = std::max(y(i, 2), _max_y);
            _max_z = std::max(z(i, 2), _max_z);
        }
    }

    void bounding_box(Float &min_x, Float &min_y, Float &min_z, 
                      Float &max_x, Float &max_y, Float &max_z) const
    {
        min_x = _min_x;
        min_y = _min_y;
        min_z = _min_z;

        max_x = _max_x;
        max_y = _max_y;
        max_z = _max_z;
    }
    
    FaceId ray_miss_id() const {return -1;}

    Float x(FaceId id, VertexId corner) const {return _x[id*3+corner];}
    Float y(FaceId id, VertexId corner) const {return _y[id*3+corner];}
    Float z(FaceId id, VertexId corner) const {return _z[id*3+corner];}
    
    Float normal_x(FaceId id) const {return _normal_x[id];}
    Float normal_y(FaceId id) const {return _normal_y[id];}
    Float normal_z(FaceId id) const {return _normal_z[id];}

    Float area(FaceId id) const {return _area[id];}

private:
    std::vector<Float> _normal_x;
    std::vector<Float> _normal_y;
    std::vector<Float> _normal_z;

    std::vector<Float> _x;
    std::vector<Float> _y;
    std::vector<Float> _z;

    std::vector<Float> _area;
};
//=============================================================================

}// end canned_meshes
}// end calico

#endif

