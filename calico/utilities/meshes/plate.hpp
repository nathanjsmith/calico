// Copyright (c) 2014-2018, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef CALICO_UTILITIES_MESHES_PLATE_HPP
#define CALICO_UTILITIES_MESHES_PLATE_HPP

#include <calico/math.hpp>

#include <vector>
#include <iostream>

namespace calico {
namespace utilities {
namespace meshes {

/**
    This Plate follows the interface requirements of the Mesh adapter.  It
    defines a flat plate located on the XY plane with normals pointing along
    +Z.
*/
template <typename Float, typename interface=calico::math::StdTypeInterface<Float>>
class Plate {
public:
    typedef int FaceId;
    typedef std::size_t VertexId;
    typedef Float FloatType;
    
    static const FaceId ray_miss_id_c = -1;

    Plate() {
        _normal_x.resize(2);
        _normal_y.resize(2);
        _normal_z.resize(2);
        _x.resize(2*3);
        _y.resize(2*3);
        _z.resize(2*3);
        _area.resize(2);
        _d.resize(2);

        Float plate_z = 3.;

        // Face 1
        _x[0] = 1;
        _y[0] = -1;
        _z[0] = plate_z;

        _x[1] = 1;
        _y[1] = 1;
        _z[1] = plate_z;

        _x[2] = -1;
        _y[2] = 1;
        _z[2] = plate_z;


        // Face 2
        _x[3] = -1;
        _y[3] = 1;
        _z[3] = plate_z;

        _x[4] = -1;
        _y[4] = -1;
        _z[4] = plate_z;

        _x[5] = 1;
        _y[5] = -1;
        _z[5] = plate_z;

        _min_x = interface::max_infinity();
        _min_y = interface::max_infinity();
        _min_z = interface::max_infinity();

        _max_x = interface::min_infinity();
        _max_y = interface::min_infinity();
        _max_z = interface::min_infinity();

        for (std::size_t i{ 0 }; i < 2u; ++i) {
            // Automatically compute the normal and d for the triangle
            math::cross(x(i,2) - x(i,1), y(i,2) - y(i,1), z(i,2) - z(i,1),
                        x(i,0) - x(i,1), y(i,0) - y(i,1), z(i,0) - z(i,1),
                        _normal_x[i], _normal_y[i], _normal_z[i]);
            math::normalize(_normal_x[i], _normal_y[i], _normal_z[i]);

            _d[i] = math::dot(x(i,0), y(i,0), z(i,0),
                             _normal_x[i], _normal_y[i], _normal_z[i]);

            _area[i] = math::area(x(i,0), y(i,0), z(i,0),
                                  x(i,1), y(i,1), z(i,1),
                                  x(i,2), y(i,2), z(i,2));

            // Now update the bounding box
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

    FaceId index_to_face_id(std::size_t const index) const {return static_cast<FaceId>(index);}

    Float x(std::size_t index, VertexId corner) const {return _x.at(index*3+corner);}
    Float y(std::size_t index, VertexId corner) const {return _y.at(index*3+corner);}
    Float z(std::size_t index, VertexId corner) const {return _z.at(index*3+corner);}
    
    Float normal_x(std::size_t index) const {return _normal_x.at(index);}
    Float normal_y(std::size_t index) const {return _normal_y.at(index);}
    Float normal_z(std::size_t index) const {return _normal_z.at(index);}

    Float d(std::size_t index) const {return _d.at(index);}

    Float area(std::size_t index) const {return _area.at(index);}

    std::size_t size() const {return _normal_x.size();}

private:
    std::vector<Float> _normal_x;
    std::vector<Float> _normal_y;
    std::vector<Float> _normal_z;

    std::vector<Float> _d;

    std::vector<Float> _x;
    std::vector<Float> _y;
    std::vector<Float> _z;

    std::vector<Float> _area;

    Float _min_x;
    Float _min_y;
    Float _min_z;

    Float _max_x;
    Float _max_y;
    Float _max_z;
};
//=============================================================================

}// end calico::utilities::meshes
}// end calico::utilities
}// end calico

#endif

