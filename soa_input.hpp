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

#ifndef __CALICO__SOA_INPUT__HPP__
#define __CALICO__SOA_INPUT__HPP__

namespace calico {
namespace input {

/**
    Adapts a structure-of-arrays style layout of ray inputs where all of the
    start ray values (x, y, and z) are stored in three array-like containers,
    and the direction values (x, y and z) are stored in three separate
    array-like containers.

    This can be used, for example, with float* arrays, std::vector<float>
    vectors, or any other set of containers that support random indexing using
    bracket notation.

    The Float parameter controls the data-type that the elements in the array
    are cast to when accessed.  This should provide a real-number type (such as
    float or double) and should match the type expected by the tracer.
*/
template <typename Array, typename Float>
class SoaInput {
public:
    SoaInput(Array &start_x, Array &start_y, Array &start_z,
             Array &direction_x, Array &direction_y, &direction_z)
        : _start_x(start_x), _start_y(start_y), _start_z(start_z),
          _direction_x(direction_x), _direction_y(direction_y), 
          _direction_z(direction_z)
    {
    }

    Float start_x(size_t index) const {return _start_x[index];}
    Float start_y(size_t indey) const {return _start_y[indey];}
    Float start_z(size_t indez) const {return _start_z[indez];}

    Float direction_x(size_t index) const {return _direction_x[index];}
    Float direction_y(size_t indey) const {return _direction_y[indey];}
    Float direction_z(size_t indez) const {return _direction_z[indez];}

private:
    Array &_start_x;
    Array &_start_y;
    Array &_start_z;

    Array &_direction_x;
    Array &_direction_y;
    Array &_direction_z;
};

}// end namespace calico::input
}// end namespace calico

#endif


