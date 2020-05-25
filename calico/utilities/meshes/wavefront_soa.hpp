// Copyright (c) 2017-2019, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef CALICO_UTILITIES_MESHES_WAVEFRONT_SOA_HPP
#define CALICO_UTILITIES_MESHES_WAVEFRONT_SOA_HPP

#include <calico/math.hpp>

#include <vector>
#include <iostream>
#include <exception>
#include <sstream>
#include <string>

namespace calico {
namespace utilities {
namespace meshes {

class ParseError : public std::exception {
public:
    ParseError(const std::string &err) : _what(err) {}
    virtual const char* what() const noexcept {return _what.c_str();}
private:
    std::string _what;
};

/// Strip all leading and trailing spaces
template <typename String>
String trim(const String &input) {
    size_t start = input.find_first_not_of(" \t");
    size_t stop = input.find_last_not_of(" \t\n\r");
    if (stop - start) {
      return input.substr(start, stop-start+1);
    }
    return input;
}


/**
    This Wavefront OBJ loader follows the interface requirements of the Mesh
    adapter. 
    
    This simplistic loader is not general purpose as it does not handle
    anything with regards to objects or groups and ignores normals specified in
    the mesh file (both on the face and on vertices). Faces with more than
    three vertices are supported, but always constructed as a triangle fan (so
    must be convex). Faces with vertex normal indices are not supported.
    
    The whole mesh is loaded into a single object and exposed through the
    standard Mesh interface. It automatically computes surface normals and
    surface area of all facets loaded. It also splits Quad elements into a pair
    of Tri-elements. Finally, it internally uses 0 based vertex counting, so
    converts Wavefront OBJ's 1 based counting to 0 based during the load
    process.

    This class lays out its internal data structures using Structure of Arrays
    (SoA). For an Array of Structures (AoS) alternative, see WavefrontAoS.
*/
template <typename Float>
class WavefrontSoA {
public:
    typedef int FaceId;
    typedef std::size_t VertexId;
    typedef Float FloatType;
    
    static const FaceId ray_miss_id_c = -1;

    WavefrontSoA(std::istream& input_stream) : _min_x{ std::numeric_limits<decltype(_max_x)>::max() },
                                               _min_y{ std::numeric_limits<decltype(_max_x)>::max() },
                                               _min_z{ std::numeric_limits<decltype(_max_x)>::max() },
                                               _max_x{ std::numeric_limits<decltype(_max_x)>::lowest() },
                                               _max_y{ std::numeric_limits<decltype(_max_x)>::lowest() },
                                               _max_z{ std::numeric_limits<decltype(_max_x)>::lowest() }
    {

        std::size_t line_number(0u);
        while (input_stream.good()) {
            std::string line;
            std::getline(input_stream, line);
            line = trim(line);
            line_number += 1u;
            if (line.substr(0,2) == "v ") {
                // Parse the rest of the line as a vertex
                std::string vertex_coordinates = line.substr(1);
                std::size_t i(0u);
                std::size_t j(0u);

                // we're using signed int as an index from the face into the
                // vertex array. Make sure that we don't run out of indexing
                // space.
                if (_x.size() == std::numeric_limits<int>::max()) {
                    std::stringstream err;
                    err << "On line " << line_number << ", the maximum vertex limit was exceeded.";
                    throw ParseError(err.str());
                }

                try {
                    _x.push_back(std::stod(vertex_coordinates, &i));
                }
                catch (const std::exception &) {
                    std::stringstream err;
                    err << "Invalid X vertex (" << vertex_coordinates << ") on line " << line_number;
                    throw ParseError(err.str());
                }

                try {
                    _y.push_back(std::stod(vertex_coordinates.substr(i), &j));
                    i += j;
                }
                catch (const std::exception &) {
                    std::stringstream err;
                    err << "Invalid Y vertex (" << vertex_coordinates << ") on line " << line_number;
                    throw ParseError(err.str());
                }

                try {
                    _z.push_back(std::stod(vertex_coordinates.substr(i), &j));
                    i += j;
                }
                catch (const std::exception &) {
                    std::stringstream err;
                    err << "Invalid Z vertex (" << vertex_coordinates << ") on line " << line_number;
                    throw ParseError(err.str());
                }
            }
            else if (line.substr(0,2) == "f ") {
                std::string face_indices = line.substr(1);
                int a{},b{},c{};
                std::size_t offset{0}, j{0};
                try {
                    a = std::stoi(face_indices, &j);
                    offset += j;
                    b = std::stoi(face_indices.substr(offset), &j);
                    offset += j;
                    c = std::stoi(face_indices.substr(offset), &j);
                    offset += j;
                }
                catch (const std::exception &) {
                    std::stringstream err;
                    err << "Invalid face on line " << line_number << ". Faces "
                        << "must contain 3 or 4 vertex IDs";
                    throw ParseError(err.str());
                }

                if (a < 0) {
                    a = static_cast<int>(_x.size()) + a;
                } else {
                    a -= 1;
                }
                if (b < 0) {
                    b = static_cast<int>(_x.size()) + b;
                } else {
                    b -= 1;
                }
                if (c < 0) {
                    c = static_cast<int>(_x.size()) + c;
                } else {
                    c -= 1;
                }

				if (a < 0 || b < 0 || c < 0) {
					throw ParseError("Invalid vertex on line " + std::to_string(line_number));
				}

                _face[0].push_back(static_cast<std::size_t>(a));
                _face[1].push_back(static_cast<std::size_t>(b));
                _face[2].push_back(static_cast<std::size_t>(c));

                // I don't know if this is true in general, but I'm assuming
                // that polygons are always convex and can be treated as
                // triangle fans.
                while (true) {
                    b = c;

                    // Are there more vertices in this face?
                    try {
                        // We could read another vetex index, so yes, it's a fan.
                        c = std::stoi(face_indices.substr(offset), &j);
                        offset += j;
                    }
                    catch (const std::exception &) {
                        // That's ok. there must not be any more vertices, so
                        // just ignore this exception.
                        break;
                    }
                    if (c < 0) {
                        c = static_cast<int>(_x.size()) + c;
                    } else {
                        c -= 1;
                    }

                    if (a < 0 || b < 0 || c < 0) {
                        throw ParseError("Invalid vertex on line " + std::to_string(line_number));
                    }

                    _face[0].push_back(static_cast<std::size_t>(a));
                    _face[1].push_back(static_cast<std::size_t>(b));
                    _face[2].push_back(static_cast<std::size_t>(c));
                }
            }
            else if (line.substr(0, 1) == "#") {
                // Ignore comment lines
            }
            else if (line.empty()) {
                // Ignore blank lines
            }
            else {
                std::string::size_type a = line.find_first_of(" ");
                std::cerr << "Warning: Unrecognized line type: " << line.substr(0, a) << std::endl;
            }
        }

        _count = _face[0].size();
        _normal_x.resize(_count);
        _normal_y.resize(_count);
        _normal_z.resize(_count);
        _d.resize(_count);
        _area.resize(_count);

        for (std::size_t i{ 0 }; i < _count; ++i) {
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

    Float x(std::size_t index, VertexId corner) const {
        return _x.at(_face[corner].at(index));
    }
    Float y(std::size_t index, VertexId corner) const {return _y.at(_face[corner].at(index));}
    Float z(std::size_t index, VertexId corner) const {return _z.at(_face[corner].at(index));}
    
    Float normal_x(std::size_t index) const {return _normal_x.at(index);}
    Float normal_y(std::size_t index) const {return _normal_y.at(index);}
    Float normal_z(std::size_t index) const {return _normal_z.at(index);}

    Float d(std::size_t index) const {return _d.at(index);}

    Float area(std::size_t index) const {return _area.at(index);}

    std::size_t size() const {return _count;}

    /// Convert a face index into a face ID.
    FaceId index_to_face_id(std::size_t const index) const {return static_cast<FaceId>(index);}

private:
    std::vector<Float> _normal_x;
    std::vector<Float> _normal_y;
    std::vector<Float> _normal_z;

    std::vector<Float> _d;

    std::vector<Float> _x;
    std::vector<Float> _y;
    std::vector<Float> _z;

    std::vector<std::size_t> _face[3];

    std::vector<Float> _area;

    Float _min_x;
    Float _min_y;
    Float _min_z;

    Float _max_x;
    Float _max_y;
    Float _max_z;

    std::size_t _count;
};
//=============================================================================

}// end calico::utilities::meshes
}// end calico::utilities
}// end calico

#endif

