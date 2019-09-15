// Copyright (c) 2017-2018, Nathan Smith <nathanjsmith@gmail.com>
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

#ifndef __CALICO__UTILITIES__MESHES__WAVEFRONT_AOS__HPP__
#define __CALICO__UTILITIES__MESHES__WAVEFRONT_AOS__HPP__

#include <calico/math.hpp>
#include <calico/utilities/meshes/wavefront_soa.hpp>

#include <vector>
#include <iostream>
#include <exception>
#include <sstream>
#include <string>
#include <cstdint>

namespace calico {
namespace utilities {
namespace meshes {

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

    This class lays out its internal data structures using Array of Structures
    (AoS). For a Structure of Arrays (SoA) alternative, see WavefrontSoA.
*/
template <typename Float, typename interface=calico::math::StdTypeInterface<Float>>
class WavefrontAoS {
public:
    typedef std::int32_t FaceId;
    typedef std::size_t VertexId;
    typedef Float FloatType;
    
    static const FaceId ray_miss_id_c = -1;

    WavefrontAoS(std::istream &input_stream) {

        _min_x = interface::max_infinity();
        _min_y = interface::max_infinity();
        _min_z = interface::max_infinity();

        _max_x = interface::min_infinity();
        _max_y = interface::min_infinity();
        _max_z = interface::min_infinity();

        std::size_t line_number(0u);
        while (input_stream.good()) {
            std::string line;
            std::getline(input_stream, line);
            line = trim(line);
            line_number += 1u;
            if (line.substr(0,2) == "v ") {
                Float x{},y{},z{};
                // Parse the rest of the line as a vertex
                std::string vertex_coordinates = line.substr(1);
                std::size_t i(0u);
                std::size_t j(0u);

                try {
                    x = std::stod(vertex_coordinates, &i);
                }
                catch (const std::exception &e) {
                    std::stringstream err;
                    err << "Invalid X vertex on line " << line_number;
                    throw ParseError(err.str());
                }

                try {
                    y = std::stod(vertex_coordinates.substr(i), &j);
                    i += j;
                }
                catch (const std::exception &e) {
                    std::stringstream err;
                    err << "Invalid Y vertex on line " << line_number;
                    throw ParseError(err.str());
                }

                try {
                    z = std::stod(vertex_coordinates.substr(i), &j);
                    i += j;
                }
                catch (const std::exception &e) {
                    std::stringstream err;
                    err << "Invalid Z vertex on line " << line_number;
                    throw ParseError(err.str());
                }

                _vertices.push_back({x,y,z});
            }
            else if (line.substr(0,2) == "f ") {
                std::string face_indices = line.substr(1);
                std::int32_t a{},b{},c{},d{};
                std::size_t offset{0}, j{0};
                try {
                    a = std::stoi(face_indices, &j);
                    offset += j;
                    b = std::stoi(face_indices.substr(offset), &j);
                    offset += j;
                    c = std::stoi(face_indices.substr(offset), &j);
                    offset += j;
                }
                catch (const std::exception &e) {
                    std::stringstream err;
                    err << "Invalid face on line " << line_number << ". Faces "
                        << "must contain 3 or 4 vertex IDs";
                    throw ParseError(err.str());
                }

                if (a < 0) {
                    a = _vertices.size() + a;
                } else {
                    a -= 1;
                }
                if (b < 0) {
                    b = _vertices.size() + b;
                } else {
                    b -= 1;
                }
                if (c < 0) {
                    c = _vertices.size() + c;
                } else {
                    c -= 1;
                }

                _faces.push_back({a,b,c});

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
                    catch (const std::exception &e) {
                        // That's ok. This isn't a quad, so just ignore this
                        // exception.
                        break;
                    }
                    if (c < 0) {
                        c = _vertices.size() + c;
                    } else {
                        c -= 1;
                    }

                    _faces.push_back({a,b,c});
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

        _count = _faces.size();
        _planes.resize(_count);
        _areas.resize(_count);

        for (FaceId i = 0u; i < _count; ++i) {
            // Automatically compute the normal and d for the triangle
            math::cross(x(i,2) - x(i,1), y(i,2) - y(i,1), z(i,2) - z(i,1),
                        x(i,0) - x(i,1), y(i,0) - y(i,1), z(i,0) - z(i,1),
                        _planes[i].x, _planes[i].y, _planes[i].z);
            math::normalize(_planes[i].x, _planes[i].y, _planes[i].z);

            _planes[i].d = math::dot(x(i,0), y(i,0), z(i,0),
                                     _planes[i].x, _planes[i].y, _planes[i].z);

            _areas[i] = math::area(x(i,0), y(i,0), z(i,0),
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

    Float x(std::size_t index, VertexId corner) const {return _vertices.at(_faces.at(index)[corner]).x;}
    Float y(std::size_t index, VertexId corner) const {return _vertices.at(_faces.at(index)[corner]).y;}
    Float z(std::size_t index, VertexId corner) const {return _vertices.at(_faces.at(index)[corner]).z;}
    
    Float normal_x(std::size_t index) const {return _planes.at(index).x;}
    Float normal_y(std::size_t index) const {return _planes.at(index).y;}
    Float normal_z(std::size_t index) const {return _planes.at(index).z;}
    Float d(std::size_t index) const {return _planes.at(index).d;}

    Float area(std::size_t index) const {return _areas.at(index);}

    std::size_t size() const {return _count;}

    /// Convert a face index into a face ID.
    FaceId index_to_face_id(std::size_t const index) const {return index;}

private:
    struct Vertex {
        Float x,y,z;
    };
    struct Plane {
        Float x,y,z,d;
    };
    struct Face {
        std::int32_t vertex_index[3];
        std::size_t const operator[](std::size_t const i) const {return vertex_index[i];}
    };

    std::vector<Plane>  _planes;
    std::vector<Vertex> _vertices;
    std::vector<Face>   _faces;
    std::vector<Float>  _areas;

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

