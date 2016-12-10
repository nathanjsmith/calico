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

#ifndef __CALICO__MATH__HPP__
#define __CALICO__MATH__HPP__

#include <cmath>

namespace calico {
namespace math {


/**
    c = A . B

    @param x1 A.x
    @param y1 A.y
    @param z1 A.z
    @param x2 B.x
    @param y2 B.y
    @param z2 B.z

    @return The dot product of A and B
*/
template <typename Float>
Float dot(const Float &x1, const Float &y1, const Float &z1, 
          const Float &x2, const Float &y2, const Float &z2)
{
    return x1*x2 + y1*y2 + z1*z2;
}
//=============================================================================


/**
    Return the length of the vector defined by (x,y,z)
*/
template <typename Float>
Float length(const Float &x, const Float &y, const Float &z) {
    return std::sqrt(dot(x,y,z,  x,y,z));
}
//=============================================================================


/**
    Right-handed C = A x B

    @param x1 A.x
    @param y1 A.y
    @param z1 A.z
    @param x2 B.x
    @param y2 B.y
    @param z2 B.z

    @param rx C.x
    @param ry C.y
    @param rz C.z
*/
template <typename FloatA, typename FloatB, typename FloatC>
void cross(const FloatA &x1, const FloatA &y1, const FloatA &z1,
           const FloatB &x2, const FloatB &y2, const FloatB &z2,
           FloatC &rx, FloatC &ry, FloatC &rz)
{
    rx = y1*z2 - z1*y2;
    ry = z1*x2 - x1*z2;
    rz = x1*y2 - y1*x2;
}
//=============================================================================


/**
  Compute the area of a triangle.

  @param x1 A.x
  @param y1 A.y
  @param z1 A.z

  @param x2 B.x
  @param y2 B.y
  @param z2 B.z

  @param x3 C.x
  @param y3 C.y
  @param z3 C.z

  area = |(B - A) x (C - A)| / 2
*/
template <typename Float>
Float area(const Float &x0, const Float &y0, const Float &z0,
           const Float &x1, const Float &y1, const Float &z1,
           const Float &x2, const Float &y2, const Float &z2)
{
  Float xx(0.), yy(0.), zz(0.);
  math::cross(x1 - x0, y1 - y0, z1 - z0,
              x2 - x0, y2 - y0, z2 - z0,
              xx, yy, zz);
  return math::length(xx, yy, zz) * Float(0.5);
}
//=============================================================================


/**
    Finds the intersection between the ray starting at (sx, sy, sz) pointing in
    the direction specified by the unit vector (dx, dy, dz) and the plane
    defined by its unit normal (nx, ny, nz) and signed distance from the origin
    d.  The distance t along the ray that the ray intersects the plane is
    returned.

    @param nx   Plane's normal vector x
    @param ny   Plane's normal vector y
    @param nz   Plane's normal vector z
    @param d    Plane's signed distance from the origin

    @param sx   Start x of the ray
    @param sy   Start y of the ray
    @param sz   Start z of the ray
    @param dx   Unit direction x of the ray
    @param dy   Unit direction y of the ray
    @param dz   Unit direction z of the ray

    @return     Distance along the ray vector that the ray and plane intersect
                one another.
*/
template <typename Float>
Float ray_plane_intersection(
        const Float &nx, const Float &ny, const Float &nz, const Float &d,
        const Float &sx, const Float &sy, const Float &sz,
        const Float &dx, const Float &dy, const Float &dz)
{
    return -(dot(sx, sy, sz, nx, ny, nz) + d) / dot(dx, dy, dz, nx, ny, nz);
}
//=============================================================================


template <typename Float, typename Mesh>
class MollerTrumboreContainmentTest
{
public:
    static
    bool contains(const Mesh &mesh, const typename Mesh::FaceId &face,
                  const Float &p_x, const Float &p_y, const Float &p_z)
    {
        // barycentric coordinates are defined by the u = area(A,B,P) / area(A,B,C)
        //                                            v = area(A,C,P) / area(A,B,C)
        //                                            w = area(B,C,P) / area(A,B,C)

        // Translate the triangle to be at the origin
        Float bx = mesh.x(face, 1) - mesh.x(face, 0);
        Float by = mesh.y(face, 1) - mesh.y(face, 0);
        Float bz = mesh.z(face, 1) - mesh.z(face, 0);

        Float cx = mesh.x(face, 2) - mesh.x(face, 0);
        Float cy = mesh.y(face, 2) - mesh.y(face, 0);
        Float cz = mesh.z(face, 2) - mesh.z(face, 0);

        Float px = p_x - mesh.x(face, 0);
        Float py = p_y - mesh.y(face, 0);
        Float pz = p_z - mesh.z(face, 0);

        // Get the area of the triangle we're testing and double it to make the
        // following calculations faster.
        Float inverse_double_area = Float(1) / (mesh.area(face) * Float(2));

        // Find the area of triangle (A, B, P) and scale it by the area of the
        // whole parent triangle.  We use double the area because the length of the
        // cross-product is actually the area of the parallelogram, which is twice
        // the triangle area.
        Float tmp_x, tmp_y, tmp_z;
        cross(bx, by, bz, px, py, pz, tmp_x, tmp_y, tmp_z);
        Float u = length(tmp_x, tmp_y, tmp_z) * inverse_double_area;

        // Find the area of triangle (A, C, P) and scale it by the area of the
        // whole parent triangle.  We use double the area because the length of the
        // cross-product is actually the area of the parallelogram, which is twice
        // the triangle area.
        cross(cx, cy, cz, px, py, pz, tmp_x, tmp_y, tmp_z);
        Float v = length(tmp_x, tmp_y, tmp_z) * inverse_double_area;

        Float w = Float(1) - u - v;
        return (Float(0) <= u && u <= Float(1) &&
                Float(0) <= v && v <= Float(1) &&
                Float(0) <= w && w <= Float(1));
    }
};
//=============================================================================



/**
    The same as MollerTrumboreContainmentTest except that it computes the
    value w from the point and the triangle corners rather than from u and v.
*/
template <typename Mesh, typename Float>
class PluckerContainmentTest
{
public:
    static
    bool contains(const Mesh &mesh, const typename Mesh::FaceId &face,
                  const Float &px, const Float &py, const Float &pz)
    {
        // barycentric coordinates are defined by the u = area(A,B,P) / area(A,B,C)
        //                                            v = area(B,C,P) / area(A,B,C)
        //                                            w = area(C,A,P) / area(A,B,C)

        // Treat the point as the origin and translate the three corners
        // accordingly
        const Float ax = mesh.corner_a_x(face) - px;
        const Float ay = mesh.corner_a_y(face) - py;
        const Float az = mesh.corner_a_z(face) - pz;

        const Float bx = mesh.corner_b_x(face) - px;
        const Float by = mesh.corner_b_y(face) - py;
        const Float bz = mesh.corner_b_z(face) - pz;

        const Float cx = mesh.corner_c_x(face) - px;
        const Float cy = mesh.corner_c_y(face) - py;
        const Float cz = mesh.corner_c_z(face) - pz;

        // Get the area of the triangle we're testing and double it to make the
        // following calculations faster.
        const Float inverse_double_area = Float(1) / (mesh.area(face) * Float(2));

        // Find the area of all three sub-triangles with P as a vertex.  Compute
        // the ratio of the sub-triangle areas against the whole triangle's area.
        // This gives us the u, v and w coordinates (i.e. barrycentric coordinates)
        // We use double the area because the length of the cross-product is
        // actually the area of the parallelogram, which is twice the triangle
        // area.
        Float tmp_x, tmp_y, tmp_z;
        cross(ax, ay, az, bx, by, bz, tmp_x, tmp_y, tmp_z);
        const Float u = length(tmp_x, tmp_y, tmp_z) * inverse_double_area;

        cross(bx, by, bz, cx, cy, cz, tmp_x, tmp_y, tmp_z);
        const Float v = length(tmp_x, tmp_y, tmp_z) * inverse_double_area;

        cross(cx, cy, cz, ax, ay, az, tmp_x, tmp_y, tmp_z);
        const Float w = length(tmp_x, tmp_y, tmp_z) * inverse_double_area;

        return (Float(0) <= u && u <= Float(1) &&
                Float(0) <= v && v <= Float(1) &&
                Float(0) <= w && w <= Float(1));
    }
};
//=============================================================================


/**
    Point in triangle test that _always_ returns True.
*/
template <typename Float>
bool true_point_in_triangle(
        const Float, const Float, const Float,
        const Float, const Float, const Float,
        const Float, const Float, const Float,
        const Float, const Float, const Float)
{
    return true;
}
//=============================================================================


/**
    Point in triangle test that _always_ returns False.
*/
template <typename Float>
bool false_point_in_triangle(
        const Float, const Float, const Float,
        const Float, const Float, const Float,
        const Float, const Float, const Float,
        const Float, const Float, const Float)
{
    return false;
}
//=============================================================================



}// end namespace math
}// end namespace calico

#endif

