// Copyright (c) 2016-2017, Nathan Smith <nathanjsmith@gmail.com>
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

#include <calico/tracer.hpp>
#include <calico/math.hpp>
#include <calico/input/soa_input.hpp>
#include <calico/result/soa_result.hpp>
#include <calico/accelerator/simple_bvh.hpp>
#include <calico/utilities/meshes/wavefront_soa.hpp>
#include <calico/utilities/meshes/wavefront_aos.hpp>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>

TEST_CASE("fire a single ray against a Wavefront OBJ loaded mesh too small to subdivide; Find the intersection using Plucker") {

  typedef double Float;
  typedef calico::utilities::meshes::WavefrontSoA<Float> Mesh;
  typedef calico::math::PluckerContainmentTest<Float, Mesh> Containment;
  // typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::SimpleBvh<Float, Mesh, Containment> Accelerator;

  std::stringstream plate_string;
  plate_string << "v -1 -1  3\n"
               << "v  1 -1  3\n"
               << "v  1  1  3\n"
               << "v -1  1  3\n"
               << "\n"
               << "f  1 2 3 4\n";

  Mesh plate(plate_string);
  Accelerator accelerator(plate);

  auto tracer = calico::make_tracer<Float>(plate, accelerator);

  Mesh::FaceId face_id[1] = {Mesh::ray_miss_id_c};
  Float start_x[1]     = { 0.};
  Float start_y[1]     = { 0.};
  Float start_z[1]     = {10.};
  Float direction_x[1] = { 0.};
  Float direction_y[1] = { 0.};
  Float direction_z[1] = {-1.};
  auto rays = 
    calico::input::make_soa_input<Float, Mesh::FaceId>(1, face_id, 
                                         start_x, start_y, start_z,
                                         direction_x, direction_y, direction_z);


  Float t[1]              = {-1000.};
  Float hit_x[1]          = {-1000.};
  Float hit_y[1]          = {-1000.};
  Float hit_z[1]          = {-1000.};
  auto results = 
    calico::result::make_soa_result<Float, Mesh::FaceId>(face_id, t, hit_x, hit_y, hit_z);


  tracer.trace_rays(rays, results);

  // Verify the ray struck the surface we expected in the location we expected.
  // It's actually ambiguous whether the ray will strike facet 0 or 1. The ray
  // should strike along the edge shared by the two facets.
  REQUIRE(bool(face_id[0] == 0u || face_id[0] == 1u));
  REQUIRE(hit_x[0] == doctest::Approx(0.));
  REQUIRE(hit_y[0] == doctest::Approx(0.));
  REQUIRE(hit_z[0] == doctest::Approx(3.));
}
//=============================================================================


TEST_CASE("fire a single ray against a Wavefront OBJ loaded mesh too small to subdivide; Find the intersection using Moller-Trumbore") {

  typedef double Float;
  typedef calico::utilities::meshes::WavefrontSoA<Float> Mesh;
  typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::SimpleBvh<Float, Mesh, Containment> Accelerator;

  std::stringstream plate_string;
  plate_string << "v -1 -1  3\n"
               << "v  1 -1  3\n"
               << "v  1  1  3\n"
               << "v -1  1  3\n"
               << "\n"
               << "f  1 2 3 4\n";

  Mesh plate(plate_string);
  Accelerator accelerator(plate);

  auto tracer = calico::make_tracer<Float>(plate, accelerator);

  Mesh::FaceId face_id[1] = {Mesh::ray_miss_id_c};
  Float start_x[1]     = { 0.};
  Float start_y[1]     = { 0.};
  Float start_z[1]     = {10.};
  Float direction_x[1] = { 0.};
  Float direction_y[1] = { 0.};
  Float direction_z[1] = {-1.};
  auto rays = 
    calico::input::make_soa_input<Float, Mesh::FaceId>(1, face_id, 
                                         start_x, start_y, start_z,
                                         direction_x, direction_y, direction_z);


  Float t[1]              = {-1000.};
  Float hit_x[1]          = {-1000.};
  Float hit_y[1]          = {-1000.};
  Float hit_z[1]          = {-1000.};
  auto results = 
    calico::result::make_soa_result<Float, Mesh::FaceId>(face_id, t, hit_x, hit_y, hit_z);


  tracer.trace_rays(rays, results);

  // Verify the ray struck the surface we expected in the location we expected.
  // It's actually ambiguous whether the ray will strike facet 0 or 1. The ray
  // should strike along the edge shared by the two facets.
  REQUIRE(bool(face_id[0] == 0u || face_id[0] == 1u));
  REQUIRE(hit_x[0] == doctest::Approx(0.));
  REQUIRE(hit_y[0] == doctest::Approx(0.));
  REQUIRE(hit_z[0] == doctest::Approx(3.));
}
//=============================================================================


TEST_CASE("fire a single ray against an AoS Wavefront OBJ loaded mesh too small to subdivide; Find the intersection using Moller-Trumbore") {

  typedef double Float;
  typedef calico::utilities::meshes::WavefrontAoS<Float> Mesh;
  typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::SimpleBvh<Float, Mesh, Containment> Accelerator;

  std::stringstream plate_string;
  plate_string << "v -1 -1  3\n"
               << "v  1 -1  3\n"
               << "v  1  1  3\n"
               << "v -1  1  3\n"
               << "\n"
               << "f  1 2 3 4\n";

  Mesh plate(plate_string);
  Accelerator accelerator(plate);

  auto tracer = calico::make_tracer<Float>(plate, accelerator);

  Mesh::FaceId face_id[1] = {Mesh::ray_miss_id_c};
  Float start_x[1]     = { 0.};
  Float start_y[1]     = { 0.};
  Float start_z[1]     = {10.};
  Float direction_x[1] = { 0.};
  Float direction_y[1] = { 0.};
  Float direction_z[1] = {-1.};
  auto rays = 
    calico::input::make_soa_input<Float, Mesh::FaceId>(1, face_id, 
                                         start_x, start_y, start_z,
                                         direction_x, direction_y, direction_z);


  Float t[1]              = {-1000.};
  Float hit_x[1]          = {-1000.};
  Float hit_y[1]          = {-1000.};
  Float hit_z[1]          = {-1000.};
  auto results = 
    calico::result::make_soa_result<Float, Mesh::FaceId>(face_id, t, hit_x, hit_y, hit_z);


  tracer.trace_rays(rays, results);

  // Verify the ray struck the surface we expected in the location we expected.
  // It's actually ambiguous whether the ray will strike facet 0 or 1. The ray
  // should strike along the edge shared by the two facets.
  REQUIRE(bool(face_id[0] == 0u || face_id[0] == 1u));
  REQUIRE(hit_x[0] == doctest::Approx(0.));
  REQUIRE(hit_y[0] == doctest::Approx(0.));
  REQUIRE(hit_z[0] == doctest::Approx(3.));
}
//=============================================================================







