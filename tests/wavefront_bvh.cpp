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
#include <calico/accelerator/brute_force.hpp>
#include <calico/utilities/meshes/wavefront_soa.hpp>
#include <calico/utilities/meshes/wavefront_aos.hpp>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>
#include <fstream>

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
    calico::result::make_soa_result<Float, std::size_t>(face_id, t, hit_x, hit_y, hit_z);


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
    calico::result::make_soa_result<Float, std::size_t>(face_id, t, hit_x, hit_y, hit_z);


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
    calico::result::make_soa_result<Float, std::size_t>(face_id, t, hit_x, hit_y, hit_z);


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


TEST_CASE("fire a single ray against an AoS Wavefront OBJ mesh loaded from disk") {

  typedef double Float;
  typedef calico::utilities::meshes::WavefrontAoS<Float> Mesh;
  typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::SimpleBvh<Float, Mesh, Containment> BvhAccelerator;
  typedef calico::accelerator::BruteForce<Float, Mesh, Containment> BruteForce;

  std::ifstream mesh_stream("../meshes/LGPL-meshes/lamp.obj");
  // std::ifstream mesh_stream("humanoid_tri.obj");
  Mesh mesh(mesh_stream);
  BvhAccelerator accelerator(mesh);
  BruteForce brute_force(mesh);

  Float aabb[6] = {0., 0., 0.,
                   0., 0., 0.};
  mesh.bounding_box(aabb[0], aabb[1], aabb[2],
                    aabb[3], aabb[4], aabb[5]);

  // Fire a ray straight through the center of the mesh
  Mesh::FaceId brute_face_id[1] = {Mesh::ray_miss_id_c};
  Mesh::FaceId bvh_face_id[1] = {Mesh::ray_miss_id_c};
  Float start_x[1]     = { aabb[3] + 5. }; // to the right of the mesh
  Float start_y[1]     = { (aabb[4] - aabb[1]) * 0.5 + aabb[1] }; // centered in Y
  Float start_z[1]     = { (aabb[5] - aabb[2]) * 0.5 + aabb[2] }; // and in Z
  Float direction_x[1] = {-1.};
  Float direction_y[1] = { 0.};
  Float direction_z[1] = { 0.};
  auto rays = 
    calico::input::make_soa_input<Float, Mesh::FaceId>(1, bvh_face_id, 
                                         start_x, start_y, start_z,
                                         direction_x, 
                                         direction_y, 
                                         direction_z);



  // Shoot a ray using the brute-force accelerator...
  Float brute_t[1]              = {-1000.};
  Float brute_hit_x[1]          = {-1000.};
  Float brute_hit_y[1]          = {-1000.};
  Float brute_hit_z[1]          = {-1000.};
  auto brute_results = 
    calico::result::make_soa_result<Float, std::size_t>(brute_face_id, brute_t, brute_hit_x, brute_hit_y, brute_hit_z);
  auto brute_tracer = calico::make_tracer<Float>(mesh, brute_force);
  brute_tracer.trace_rays(rays, brute_results);

  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "Brute force result: Face " << brute_face_id[0] << " at (" << brute_hit_x[0] << ", " << brute_hit_y[0] << ", " << brute_hit_z[0] << ") at t=" << brute_t[0] << std::endl;

  // ... and the *same* ray using the brute-force accelerator. This _should_ be
  // *much* faster, but should match the brute-force accelerator.
  Float bvh_t[1]              = {-1000.};
  Float bvh_hit_x[1]          = {-1000.};
  Float bvh_hit_y[1]          = {-1000.};
  Float bvh_hit_z[1]          = {-1000.};
  auto bvh_results = 
    calico::result::make_soa_result<Float, std::size_t>(bvh_face_id, bvh_t, bvh_hit_x, bvh_hit_y, bvh_hit_z);
  auto bvh_tracer = calico::make_tracer<Float>(mesh, accelerator);
  bvh_tracer.trace_rays(rays, bvh_results);

  std::cout << "BVH result: Face " << bvh_face_id[0] << " at (" << bvh_hit_x[0] << ", " << bvh_hit_y[0] << ", " << bvh_hit_z[0] << ") at t=" << bvh_t[0] << std::endl;

  // Verify that both tracers hit the same facet at the same position
  REQUIRE(bool(bvh_face_id[0] == brute_face_id[0]));
  REQUIRE(bvh_hit_x[0] == doctest::Approx(brute_hit_x[0]));
  REQUIRE(bvh_hit_y[0] == doctest::Approx(brute_hit_y[0]));
  REQUIRE(bvh_hit_z[0] == doctest::Approx(brute_hit_z[0]));
}
//=============================================================================







