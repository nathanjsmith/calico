// Copyright (c) 2016, Nathan Smith <nathanjsmith@gmail.com>
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
#include <calico/accelerator/brute_force.hpp>

#include <calico/utilities/meshes/plate.hpp>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>

TEST_CASE("fire a single ray against a flat 1x1 plate") {

  typedef double Float;
  typedef calico::utilities::meshes::Plate<Float> Mesh;
  typedef calico::math::PluckerContainmentTest<Float, Mesh> Containment;
  // typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::BruteForce<Float, Mesh, Containment> Accelerator;

  Mesh plate;
  Accelerator accelerator(plate);

  auto tracer = calico::make_tracer(accelerator);

  Mesh::FaceId start_face_id[1] = {Mesh::ray_miss_id_c};
  Float start_x[1]     = {0.};
  Float start_y[1]     = {0.};
  Float start_z[1]     = {10.}; // start above the plate
  Float direction_x[1] = {0.};
  Float direction_y[1] = {0.};
  Float direction_z[1] = {-1.}; // fire down towards the plate
  auto rays = 
    calico::input::make_soa_input<Float, Mesh::FaceId>
            (1, start_face_id, 
             start_x, start_y, start_z,
             direction_x, direction_y, direction_z);


  Mesh::FaceId strike_face_id[1] = {Mesh::ray_miss_id_c};
  Float t[1]                     = {-1000.};
  Float hit_x[1]                 = {-1000.};
  Float hit_y[1]                 = {-1000.};
  Float hit_z[1]                 = {-1000.};

  // Create a SoaResults object that adapts (provides a view into) the
  // strike_face_id, t, hit_x, hit_y, and hit_z arrays.
  auto results = 
    calico::result::make_soa_result<Float, Mesh::FaceId>(strike_face_id, t, hit_x, hit_y, hit_z);


  // Trace our 1 ray against the plate and record the results through the
  // "results" adapter into the strike_face_id, t, and hit_* arrays.
  tracer.trace_rays(rays, results);

  // Verify the ray struck the surface we expected in the location we expected.
  // It's actually ambiguous whether the ray will strike facet 0 or 1. The ray
  // should strike along the edge shared by the two facets.
  REQUIRE(bool(strike_face_id[0] == 0u || strike_face_id[0] == 1u));
  REQUIRE(hit_x[0] == doctest::Approx(0.));
  REQUIRE(hit_y[0] == doctest::Approx(0.));
  REQUIRE(hit_z[0] == doctest::Approx(3.));
}
//=============================================================================


