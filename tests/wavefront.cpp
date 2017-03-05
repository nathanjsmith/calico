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
#include <calico/accelerator/brute_force.hpp>

#include <canned_meshes/plate.hpp>

#include <iostream>

int main(int argc, const char *argv[]) {

  typedef double Float;
  typedef calico::canned_meshes::Plate<Float> Mesh;
  typedef calico::math::PluckerContainmentTest<Float, Mesh> Containment;
  // typedef calico::math::MollerTrumboreContainmentTest<Float, Mesh> Containment;
  typedef calico::accelerator::BruteForce<Float, Mesh, Containment> Accelerator;

  Mesh plate;
  Accelerator accelerator(plate);

  auto tracer = calico::make_tracer<Float>(plate, accelerator);

  Float start_x[1]     = {0.};
  Float start_y[1]     = {0.};
  Float start_z[1]     = {10.};
  Float direction_x[1] = {0.};
  Float direction_y[1] = {0.};
  Float direction_z[1] = {-1.};
  auto rays = 
    calico::input::make_soa_input<Float>(1, start_x, start_y, start_z,
                                         direction_x, direction_y, direction_z);


  Mesh::FaceId face_id[1] = {0};
  Float t[1]              = {-1000.};
  Float hit_x[1]          = {-1000.};
  Float hit_y[1]          = {-1000.};
  Float hit_z[1]          = {-1000.};
  auto results = 
    calico::result::make_soa_result<Float, std::size_t>(face_id, t, hit_x, hit_y, hit_z);


  tracer.trace_rays(rays, results);


  std::cerr << "Hit facet " << face_id[0] << " at (" 
            << hit_x[0] << ", " << hit_y[0] << ", " << hit_z[0] 
            << ")" << std::endl;

  return 0;
}
//=============================================================================


