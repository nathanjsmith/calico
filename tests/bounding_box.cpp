// Copyright (c) 2017, Nathan Smith <nathanjsmith@gmail.com>
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

#include <calico/accelerator/bounding_box.hpp>

#include <random>
#include <iostream>

#define BOOST_TEST_MODULE Bounding Box Tests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(compare_safe_and_unsafe_intersection_routines) {
  typedef double Float;

  std::mt19937 rng;
  std::uniform_real_distribution<Float> uniform_start(-2., 2.);
  std::uniform_real_distribution<Float> uniform_direction(-1., 1.);

  auto generate_start = std::bind(uniform_start, rng);
  auto generate_direction = std::bind(uniform_direction, rng);

  for (size_t round = 0u; round < 1000u; ++round) {
    const size_t ray_count(1000);
    Float start_x[ray_count];
    Float start_y[ray_count];
    Float start_z[ray_count];

    Float direction_x[ray_count];
    Float direction_y[ray_count];
    Float direction_z[ray_count];

    Float inverse_direction_x[ray_count];
    Float inverse_direction_y[ray_count];
    Float inverse_direction_z[ray_count];

    // Make a few pathological rays
    size_t pathological_ray_count(4);
    start_x[0]     = Float(0);
    start_y[0]     = Float(0);
    start_z[0]     = Float(0);
    direction_x[0] = Float(1);
    direction_y[0] = Float(0);
    direction_z[0] = Float(0);

    start_x[1]     = Float(-1);
    start_y[1]     = Float(0);
    start_z[1]     = Float(0);
    direction_x[1] = Float(1);
    direction_y[1] = Float(0);
    direction_z[1] = Float(0);

    start_x[2]     = Float(2);
    start_y[2]     = Float(0);
    start_z[2]     = Float(0);
    direction_x[2] = Float(-1);
    direction_y[2] = Float(0);
    direction_z[2] = Float(0);

    start_x[3]     = Float(-1);
    start_y[3]     = Float(0.5);
    start_z[3]     = Float(0);
    direction_x[3] = Float(-1);
    direction_y[3] = Float(0);
    direction_z[3] = Float(0);

    for (size_t i = 0u; i < ray_count; ++i) {
      if (i >= pathological_ray_count) {
        start_x[i]     = generate_start();
        start_y[i]     = generate_start();
        start_z[i]     = generate_start();

        direction_x[i] = generate_direction();
        direction_y[i] = generate_direction();
        direction_z[i] = generate_direction();
      }

      // normalize the direction
      calico::math::normalize(direction_x[i], direction_y[i], direction_z[i]);

      inverse_direction_x[i] = Float(1) / direction_x[i];
      inverse_direction_y[i] = Float(1) / direction_y[i];
      inverse_direction_z[i] = Float(1) / direction_z[i];
    }

    const size_t aabb_count(1);
    Float min_x[1] = {0.};
    Float min_y[1] = {0.};
    Float min_z[1] = {0.};

    Float max_x[1] = {1.};
    Float max_y[1] = {1.};
    Float max_z[1] = {1.};

    for (size_t j = 0u; j < aabb_count; ++j) {
      for (size_t i = 0u; i < ray_count; ++i) {
        bool safe = 
          calico::accelerator::intersects(start_x[i], start_y[i], start_z[i],
                                          direction_x[i], direction_y[i], direction_z[i],
                                          min_x[0], min_y[0], min_z[0],
                                          max_x[0], max_y[0], max_z[0]);
        bool unsafe = 
          calico::accelerator::unsafe_intersects(
                                          start_x[i], start_y[i], start_z[i],
                                          inverse_direction_x[i], inverse_direction_y[i], inverse_direction_z[i],
                                          min_x[0], min_y[0], min_z[0],
                                          max_x[0], max_y[0], max_z[0]);
        if (safe != unsafe) {
            std::cerr << "Ray (" << start_x[i] << ", " << start_y[i] << ", " << start_z[i] << ") -> ("
                                 << direction_x[i] << ", " << direction_y[i] << ", " << direction_z[i] << ") "
                      << "intersecting AABB [(" << min_x[j] << ", " << min_y[j] << ", " << min_z[j] << "), ("
                      << max_x[j] << ", " << max_y[j] << ", " << max_z[j] << ")] disagreed between safe and "
                      << "unsafe implementations (" << safe << " vs " << unsafe << ")" << std::endl;
        }
        BOOST_REQUIRE_EQUAL(safe, unsafe);
      }
    }
  }
}
//=============================================================================


BOOST_AUTO_TEST_CASE(misses) {
  typedef double Float;

  std::mt19937 rng;
  std::uniform_real_distribution<Float> uniform_start(-2., 2.);
  std::uniform_real_distribution<Float> uniform_direction(-1., 1.);

  auto generate_start = std::bind(uniform_start, rng);
  auto generate_direction = std::bind(uniform_direction, rng);

  const size_t ray_count(6);
  Float start_x[ray_count];
  Float start_y[ray_count];
  Float start_z[ray_count];

  Float direction_x[ray_count];
  Float direction_y[ray_count];
  Float direction_z[ray_count];

  Float inverse_direction_x[ray_count];
  Float inverse_direction_y[ray_count];
  Float inverse_direction_z[ray_count];

  bool expected[ray_count];

  // Start inside the box. This should return true
  start_x[0]     = Float(0.5);
  start_y[0]     = Float(0.5);
  start_z[0]     = Float(0.5);
  direction_x[0] = Float(1);
  direction_y[0] = Float(0);
  direction_z[0] = Float(0);
  expected[0]    = true;

  // Start on the near wall of a bounding box.
  start_x[1]     = Float(0);
  start_y[1]     = Float(0.5);
  start_z[1]     = Float(0.5);
  direction_x[1] = Float(1);
  direction_y[1] = Float(0);
  direction_z[1] = Float(0);
  expected[1]    = true;

  // Start on the far wall of the box
  start_x[2]     = Float(1);
  start_y[2]     = Float(0.5);
  start_z[2]     = Float(0.5);
  direction_x[2] = Float(1);
  direction_y[2] = Float(0);
  direction_z[2] = Float(0);
  expected[2]    = true;

  // Start outside the box
  start_x[3]     = Float(-1);
  start_y[3]     = Float(0.5);
  start_z[3]     = Float(0.5);
  direction_x[3] = Float(1);
  direction_y[3] = Float(0);
  direction_z[3] = Float(0);
  expected[3]    = true;

  // Start outside the near-side of the box, and point away from the box
  start_x[4]     = Float(-1);
  start_y[4]     = Float(0.5);
  start_z[4]     = Float(0.5);
  direction_x[4] = Float(-1);
  direction_y[4] = Float(0);
  direction_z[4] = Float(0);
  expected[4]    = false;

  // Start outside the far-side of the box, and point away from the box
  start_x[5]     = Float(1.01);
  start_y[5]     = Float(0.5);
  start_z[5]     = Float(0.5);
  direction_x[5] = Float(1);
  direction_y[5] = Float(0);
  direction_z[5] = Float(0);
  expected[5]    = false;

  for (size_t i = 0u; i < ray_count; ++i) {
    // normalize the direction
    calico::math::normalize(direction_x[i], direction_y[i], direction_z[i]);

    // pre-compute the inverse of the direction
    inverse_direction_x[i] = Float(1) / direction_x[i];
    inverse_direction_y[i] = Float(1) / direction_y[i];
    inverse_direction_z[i] = Float(1) / direction_z[i];
  }

  const size_t aabb_count(1);
  Float min_x[1] = {0.};
  Float min_y[1] = {0.};
  Float min_z[1] = {0.};

  Float max_x[1] = {1.};
  Float max_y[1] = {1.};
  Float max_z[1] = {1.};

  const size_t j(0);
  for (size_t i = 0u; i < ray_count; ++i) {
    bool safe = 
      calico::accelerator::intersects(start_x[i], start_y[i], start_z[i],
                                      direction_x[i], direction_y[i], direction_z[i],
                                      min_x[0], min_y[0], min_z[0],
                                      max_x[0], max_y[0], max_z[0]);
    bool unsafe = 
      calico::accelerator::unsafe_intersects(
                                      start_x[i], start_y[i], start_z[i],
                                      inverse_direction_x[i], inverse_direction_y[i], inverse_direction_z[i],
                                      min_x[0], min_y[0], min_z[0],
                                      max_x[0], max_y[0], max_z[0]);
    if (safe != unsafe) {
        std::cerr << "Ray (" << start_x[i] << ", " << start_y[i] << ", " << start_z[i] << ") -> ("
                             << direction_x[i] << ", " << direction_y[i] << ", " << direction_z[i] << ") "
                  << "intersecting AABB [(" << min_x[j] << ", " << min_y[j] << ", " << min_z[j] << "), ("
                  << max_x[j] << ", " << max_y[j] << ", " << max_z[j] << ")] disagreed between safe and "
                  << "unsafe implementations (" << safe << " vs " << unsafe << ")" << std::endl;
    }
    BOOST_REQUIRE_EQUAL(expected[i], safe);
    BOOST_REQUIRE_EQUAL(safe, unsafe);
  }
}
