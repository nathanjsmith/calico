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

BOOST_AUTO_TEST_CASE(hard_coded_cases) {
  typedef double Float;

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

    // pre-compute the inverse of the direction -- internals of intersects
    // should be robust to divide-by-zero
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
    Float min_t{0};
    Float max_t = calico::accelerator::StdTypeInterface<Float>::max_infinity();
    bool safe = 
      calico::accelerator::intersects(start_x[i], start_y[i], start_z[i],
                                      direction_x[i], direction_y[i], direction_z[i],
                                      inverse_direction_x[i], 
                                      inverse_direction_y[i], 
                                      inverse_direction_z[i],
                                      min_x[0], min_y[0], min_z[0],
                                      max_x[0], max_y[0], max_z[0],
                                      min_t, max_t);

    if (safe != expected[i]) {
        std::cerr << "Error: We expected that ray " << i << " would "
                  << (expected[i] ? "hit" : "miss") << " the bounding "
                  << "box, but it didn't." << std::endl;
        std::cerr << "Ray (" << start_x[i] << ", " << start_y[i] << ", " << start_z[i] << ")\n"
                  << "    (" << direction_x[i] << ", " << direction_y[i] << ", " << direction_z[i] << ")\n"
                  << "AABB [(" << min_x[j] << ", " << min_y[j] << ", " << min_z[j] << "),\n"
                  << "      (" << max_x[j] << ", " << max_y[j] << ", " << max_z[j] << ")\n"
                  << "     ]" << std::endl;
    }
    BOOST_REQUIRE_EQUAL(expected[i], safe);
  }
}


BOOST_AUTO_TEST_CASE(axis_aligned_hits) {
  typedef double Float;

  std::mt19937 rng;
  std::uniform_real_distribution<Float> uniform_backoff( 1., 10000.);
  std::uniform_real_distribution<Float> uniform_slide  (-1., 1.);

  auto generate_backoff = std::bind(uniform_backoff, rng);
  auto generate_slide = std::bind(uniform_slide, rng);

  const Float min_x{-1};
  const Float min_y{-1};
  const Float min_z{-1};

  const Float max_x{1};
  const Float max_y{1};
  const Float max_z{1};
  
  // Try rays in the x-axis
  const std::size_t axes[] = {0, 1, 2};
  const Float directions[] = {1., -1.};
  bool hit;
  for (const std::size_t axis : axes) {
    for (Float direction : directions) {
      for (std::size_t i = 0u; i < 10000u; ++i) {
        Float start_x = -direction * generate_backoff();
        Float start_y = generate_slide();
        Float start_z = generate_slide();

        Float direction_x{0};
        Float direction_y{0};
        Float direction_z{0};

        if (axis == 0) {
          direction_x = direction;
          start_x = -direction * generate_backoff();
          start_y = generate_slide();
          start_z = generate_slide();
        }
        else if (axis == 1) {
          direction_y = direction;
          start_x = generate_slide();
          start_y = -direction * generate_backoff();
          start_z = generate_slide();
        }
        else if (axis == 2) {
          direction_z = direction;
          start_x = generate_slide();
          start_y = generate_slide();
          start_z = -direction * generate_backoff();
        }

        const Float inverse_direction_x{Float(1)/direction_x};
        const Float inverse_direction_y{Float(1)/direction_y};
        const Float inverse_direction_z{Float(1)/direction_z};

        Float min_t{0};
        Float max_t{calico::accelerator::StdTypeInterface<Float>::max_infinity()};

        hit = calico::accelerator::intersects(start_x, start_y, start_z,
                                              direction_x, direction_y, direction_z,
                                              inverse_direction_x, inverse_direction_y, inverse_direction_z,
                                              min_x, min_y, min_z, 
                                              max_x, max_y, max_z,
                                              min_t, max_t);
        if (!hit) {
          std::cerr << "Ray (" << start_x << ", " << start_y << ", " << start_z << ") -> (" << direction_x << ", " << direction_y << ", " << direction_z << ") missed" << std::endl;
        }
        BOOST_REQUIRE(hit && "Failed intersection for index");
      }
    }
  }

}


/**
  This unit test finds a ponit on each of a bounding box's six faces and then
  fires a ray towards that point. All of the fired rays should then have an
  intersection with the bounding box.
*/
BOOST_AUTO_TEST_CASE(generalized_hits) {
  typedef double Float;

  std::mt19937 rng;
  std::uniform_real_distribution<Float> uniform_neg_pos(-1., 1.);
  std::uniform_real_distribution<Float> uniform_scalar(-1., 1.);

  auto generate_relative_position = std::bind(uniform_neg_pos, rng);
  auto generate_source_position = std::bind(uniform_scalar, rng);

  const Float min_x{-1};
  const Float min_y{-1};
  const Float min_z{-1};

  const Float max_x{1};
  const Float max_y{1};
  const Float max_z{1};
  
  // Try rays in the x-axis
  const std::size_t sides[] = {0, 1, 2, 3, 4, 5};
  for (const std::size_t side : sides) {
    for (std::size_t i = 0u; i < 10000u; ++i) {
      // Generate a point on the surface of this side. We'll fire a ray towards
      // that point.
      Float target_u = generate_relative_position();
      Float target_v = generate_relative_position();

      Float target_x,target_y,target_z,start_x,start_y,start_z;
      if (side == 0) {
        target_x = 1.;
        target_y = target_u;
        target_z = target_v;
        start_x = target_x + generate_source_position() * 100.;
        start_y = generate_source_position() * 100.;
        start_z = generate_source_position() * 100.;
      }
      else if (side == 1) {
        target_x = -1.;
        target_y = target_u;
        target_z = target_v;
        start_x = target_x - generate_source_position() * 100.;
        start_y = generate_source_position() * 100.;
        start_z = generate_source_position() * 100.;
      }
      else if (side == 2) {
        target_x = target_u;
        target_y = 1;
        target_z = target_v;
        start_x = generate_source_position() * 100.;
        start_y = target_y + generate_source_position() * 100.;
        start_z = generate_source_position() * 100.;
      }
      else if (side == 3) {
        target_x = target_u;
        target_y = -1;
        target_z = target_v;
        start_x = generate_source_position() * 100.;
        start_y = target_y - generate_source_position() * 100.;
        start_z = generate_source_position() * 100.;
      }
      else if (side == 4) {
        target_x = target_u;
        target_y = target_v;
        target_z = 1;
        start_x = generate_source_position() * 100.;
        start_y = generate_source_position() * 100.;
        start_z = target_z + generate_source_position() * 100.;
      }
      else {
        target_x = target_u;
        target_y = target_v;
        target_z = -1;
        start_x = generate_source_position() * 100.;
        start_y = generate_source_position() * 100.;
        start_z = target_z - generate_source_position() * 100.;
      }

      Float direction_x{target_x - start_x};
      Float direction_y{target_y - start_y};
      Float direction_z{target_z - start_z};

      calico::math::normalize(direction_x, direction_y, direction_z);

      const Float inverse_direction_x{Float(1)/direction_x};
      const Float inverse_direction_y{Float(1)/direction_y};
      const Float inverse_direction_z{Float(1)/direction_z};

      Float min_t{0};
      Float max_t{calico::accelerator::StdTypeInterface<Float>::max_infinity()};

      bool hit =
            calico::accelerator::intersects(start_x, start_y, start_z,
                                            direction_x, direction_y, direction_z,
                                            inverse_direction_x, inverse_direction_y, inverse_direction_z,
                                            min_x, min_y, min_z, 
                                            max_x, max_y, max_z,
                                            min_t, max_t);
      if (!hit) {
        std::cerr << "Ray (" << start_x << ", " << start_y << ", " << start_z << ") -> (" << direction_x << ", " << direction_y << ", " << direction_z << ") missed" << std::endl;
      }
      BOOST_REQUIRE(hit && "Failed intersection for index");
    }
  }

}



/**
  Generates points that lie above, to the side-of, or below the bounding box,
  fires a ray towards those points, and verifies that they missed.
*/
BOOST_AUTO_TEST_CASE(generalized_misses) {
  typedef double Float;

  std::mt19937 rng;
  std::uniform_real_distribution<Float> uniform_before(-1., -2.);
  std::uniform_real_distribution<Float> uniform_within( 0.,  1.);
  std::uniform_real_distribution<Float> uniform_beyond( 1.,  2.);
  std::uniform_real_distribution<Float> uniform_scalar(-1., 1.);
  auto generate_source_position = std::bind(uniform_scalar, rng);

  const Float min_x{-1};
  const Float min_y{-1};
  const Float min_z{-1};

  const Float max_x{1};
  const Float max_y{1};
  const Float max_z{1};
  
  // Try rays in the x-axis
  const std::size_t sides[] = {0, 1, 2, 3, 4, 5};
  const std::size_t above_below_left_right[] = {0, 1, 2, 3};

  for (const std::size_t side : sides) {
    for (const std::size_t ablr : above_below_left_right) {
      std::uniform_real_distribution<Float> *u_generator;
      std::uniform_real_distribution<Float> *v_generator;
      if (ablr == 0u) {
        // above
        u_generator = &uniform_within;
        v_generator = &uniform_beyond;
      }
      else if (ablr == 1u) {
        // left
        u_generator = &uniform_before;
        v_generator = &uniform_within;
      }
      else if (ablr == 2u) {
        // right
        u_generator = &uniform_beyond;
        v_generator = &uniform_within;
      }
      else {
        // below
        u_generator = &uniform_within;
        v_generator = &uniform_before;
      }

      for (std::size_t i = 0u; i < 10000u; ++i) {
        // Generate a point on the surface of this side. We'll fire a ray towards
        // that point.
        Float target_u = (*u_generator)(rng);
        Float target_v = (*v_generator)(rng);

        Float target_x,target_y,target_z,start_x,start_y,start_z;
        if (side == 0) {
          target_x = 1.;
          target_y = target_u;
          target_z = target_v;
          start_x = target_x + 1 + generate_source_position();
          start_y = generate_source_position();
          start_z = generate_source_position();
        }
        else if (side == 1) {
          target_x = -1.;
          target_y = target_u;
          target_z = target_v;
          start_x = target_x - 1 - generate_source_position();
          start_y = generate_source_position();
          start_z = generate_source_position();
        }
        else if (side == 2) {
          target_x = target_u;
          target_y = 1;
          target_z = target_v;
          start_x = generate_source_position();
          start_y = target_y + 1 + generate_source_position();
          start_z = generate_source_position();
        }
        else if (side == 3) {
          target_x = target_u;
          target_y = -1;
          target_z = target_v;
          start_x = generate_source_position();
          start_y = target_y - 1 - generate_source_position();
          start_z = generate_source_position();
        }
        else if (side == 4) {
          target_x = target_u;
          target_y = target_v;
          target_z = 1;
          start_x = generate_source_position();
          start_y = generate_source_position();
          start_z = target_z + 1 + generate_source_position();
        }
        else {
          target_x = target_u;
          target_y = target_v;
          target_z = -1;
          start_x = generate_source_position();
          start_y = generate_source_position();
          start_z = target_z - 1 - generate_source_position();
        }

        Float direction_x{target_x - start_x};
        Float direction_y{target_y - start_y};
        Float direction_z{target_z - start_z};

        calico::math::normalize(direction_x, direction_y, direction_z);

        const Float inverse_direction_x{Float(1)/direction_x};
        const Float inverse_direction_y{Float(1)/direction_y};
        const Float inverse_direction_z{Float(1)/direction_z};

        Float min_t{0};
        Float max_t{calico::accelerator::StdTypeInterface<Float>::max_infinity()};

        bool hit =
              calico::accelerator::intersects(start_x, start_y, start_z,
                                              direction_x, direction_y, direction_z,
                                              inverse_direction_x, inverse_direction_y, inverse_direction_z,
                                              min_x, min_y, min_z, 
                                              max_x, max_y, max_z,
                                              min_t, max_t);
        if (hit) {
          std::cerr << "Ray (" << start_x << ", " << start_y << ", " << start_z << ") -> (" 
                               << direction_x << ", " << direction_y << ", " << direction_z 
                               << ") hit, but should have missed" << std::endl;
          std::cerr << "Target (" << target_x << ", " << target_y << ", " << target_z << ")" << std::endl;
          std::cerr << "Above/below/left/right == " << ablr << std::endl;
          std::cerr << "Side: " << side << std::endl;
        }
        BOOST_REQUIRE(!hit && "Expected ray to miss the bounding box");
      }
    }
  }

}



