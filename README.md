# Calico

**Calico is under development and does not currently function, except in the most basic manner.**

Calico is a ray/geometry intersection library for building ray-traced applications.
Its focus is on presenting an API that works with your existing datastructures through use of the [adapter design pattern](http://sourcemaking.com/design_patterns/adapter).
Because it is implemented using templates, Calico is a header-only library with nothing to compile.
You may include Calico in your project by adding it to your include path, and including its headers.

Calico is a ray-casting engine, which can be used to build a ray-tracer.
For example, it focuses on answering the question of "what surfaces did this ray hit?".
It does not provide reflection or transmission calculations, lighting calculations, etc.
The physics simulation is left to the implementer as they are the pieces that are application specific.

## Uses

Ideas for how Calico can be used include:

1. the basis of a ray-physics simulation such as:
   1. rendering images (e.g. POV-Ray)
   2. particle simulations
2. geometry picking in an application
3. acoustic ray tracing
4. electromagnetic shooting-and-bouncing-rays (SBR)


## Adapter interface

While the canonical implementation of the Adapter design pattern uses polymorphism, Calico implements the pattern through C++ templates; an approach inspired by [nanoflann](https://github.com/jlblancoc/nanoflann) by Jose-Luis Blanco-Claraco.

The disadvantage of this template based implementation is that errors from the compiler are considerably more cryptic than when using polymorphism.
However, the performance should be substantially higher for two reasons:

1. there are no virtual functions and hence no vtable lookups/branch mispredictions
2. your adapter code is available to the compiler for inlining.

The second point means that a really smart compiler could produce an implementation that is just as fast as if I had designed Calico to use your datastructures directly rather than through an interface.
I.e. the adapter overhead is reduced at compile time to near zero.
However, the performance of the ray casting is dependent on how memory/cache friendly your adapted data structures are.
I recommend using a structure-of-arrays approach to designing your data structures as that generally provides a cache friendly data layout, and sets up conditions to best allow the compiler to auto-vectorize the implementation.
See Intel's [Optimization Reference Manual](http://www.intel.com/content/www/us/en/architecture-and-technology/64-ia-32-architectures-optimization-manual.html) for more best practices.

## Usage

As previously stated, Calico traverses your existing data structures by using adapters that you provide.
Details of the various adapters are provided in the following sections.
In its most basic usage, you must provide a Mesh adapter that allows Calico to walk your triangular mesh.
The rest of the system can utilize adapters provided within the Calico library itself, or you can provide your own structures for advanced control.

Calico can be used by defining a Tracer object and assigning to that Tracer a Mesh Adapter and Acceleration Structure.
Below is a simple example that traces a single ray against a compile-time defined square mesh object named Plate.

```C++
  #include <calico/tracer.hpp>
  #include <calico/math.hpp>
  #include <calico/input/soa_input.hpp>
  #include <calico/result/soa_result.hpp>
  #include <calico/accelerator/brute_force.hpp>

  #include <calico/utilities/meshes/plate.hpp>

  #include <iostream>

  void example() {
    // Choose the system internal precision
    typedef double Float;

    // Typedef the adapters using utilities that are provided with Calico.  Many
    // of the adapters are based on the Mesh adapter, so that's where we must
    // start.
    typedef calico::utilities::meshes::Plate<Float> Mesh;
    typedef calico::math::PluckerContainmentTest<Float, Mesh> Containment;
    typedef calico::accelerator::BruteForce<Float, Mesh, Containment> Accelerator;

    // The Plate mesh defines a constant unit square mesh object at compile
    // time. This could also load a mesh from disk, or provide a view into
    // your existing mesh data.
    Mesh plate(plate_string);
    Accelerator accelerator(plate);

    // We can finally create a tracer object. This object will traverse the
    // accelerator and mesh to determine what each ray strikes when we call
    // trace_rays. It can be used over and over again. 
    auto tracer = calico::make_tracer(accelerator);

    // Create one ray pointing towards the plate using the Structure-of-Arrays
    // (SoA) SoaInput adapter provided with Calico. We're using arrays, but these
    // could just as easily be std::vectors, they just need to be "array-like"
    // using array access syntax.
    const std::size_t ray_count{1};
    Mesh::FaceId start_face_id[1] = {Mesh::ray_miss_id_c};
    Float start_x[1]              = {0.};
    Float start_y[1]              = {0.};
    Float start_z[1]              = {10.};
    Float direction_x[1]          = {0.};
    Float direction_y[1]          = {0.};
    Float direction_z[1]          = {-1.};
    auto rays =
      calico::input::make_soa_input<Float, Mesh::FaceId>
          (ray_count, start_face_id, start_x, start_y, start_z,
           direction_x, direction_y, direction_z);

    // Make a SoaResult structure with enough storage space to hold the results
    // from tracing all of the rays (one ray, in our example). You can write your 
    // result structures, but SoaResult is provided with Calico.
    Mesh::FaceId strike_face_id[1] = {Mesh::ray_miss_id_c};
    Float t[1]                     = {-1000.};
    Float hit_x[1]                 = {-1000.};
    Float hit_y[1]                 = {-1000.};
    Float hit_z[1]                 = {-1000.};
    auto results =
      calico::result::make_soa_result<Float, std::size_t>
          (strike_face_id, t, hit_x, hit_y, hit_z);

    // Trace the rays. The strike object (if any) and its distance are stored in
    // results.
    tracer.trace_rays(rays, results);

    if (strike_id[0] != Mesh::ray_miss_id_c) {
      std::cerr << "Hit facet " << strike_id[0] << " at ("
                << hit_x[0] << ", " << hit_y[0] << ", " << hit_z[0]
                << ")" << std::endl;
    }
    else {
      std::cerr << "The ray missed!" << std::endl;
    }
  }
```

## Mesh

All usages of Calico will involve defining an adapter to your mesh data structures.
The adapter provides Calico with a way to traverse the mesh by investigating it one triangle at a time.
Internally, Calico views the world through Structure-of-Arrays style interfaces, but your data does not need to be laid out that way.
The interface is the only thing that is important to Calico.
Each mesh adapter must provide at least the following interface, but may provide additional convenience methods for your own use.

  `FloatType typedef`
   This typedef provides access to the floating-point type used by the mesh. A good choice is `float` or `double`.

  `FaceId typedef`
   This typedef defines a mechanism for naming each face. A good choice is an `std::int32_t` (aka int on most platforms).
   Some of the algorithms iterate over the face IDs using a ++ operator, so this type needs to be able to increment as an iterator.

  `VertexId typedef`
   This typedef defines a mechanism for naming each vertex. A good choice is an `std::int32_t` (aka int on most platforms).
   Some of the algorithms iterate over the vertex IDs using a ++ operator, so this type needs to be able to increment as an iterator.

  `static const FaceId ray_miss_id_c`
   This constant defines the ID that Calico can use to indicate that a ray did not strike any face when tracing the scene.
   It should be a value that will never be assigned as the identifier for a face in the mesh.
   A common choice is -1 when a signed data type is used.

  `bounding_box(Float &min_x, Float &min_y, Float &min_z,
                Float &max_x, Float &max_y, Float &max_z) const`
   The adapter should be able to compute a bounding box for the mesh and set the minimum and maximum values into the variables passed by reference into this routine.
   It should also not modify the underlying mesh in the process (hence its const attribute).
   Each of Float should be of the same floating point precision type that was chosen for Calico's internal use (see example above).

  `Float x(FaceId id, VertexId corner) const`
  `Float y(FaceId id, VertexId corner) const`
  `Float z(FaceId id, VertexId corner) const`
   This method of the adapter should return the x, y or z component of the `corner`'th vertex of triangle `id`.
   For example, suppose triangle 2 of the mesh is defined with the following three vertices: `(0, 0, 0) (6, 2, 3), (-4, 5, 1)`.
   Calling `x(2, 0)` would return 0; `x(2, 1)` would return 6; and `x(2, 2)` would return -4.
   This routine should return a floating point value.

  `Float normal_x(FaceId id) const`
  `Float normal_y(FaceId id) const`
  `Float normal_z(FaceId id) const`
   Returns the x, y or z component of face `id`'s unit normal.
   Because this routine will be called frequently and because we query the x, y and z components separately, it is a good idea to compute the unit normal a priori and store it for future accesses.

  `Float d(FaceId id) const`
   Returns the distance that the triangle lies along the normal vector from the origin.

  `Float area(FaceId id) const`
   Returns the surface area of the triangle in model units.

  `FaceId size() const`
   Number of triangles in the mesh.

# Plan

**Calico is under development and does not currently function, except in the most basic manner.**

- [X] Design mesh adapter interface
  - [ ] Explore refactoring the mesh adapter interface to use iterators for face IDs and vertex IDs.
- [X] Design ray adapter interface
  - [X] Provide helper adapter for structure-of-arrays ray inputs
- [X] Implement a math library to support ray-tracing routines
- [X] Implement ray/bounding-box intersection routines
  - [X] Provide unit tests verifying ray/bounding-box intersection routines work in a variety of cases
- [X] Write ray/triangle intersection routines in an abstract manner
  - [ ] Write ray/triangle intersection unit tests to verify routines work in a variety of cases
- [ ] Design acceleration structure interface
  - [X] Implement brute-force acceleration structure that checks each ray against all triangles; complexity is O(n)
    - [ ] Write unit tests to verify acceleration structure works in a variety of cases
  - [ ] Implement bounding-volume-hierarchy acceleration structure to speed ray/triangle checks; complexity is O(log(n))
    - [ ] Write unit tests to verify acceleration structure works in a variety of cases
- [ ] Write demonstration image rendering application using framework to render Wavefront OBJ meshes

# License

Calico is provided under the simplified BSD license.
See LICENSE for details.
