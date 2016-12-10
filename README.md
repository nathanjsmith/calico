------
Calico
------

**Calico is under development and does not currently function.**

Calico is a ray/geometry intersection library for building ray-traced
applications.  Its focus is on presenting an API that works with your existing
datastructures through use of the [Adapter design
pattern](http://sourcemaking.com/design_patterns/adapter).

While some may call Calico a ray-trace engine, it doesn't quite provide
everything needed to make it a ray-tracer.  Instead, it focuses on providing
the tools necessary to _implement_ a ray-tracer.  For example, it focuses on
answering the question of "what surfaces did this ray hit?".  It does not provide
reflection or transmission calculations, lighting calculations, etc.  The
physics simulation is left to the implementer as they are the pieces that are
application specific.


Uses
====

Ideas for how Calico can be used include:

1. the basis of a ray-physics simulation such as:
  1. rendering images (e.g. POV-Ray)
  2. particle simulations
2. geometry picking in an application


Adapter interface
=================

While the canonical implementation of the Adapter design pattern uses
polymorphism, Calico implements the pattern through C++ templates; an approach
inspired by [nanoflann](https://github.com/jlblancoc/nanoflann) by Jose-Luis
Blanco-Claraco.

The disadvantage of this template based implementation is that errors from the
compiler are considerably more cryptic than when using polymorphism.  However,
the performance should be substantially higher for two reasons:

1. there are no virtual functions and hence no vtable lookups/branch
   mispredictions
2. your adapter code is available to the compiler for inlining.

The second point means that a really smart compiler could produce an
implementation that is just as fast as if I had designed Calico to use your
datastructures directly rather than through an interface.  I.e. the adapter
overhead is reduced at compile time to near zero.

License
=======

Calico is provided under the simplified BSD license.  See LICENSE for details.
