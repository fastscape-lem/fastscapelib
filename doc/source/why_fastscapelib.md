# Why Fastscapelib?

There are already many Landscape Evolution Models (LEMs) available as well as a
few great open-source frameworks for landscape evolution modeling such as
[Landlab](https://landlab.readthedocs.io) and
[Badlands](https://badlands.readthedocs.io/), so why do we need yet another
tool?

Fastscapelib's philosophy can be summarized as *"write algorithms once, reuse it
everywhere"*.

The last decades have seen a great number of LEM alternatives proposed (see,
e.g., {cite:t}`Barnhart2019`), yet most alternatives depend on a few, common
concepts and algorithms, e.g., for routing flow on the topographic surface. It
is quite usual to see the same -- sometimes complex -- algorithms re-implemented
over and over, with only minor differences from one implementation to another.

## Flexibility, Reusability and Efficiency

Fastscapelib's goal is to provide a core, "bare-bones" library that is focused
on the three following aspects: *flexibility, reusability and efficiency*.

Maximizing reusability is made via a flexible grid system that is decoupled from
flow routing and through first-class APIs in both C++ and Python. Support could
be easily extended to other languages like R or Julia thanks to the
[Xtensor](https://xtensor.readthedocs.io) library (more details {ref}`here
<internals-xtensor-language-bindings>`). In addition, Fastscapelib makes as few
as possible assumptions about how a LEM should be organized (e.g., variable
names, units, etc.).

Optimizing for both flexibility and efficiency is hard (there are always
trade-offs) but we did our best trying to achieve it at the expense of a fairly
complex {doc}`internal C++ implementation <internals>` although with extensive
unit tests and carefully designed API and abstractions.

## What is Fastscapelib for?

Fastscapelib provides some core structures (grids, meshes, flow graphs) and
solvers (diffusion, Stream Power Law) for building custom LEMs and/or for
integrating it with other existing codes, models (e.g., geodynamics, climate) or
frameworks written in various languages.

Fastscapelib will eventually replace the current [Fortran
backend](https://fastscape.org/fastscapelib-fortran/) of the
[Fastscape](https://fastscape.readthedocs.io) high-level Python package.

## What is Fastscapelib *not* for?

Fastscapelib is a bare-bones library. It doesn't provide any ready-to-use LEM
nor anything for setting up simulations, driving it and processing or
visualizing the outputs. Instead, you can freely use your tools of choice for
those tasks.

Fastscapelib will never provide a list of model components that is as exhaustive
as in [Landlab](https://landlab.readthedocs.io).

## Roadmap

Fastscapelib is still at an early stage of development and only provides a
limited number of features regarding flow routing, bedrock erosion, sediment
transport, etc.

TODO.
