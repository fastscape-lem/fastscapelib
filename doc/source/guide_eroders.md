# Eroders

Fastscapelib implements efficient numerical solutions of few common erosion
equation terms used in Landscape Evolution Models (LEMs).

:::{note}

This is currently limited to bedrock channel erosion (Stream Power Law) on any grid and
hillslope erosion / transport (linear diffusion) on raster grids only.

We plan to add more processes in the future, like the continental and marine
sediment transport / deposition modules available in the [Fastscapelib
(Fortran)](https://fastscape.org/fastscapelib-fortran/) library, glacial
erosion, etc.

:::

## Eroder Class Interface

Each erosion process is exposed in Fastscapelib via an "Eroder" class with the
following conventions:

- the class constructor may take several arguments including a grid object, a
  graph object and/or values for the eroder input parameters.

- the eroder parameters are exposed as properties or getter (setter) methods so
  it is possible to read (update) their values after the eroder instance has
  been created (e.g., external forcing).

- the class exposes an ``erode()`` method that usually takes an input
  topographic elevation and a time step duration as arguments (+ maybe other
  arguments) and that return the amount of erosion (positive values) or
  accumulation (negative values) for one time step computed at
  each grid node.

See the {ref}`C++ <cpp-api-eroders>` and {ref}`Python <py-api-eroders>` API
reference for more details on each available erosion process.

## Examples

See {doc}`examples/index`.
