(release-notes)=

# Release Notes

## v0.3.0 (Unreleased)

## v0.2.1 (20 October 2023)

### Bug fixes

- Fixed ``pflood_sink_resolver`` that was "flooding too much" ({issue}`145`,
  {pull}`146`).
- Fixed ``RuntimeWarning`` (invalid value in cast) issued when calling
  ``FlowGraph.basins()`` in Python ({pull}`147`).
- Python bindings: fixed Python interpreter crash (segfault) when accessing data
  members of the implementation of a graph snapshot via ``FlowGraph.impl()``
  ({pull}`148`).

## v0.2.0 (9 October 2023)

A complete re-write of Fastscapelib with brand new features and API (note: the
API of the previous version has been mostly removed). Some highlights:

- A flexible grid system, including 1D profile grid, 2D raster grid (with
  support of diagonals vs. non-diagonals connectivity) and 2D triangular mesh
- Full support of looped (periodic) boundary conditions for uniform rectangular
  grids
- Flexible flow routing using a "flow graph" and "flow operators"
- Support for both single direction and multiple direction flow (grid-agnostic
  implementation)
- Efficient resolution of closed depressions in the topography while routing the
  flow paths, based either on explicit computation of a graph of basins or on
  the priority flood algorithm (grid-agnostic implementation)
- The current flow graph implementation graph based on fixed-size arrays is
  extensible to alternative representations (e.g., sparse matrix, linked-lists)
- First-class C++ and Python APIs
- Detailed documentation (examples, user-guide, API reference, etc.)

Thanks to the contributors to this release: Benoît Bovy, Adrien Delsalle,
Guillaume Cordonnier and Hannah Sophia Davies. Thanks also to Johan Mabille and
QuantStack (https://quantstack.net/) for their contribution and advice on the
design of the library.

## v0.1.3 (5 November 2018)

- Update to xtensor latest version (0.18) ({issue}`35`).

## v0.1.2 (23 August 2018)

- More robust and cleaner version management ({issue}`32`).
- Support older versions of apple-clang ({issue}`30`, {issue}`33`).

## v0.1.1 (22 August 2018)

- Fix version ({issue}`29`).

## v0.1.0 (21 August 2018)

Initial release.
