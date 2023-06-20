# Fastscapelib Documentation

Fastscapelib is a C++/Python library of efficient and reusable algorithms for
landscape evolution modeling. It aims at providing a (bare bones) toolkit for
building your own Landscape Evolution Models (LEMs) or integrating it with other
existing codes or models.

**Useful links**:
[Home](http://fastscapelib.readthedocs.io/) |
[Fastscape Home](http://fastscape.org) |
[Code Repository](https://github.com/fastscape-lem/fastscapelib) |
[Issues](https://github.com/fastscape-lem/fastscapelib/issues) |
[Discussions](https://github.com/fastscape-lem/fastscapelib/discussions) |
[Releases](https://github.com/fastscape-lem/fastscapelib/releases)

::::{grid} 2
:gutter: 3

:::{grid-item-card} Getting Started
:img-top: _static/index_getting_started.svg

New to Fastscapelib? Check the guide on how to install the library and get
started with a few examples.

- {doc}`why_fastscapelib`
- {doc}`install`
- {doc}`examples/index`
:::
:::{grid-item-card} User Guide
:img-top: _static/index_user_guide.svg

The user guide provides in-depth information on the key concepts used in
Fastscapelib.

- {doc}`guide_grids`
- {doc}`guide_flow`
- {doc}`guide_eroders`
:::
:::{grid-item-card} API Reference
:img-top: _static/index_api.svg

The reference guide describes in detail all the functions and classes that are
part of Fastscapelib's public API.

- {doc}`api_cpp/index`
- {doc}`api_python/index`
:::
:::{grid-item-card} Developer Guide
:img-top: _static/index_contribute.svg

The developer guide is for anyone who would like to customize, reuse or
contribute improving the Fastscapelib library!

- {doc}`build_options`
- {doc}`internals`
- {doc}`contributing`
:::
::::

```{toctree}
:caption: Getting Started
:hidden: true
:maxdepth: 1

why_fastscapelib
install
examples/index
references
release_notes
```

```{toctree}
:caption: User Guide
:hidden: true
:maxdepth: 1

guide_grids
guide_flow
guide_eroders
```

```{toctree}
:caption: API Reference
:hidden: true
:maxdepth: 1

api_cpp/index
api_python/index
```

```{toctree}
:caption: Developer Guide
:hidden: true
:maxdepth: 1

build_options
internals
contributing
```

## Citing Fastscapelib

{ref}`How to cite Fastscapelib? <citation>`

## Acknowledgment

This project is supported by the [Earth Surface Process Modelling]
group of the GFZ Helmholtz Centre Potsdam.

[earth surface process modelling]: http://www.gfz-potsdam.de/en/section/earth-surface-process-modelling/
