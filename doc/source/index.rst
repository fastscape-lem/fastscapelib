Fastscapelib Documentation
==========================

Fastscapelib is a C++/Python library of efficient and reusable algorithms for
processing topographic data and landscape evolution modeling. It aims at providing a
toolkit for building your own Landscape Evolution Models (LEMs).

**Useful links**:
`Home <http://fastscapelib.readthedocs.io/>`__ |
`Fastscape Home <http://fastscape.org>`__ |
`Code Repository <https://github.com/fastscape-lem/fastscapelib>`__ |
`Issues <https://github.com/fastscape-lem/fastscapelib/issues>`__ |
`Discussions <https://github.com/fastscape-lem/fastscapelib/discussions>`__ |
`Releases <https://github.com/fastscape-lem/fastscapelib/releases>`__

|

.. grid:: 2
   :gutter: 3

   .. grid-item-card:: Getting Started
      :img-top: _static/index_getting_started.svg

      New to Fastscapelib? Check the guide on how to install the library and get
      started with a few examples.

      - :doc:`install`
      - :doc:`examples/index`

   .. grid-item-card::  User Guide
      :img-top: _static/index_user_guide.svg

      The user guide provides in-depth information on the key concepts used in
      Fastscapelib.

      - :doc:`guide_grids`
      - :doc:`guide_flow`
      - :doc:`guide_eroders`

   .. grid-item-card:: API Reference
      :img-top: _static/index_api.svg

      The reference guide describes in detail all the functions and classes that
      are part of Fastscapelib's public API.

      - :doc:`api_cpp/index`
      - :doc:`api_python/index`

   .. grid-item-card:: Developer Guide
      :img-top: _static/index_contribute.svg

      The developer guide is for anyone who would like to customize, reuse or
      contribute improving the Fastscapelib library!

      - :doc:`build_options`
      - :doc:`internals`


.. toctree::
   :caption: Getting Started
   :hidden:
   :maxdepth: 1

   install
   examples/index
   references
   release_notes

.. toctree::
   :caption: User Guide
   :hidden:
   :maxdepth: 1

   guide_grids
   guide_flow
   guide_eroders

.. toctree::
   :caption: API Reference
   :hidden:
   :maxdepth: 1

   api_cpp/index
   api_python/index

.. toctree::
   :caption: Developer Guide
   :hidden:
   :maxdepth: 1

   build_options
   internals

Citing Fastscapelib
-------------------

:ref:`How to cite Fastscapelib? <citation>`

Acknowledgment
--------------

This project is supported by the `Earth Surface Process Modelling`_
group of the GFZ Helmholtz Centre Potsdam.

.. _`Earth Surface Process Modelling`: http://www.gfz-potsdam.de/en/section/earth-surface-process-modelling/
