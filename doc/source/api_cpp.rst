.. _api_cpp:

C++ API Reference
=================

Fastscapelib heavily relies on xtensor_ for handling n-dimensional
arrays. Almost all arrays defined in the C++ interface have the
``xtensor_t`` type, which refers to any xtensor-compatible array
such as :cpp:type:`xt::xtensor <xt::xtensor>`, :cpp:type:`xt::xarray
<xt::xtensor>`, :cpp:type:`xt::pytensor <xt::pytensor>` or
:cpp:type:`xt::pyarray <xt::pyarray>` (this list is extensible).

In fact, ``xtensor_t`` is just an alias of xtensor's class
:cpp:type:`xt::xexpression <xt::xexpression>`. It has a unique
template parameter which refers to the derived array type like one of
those above. By convention, the parameter name is chosen using the
first character(s) of the corresponding argument in function
signatures.

.. _xtensor: https://xtensor.readthedocs.io/en/latest/

The fastscapelib library is organized into different topics listed
here below. The namespace ``fastscapelib`` is used for the public
API. You can either include all functions in your project at once:

.. code-block:: cpp

   #include "fastscapelib/fastscapelib.hpp"

or include just a subset of the functions (by topic), e.g.,

.. code-block:: cpp

   #include "fastscapelib/flow_routing.hpp"

Sinks (depressions)
-------------------

Functions used for depression filling or pit resolving.

Defined in ``fastscapelib/sinks.hpp``.

Sink filling
~~~~~~~~~~~~

.. doxygenfunction:: fastscapelib::fill_sinks_flat
   :project: fastscapelib

.. doxygenfunction:: fastscapelib::fill_sinks_sloped
   :project: fastscapelib

Flow routing
------------

Functions used to route (water) flow on a topographic surface and
compute flow path-related features or structures.

Defined in ``fastscapelib/flow_routing.hpp``.

Flow topology
~~~~~~~~~~~~~

.. doxygenfunction:: fastscapelib::compute_receivers_d8
   :project: fastscapelib

.. doxygenfunction:: fastscapelib::compute_donors
   :project: fastscapelib

Flow tree sorting
~~~~~~~~~~~~~~~~~

.. doxygenfunction:: fastscapelib::compute_stack
   :project: fastscapelib

Drainage area, basins, outlets & pits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: fastscapelib::compute_basins
   :project: fastscapelib

.. doxygenfunction:: fastscapelib::find_pits
   :project: fastscapelib

.. doxygenfunction:: fastscapelib::compute_drainage_area(xtensor_t<D>&, const xtensor_t<C>&, const xtensor_t<S>&, const xtensor_t<R>&)
   :project: fastscapelib

.. doxygenfunction:: fastscapelib::compute_drainage_area(xtensor_t<D>&, const xtensor_t<S>&, const xtensor_t<R>&, double, double)
   :project: fastscapelib

Bedrock channel erosion
-----------------------

.. doxygenfunction:: fastscapelib::erode_stream_power
   :project: fastscapelib
