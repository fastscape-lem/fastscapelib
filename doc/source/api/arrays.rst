N-d dimensional arrays
======================

Fastscapelib heavily relies on xtensor_ for handling n-dimensional
arrays. All arrays defined in the C++ interface have the
``xtensor_t<D>`` type, which refers to any xtensor-compatible array
such as :cpp:type:`xt::xtensor <xt::xtensor>`, :cpp:type:`xt::xarray
<xt::xtensor>`, :cpp:type:`xt::pytensor <xt::pytensor>` or
:cpp:type:`xt::pyarray <xt::pyarray>` (this list is extensible).

In fact, ``xtensor_t<D>`` is just an alias of the base class
``xt::xcontainer<D>``. The unique template parameter (here noted
``D``) refers to the derived type like one of those above. By
convention, the parameter name is chosen using the first character(s)
of the corresponding argument in function signatures.

The lazy expression capabilities of xtensor (i.e.,
:cpp:type:`xt::xexpression <xt::xexpression>`) are not used here.

.. _xtensor: https://xtensor.readthedocs.io/en/latest/
