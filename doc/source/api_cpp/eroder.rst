.. _cpp-api-eroders:

Eroders
=======

"Eroder" classes implement various erosion processes. They all provide an
``erode`` method for computing erosion and/or updating the elevation of
the topographic surface during one time step.

Channel erosion
---------------

Defined in ``fastscapelib/eroders/spl.hpp``.

.. doxygenclass:: fastscapelib::spl_eroder
   :members:
   :undoc-members:

.. doxygenfunction:: fastscapelib::make_spl_eroder

Hillslope erosion
-----------------

Defined in ``fastscapelib/eroders/diffusion_adi.hpp``.

.. doxygenclass:: fastscapelib::diffusion_adi_eroder
   :members:
   :undoc-members:

.. doxygenfunction:: fastscapelib::make_diffusion_adi_eroder
