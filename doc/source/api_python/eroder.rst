.. _py-api-eroders:

Eroders
=======

"Eroder" classes implement various erosion processes. They all provide an
``erode`` method for computing erosion and/or updating the elevation of
the topographic surface during one time step.

.. currentmodule:: fastscapelib

Channel erosion
---------------

.. autosummary::
   :toctree: _api_generated/
   :template: fastscapelib-class-template.rst

   SPLEroder

Hillslope erosion
-----------------

.. autosummary::
   :toctree: _api_generated/
   :template: fastscapelib-class-template.rst

   DiffusionADIEroder
