Examples
========

A gallery of a few high-level examples in Python or C++ showcasing
Fastscapelib's features (multiple grids, flexible boundary conditions, flexible
flow routing, eroder classes, etc.)

|

.. |python-logo| image:: ../_static/python_logo.svg
   :width: 7%
   :height: 2ex
   :class: no-scaled-link

.. |cpp-logo| image:: ../_static/cpp_logo.svg
   :width: 7%
   :height: 2ex
   :class: no-scaled-link

.. grid:: 2
   :gutter: 3

   .. grid-item-card:: River Profile
       :link: river_profile_py
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_river_profile.png

       |python-logo|
       ^^^

       1-d profile grid, block-uplift + stream-power law.

   .. grid-item-card:: River Profile
       :link: river_profile_cpp
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_river_profile.png

       |cpp-logo|
       ^^^

       1-d profile grid, block-uplift + stream-power law.

   .. grid-item-card:: Mountain
       :link: mountain_py
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_mountain.png

       |python-logo|
       ^^^

       2-d raster grid, block uplift + stream-power law + hillslope diffusion.

   .. grid-item-card:: Mountain
       :link: mountain_cpp
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_mountain.png

       |cpp-logo|
       ^^^

       2-d raster grid, block uplift + stream-power law + hillslope diffusion.

   .. grid-item-card:: Escarpment
       :link: escarpment_py
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_escarpment.png

       |python-logo|
       ^^^

       2-d raster grid, semi-infinite domain (half-periodic boundaries),
       multiple direction flow, stream-power law + hillslope diffusion.

   .. grid-item-card:: Inner Base Levels
       :link: inner_base_levels_py
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_inner_base_levels.png

       |python-logo|
       ^^^

       2-d raster grid, infinite domain (full-periodic boundaries), custom base
       level nodes, uplift + stream-power law + hillslope diffusion.

   .. grid-item-card:: Catchment (TriMesh)
       :link: catchment_py
       :link-type: doc
       :text-align: center
       :img-bottom: ../_static/fig_catchment.png

       |python-logo|
       ^^^

       2-d triangular (irregular) mesh, one catchment with one base level node
       (outlet), stream-power law + base level lowering.

.. toctree::
   :hidden:
   :maxdepth: 1

   river_profile_py
   river_profile_cpp
   mountain_py
   mountain_cpp
   escarpment_py
   inner_base_levels_py
   catchment_py
