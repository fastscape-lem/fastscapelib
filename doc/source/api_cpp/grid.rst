Grids
=====

Grid elements
-------------

- :cpp:enum:`~fastscapelib::node_status`
- :cpp:struct:`~fastscapelib::node`
- :cpp:struct:`~fastscapelib::raster_node`
- :cpp:struct:`~fastscapelib::neighbor`
- :cpp:struct:`~fastscapelib::raster_neighbor`

|

Grid boundaries
---------------

- :cpp:class:`~fastscapelib::boundary_status`
- :cpp:class:`~fastscapelib::profile_boundary_status`
- :cpp:class:`~fastscapelib::raster_boundary_status`

|

Grid inner types
----------------

.. doxygenstruct:: fastscapelib::grid_inner_types

|

Base classes and structures
---------------------------

Defined in ``fastscapelib/grid/base.hpp``.

- :cpp:enum:`~fastscapelib::node_status`
- :cpp:struct:`~fastscapelib::node`
- :cpp:struct:`~fastscapelib::neighbor`
- :cpp:class:`~fastscapelib::boundary_status`
- :cpp:class:`~fastscapelib::grid`
- :cpp:class:`~fastscapelib::structured_grid`

|

.. doxygenenum:: fastscapelib::node_status

.. doxygenstruct:: fastscapelib::node
   :members:

.. doxygenstruct:: fastscapelib::neighbor
   :members:

.. doxygenclass:: fastscapelib::boundary_status
   :members:

.. doxygenclass:: fastscapelib::grid
   :members:
   :undoc-members:

Defined in ``fastscapelib/grid/structured_grid.hpp``

.. doxygenclass:: fastscapelib::structured_grid
   :members:
   :undoc-members:

|

Caching neighbor indices
------------------------

Defined in ``fastscapelib/grid/base.hpp``.

A very common repetitive task in Fastscapelib is getting the node indices of all
the node neighbors at a given grid node. Using a cache may speed-up this task
(at the cost of memory usage), especially for structured grids with looped
boundaries. In other cases using a cache won't provide any benefit, like for
unstructured meshes where the topology is already fully contained.

- :cpp:class:`~fastscapelib::neighbors_cache`
- :cpp:class:`~fastscapelib::neighbors_no_cache`

|

.. doxygenclass:: fastscapelib::neighbors_cache
   :members:
   :undoc-members:

.. doxygenclass:: fastscapelib::neighbors_no_cache
   :members:
   :undoc-members:

|

Profile grid
------------

Defined in ``fastscapelib/grid/profile_grid.hpp``.

- :cpp:class:`~fastscapelib::profile_boundary_status`
- :cpp:class:`~template\<class S, class C> fastscapelib::grid_inner_types\<profile_grid_xt\<S, C>>`,
- :cpp:class:`~template\<class S, class C = neighbors_cache\<2>> fastscapelib::profile_grid_xt`
- :cpp:type:`~fastscapelib::profile_grid`

|

.. doxygenclass:: fastscapelib::profile_boundary_status
   :members:

.. doxygenstruct:: fastscapelib::grid_inner_types< profile_grid_xt< S, C > >
   :members:
   :undoc-members:

.. doxygenclass:: fastscapelib::profile_grid_xt
   :members:
   :undoc-members:

.. doxygentypedef:: fastscapelib::profile_grid

|

Raster grid
-----------

Defined in ``fastscapelib/grid/raster_grid.hpp``.

- :cpp:struct:`~fastscapelib::raster_node`
- :cpp:struct:`~fastscapelib::raster_neighbor`
- :cpp:class:`~fastscapelib::raster_boundary_status`
- :cpp:enum:`~fastscapelib::raster_connect`
- :cpp:class:`~template\<class S, raster_connect RC, class C> fastscapelib::grid_inner_types\<raster_grid_xt\<S, RC, C>>`
- :cpp:class:`~template\<class S, raster_connect RC, class C = neighbors_cache\<raster_neighbors\<RC>::_n_neighbors_max>> fastscapelib::raster_grid_xt`
- :cpp:type:`~fastscapelib::raster_grid`

|

.. doxygenstruct:: fastscapelib::raster_node
   :members:

.. doxygenstruct:: fastscapelib::raster_neighbor
   :members:

.. doxygenclass:: fastscapelib::raster_boundary_status
   :members:

.. doxygenenum:: fastscapelib::raster_connect

.. doxygenstruct:: fastscapelib::grid_inner_types< raster_grid_xt< S, RC, C > >
   :members:
   :undoc-members:

.. doxygenclass:: fastscapelib::raster_grid_xt
   :members:
   :undoc-members:

.. doxygentypedef:: fastscapelib::raster_grid

|

Triangular mesh
---------------

Defined in ``fastscapelib/grid/trimesh.hpp``.

- :cpp:class:`~template\<class S, unsigned int N> fastscapelib::grid_inner_types\<trimesh_xt\<S, N>>`
- :cpp:class:`~template\<class S, unsigned int N = 30> fastscapelib::trimesh_xt`
- :cpp:type:`~fastscapelib::trimesh`

|

.. doxygenstruct:: fastscapelib::grid_inner_types< trimesh_xt< S, N > >
   :members:
   :undoc-members:

.. doxygenclass:: fastscapelib::trimesh_xt
   :members:
   :undoc-members:

.. doxygentypedef:: fastscapelib::trimesh
