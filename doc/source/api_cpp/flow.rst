Flow Routing
============

Flow graph
----------

Defined in ``fastscapelib/flow/flow_graph.hpp``

.. doxygenclass:: fastscapelib::flow_graph
   :members:
   :undoc-members:

Implementations
~~~~~~~~~~~~~~~

Defined in ``fastscapelib/flow/flow_graph_impl.hpp``

.. doxygenstruct:: fastscapelib::flow_graph_fixed_array_tag

.. _cpp-api-flow-operators:

|

Flow operators
--------------

Defined in ``fastscapelib/flow/flow_operator.hpp``

.. doxygenenum:: fastscapelib::flow_direction

.. doxygenclass:: fastscapelib::flow_operator
   :members:

.. doxygenclass:: fastscapelib::flow_operator_sequence
   :members:

|

Flow routers
~~~~~~~~~~~~

Defined in ``fastscapelib/flow/flow_router.hpp``

- :cpp:class:`~fastscapelib::single_flow_router`
- :cpp:class:`~fastscapelib::multi_flow_router`

|

.. doxygenclass:: fastscapelib::single_flow_router
   :members:
   :undoc-members:

.. doxygenclass:: fastscapelib::multi_flow_router
   :members:
   :undoc-members:

|

Sink resolvers
~~~~~~~~~~~~~~

Defined in ``fastscapelib/flow/sink_resolver.hpp``

- :cpp:struct:`~fastscapelib::pflood_sink_resolver`
- :cpp:class:`~fastscapelib::mst_sink_resolver`

|

.. doxygenstruct:: fastscapelib::pflood_sink_resolver
   :members:
   :undoc-members:

.. doxygenenum:: fastscapelib::mst_route_method

.. doxygenclass:: fastscapelib::mst_sink_resolver
   :members:
   :undoc-members:

|

Flow snapshots
~~~~~~~~~~~~~~

Defined in ``fastscapelib/flow/flow_snapshot.hpp``

- :cpp:class:`~fastscapelib::flow_snapshot`

|

.. doxygenclass:: fastscapelib::flow_snapshot
   :members:
   :undoc-members:

|

Basin graph
-----------

Defined in ``fastscapelib/flow/basin_graph.hpp``

.. doxygenenum:: fastscapelib::mst_method

.. doxygenclass:: fastscapelib::basin_graph
   :members:
   :undoc-members:
