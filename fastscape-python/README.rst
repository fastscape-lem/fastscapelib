Fastscape (Python)
==================

A collection of efficient algorithms for processing topographic data
and landscape evolution modeling, available from within Python.

Installation
------------

Building the Python bindings requires numpy, pybind11_ and xtensor-python_.

.. _pybind11: https://github.com/pybind/pybind11
.. _xtensor-python: https://github.com/QuantStack/xtensor-python

Installation from source in a conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The procedure below assumes that you have first installed the
Fastscape C++ library in a conda envrionment by following the steps in
the main README file of this repository.

To install the required dependencies in the activated conda
environment::

  $ conda install python=3.6 numpy pybind11 xtensor-python -c conda-forge

To install the Python package, from the source directory::

  $ pip install .

Import the package from within Python
-------------------------------------

In a Python console::

  >>> import fastscape
