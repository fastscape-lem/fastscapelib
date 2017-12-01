Fastscape Core
==============

A C++ library of efficient algorithms for processing topographic data
and landscape evolution modeling.

Installation
------------

This library is header only and uses C++14 standards. It depends on
xtensor_, which can be installed using the conda package manager::

  $ conda install xtensor -c conda-forge

Currently, this repository also provides Python bindings, which
requires xtensor-python_ and pybind11_::

  $ conda install xtensor-python pybind11 -c conda-forge

There is no ``setup.py`` file yet to install the Python package, but
you can use cppimport_ to compile and import the python module in one
step. For this you can install cppimport using pip::

  $ pip install cppimport

Then in a Python console (launched from within this directory)::

  >>> import cppimport
  >>> fscape = cppimport.imp('pyfastscape')

.. _xtensor: https://github.com/QuantStack/xtensor
.. _xtensor-python: https://github.com/QuantStack/xtensor-python
.. _pybind11: https://github.com/pybind/pybind11
.. _cppimport: https://github.com/tbenthompson/cppimport
