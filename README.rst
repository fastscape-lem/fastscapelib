Fastscapelib
============

|Build Status| |Doc Status|

A C++ library of efficient algorithms for processing topographic data
and landscape evolution modeling.

This library has also Python bindings: see the README file in the
``fastscapelib-python`` subfolder of this repository for more information
on how to install and use it.

.. |Build Status| image:: https://travis-ci.org/fastscape-lem/fastscapelib.svg?branch=master
   :target: https://travis-ci.org/fastscape-lem/fastscapelib
   :alt: Build Status
.. |Doc Status| image:: http://readthedocs.org/projects/fastscapelib/badge/?version=latest
   :target: http://fastscapelib.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Installation
------------

This library is header only and uses C++14 standards. It depends on
xtensor_.

Installation from source in a conda environment, using cmake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend installing Fastscapelib in its own conda_ environment. To
create a new environment named 'fastscape' and install the required
dependencies::

  $ conda create -n fastscape xtensor cmake -c conda-forge

To activate the environment (on Unix platforms)::

  $ source activate fastscape

To activate the environment (on Windows platforms)::

  $ activate fastscape

Run the commands below from the source directory to install the
fastscapelib's header files using cmake (on Unix platforms)::

  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ make install

Where ``/path/to/prefix`` is the path to the newly created conda environment.

On Windows platforms::

  $ mkdir build
  $ cd build
  $ cmake -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ nmake
  $ nmake install

.. _xtensor: https://github.com/QuantStack/xtensor
.. _conda: https://conda.io/docs/

Testing
-------

Fastscapelib has a test suite based on GTest_. To install it you can
use conda too::

  $ conda install gtest -c conda-forge

To build and run the test suite (Unix platforms), from the build
directory created above::

  $ cmake -DBUILD_TESTS=ON ..
  $ make run_tests

.. _GTest: https://github.com/google/googletest

Documentation
-------------

Documentation is hosted on ReadTheDocs:
http://fastscapelib.readthedocs.io

Acknowledgment
--------------

This project is supported by the `Earth Surface Process Modelling`_
group of the GFZ Helmholtz Centre Potsdam.

.. _`Earth Surface Process Modelling`: http://www.gfz-potsdam.de/en/section/earth-surface-process-modelling/
