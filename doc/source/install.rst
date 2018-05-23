Install fastscapelib
====================

This library is header only and uses C++14 standards. It depends on
xtensor_. You can install xtensor, e.g., using conda_::

  $ conda install xtensor -c conda-forge

Installation from source using cmake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the commands below from the source directory to install the
fastscapelib's header files using cmake (on Unix platforms)::

  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ make install

Where ``/path/to/prefix`` is the path where the header files will be
installed.

On Windows platforms::

  $ mkdir build
  $ cd build
  $ cmake -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ nmake
  $ nmake install

.. _xtensor: https://github.com/QuantStack/xtensor
.. _conda: https://conda.io/docs/
