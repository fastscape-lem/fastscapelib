.. _install:

Install Fastscapelib
====================

This library is header only and uses C++17 standards. It depends on
xtensor_ (>0.18).

.. _xtensor: https://github.com/QuantStack/xtensor

Install the C++ library
-----------------------

Using conda
~~~~~~~~~~~

Fastscapelib is available as a conda_ package in the conda-forge channel.

.. code-block:: bash

   $ conda install fastscapelib -c conda-forge

.. _conda: https://conda.io/docs/

From source using cmake
~~~~~~~~~~~~~~~~~~~~~~~

You need to first install the dependencies. You can install xtensor,
e.g., using conda_:

.. code-block:: bash

  $ conda install xtensor -c conda-forge

Run the commands below from the source directory to install the
fastscapelib's header files using cmake (on Unix platforms):

.. code-block:: bash

  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ make install

Where ``/path/to/prefix`` is the path where the header files will be
installed.

On Windows platforms:

.. code-block:: batch

  $ mkdir build
  $ cd build
  $ cmake -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ nmake
  $ nmake install

Install the Python library
--------------------------

Using conda
~~~~~~~~~~~

Fastscapelib's Python bindings is available as a separate conda_
package, still in the conda-forge channel.

.. code-block:: bash

   $ conda install fastscapelib-python -c conda-forge

From source using pip
~~~~~~~~~~~~~~~~~~~~~

Fastscapelib's Python bindings requires Python (2.7.x or 3.4+), numpy,
pybind11_ and xtensor-python_, which are also available through conda:

.. code-block:: bash

  $ conda install python numpy pybind11 xtensor-python -c conda-forge

The ``fastscapelib`` Python package can then be installed locally
using ``pip`` in editable mode:

.. code-block:: bash

   $ cd python
   $ pip install -e .

.. _pybind11: https://github.com/pybind/pybind11
.. _xtensor-python: https://github.com/QuantStack/xtensor-python
