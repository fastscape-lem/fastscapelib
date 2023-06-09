.. _install:

Install Fastscapelib
====================

Choose the way to install Fastscapelib that best suits your needs. In general we
recommend installing either the C++ or Python library using conda_ (or mamba_).

.. |conda-logo| image:: _static/conda_logo.svg
   :width: 20%
   :height: 2ex
   :class: no-scaled-link

.. grid:: 2
   :gutter: 3

   .. grid-item-card:: C++ / |conda-logo|
       :img-top: _static/cpp_logo.svg
       :link: install-cpp-conda
       :link-type: ref
       :text-align: center

       Install the C++ (header-only) library with conda_.

   .. grid-item-card:: C++ / from source
       :img-top: _static/cpp_logo.svg
       :link: install-cpp-source
       :link-type: ref
       :text-align: center

       Download and install the C++ (header-only) library with CMake.

   .. grid-item-card:: Python / |conda-logo|
       :img-top: _static/python_logo.svg
       :link: install-python-conda
       :link-type: ref
       :text-align: center

       Install the pre-compiled Python library with conda_.

   .. grid-item-card:: Python / from source
       :img-top: _static/python_logo.svg
       :link: install-python-source
       :link-type: ref
       :text-align: center

       Download, build and install the Python library from source using Pip.


.. _download-fastscapelib:

Download the Fastscapelib Source
--------------------------------

If you want to build Fastscapelib from source, you can either download the last
release `here <https://github.com/fastscape-lem/fastscapelib/releases/latest>`_
or you can clone the source repository using git (by default it will clone the
active development ``main`` branch):

.. code-block:: bash

   $ git clone https://github.com/fastscape-lem/fastscapelib

.. _install-cpp:

Install the C++ Library
-----------------------

Fastscapelib C++ is header only and requires a recent compiler that supports the
C++17 standards. It also depends on:

- xtensor_ (C++ array computing)
- cmake_ (build and configuration)

See Section :doc:`build_options` for more details on how to include Fastscapelib
in a (CMake) project.

.. _install-cpp-conda:

Using Conda
~~~~~~~~~~~

Fastscapelib conda_ packages are available via the conda-forge channel. You can
install Fastscapelib's C++ headers all dependencies using the following command
(alternatively you can use mamba_):

.. code-block:: bash

   $ conda install fastscapelib -c conda-forge

.. _install-cpp-source:

From Source Using CMake
~~~~~~~~~~~~~~~~~~~~~~~

In addition to a C++ compiler supporting C++17, you need xtensor_ and cmake_
that you can install, e.g., using conda_:

.. code-block:: bash

  $ conda install xtensor cmake -c conda-forge

After :ref:`downloading the Fastscapelib source <download-fastscapelib>`, run
the commands below from the source root directory to install the fastscapelib's
header files using CMake:

.. code-block:: bash

  $ cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/prefix ..
  $ cmake --build build
  $ cmake --install build

Where ``/path/to/prefix`` is the path where the header files will be installed
(skip this option if you want to install Fastscapelib in a default location).

See Section :doc:`build_options` for more information on the available build
options.

.. _install-python:

Install the Python Library
--------------------------

.. _install-python-conda:

Using Conda
~~~~~~~~~~~

Fastscapelib's Python bindings are available as binary conda_ packages for
Linux, MacOS and Windows via the conda-forge channel. You can install it using
the following command (alternatively you can use mamba_):

.. code-block:: bash

   $ conda install fastscapelib-python -c conda-forge

.. _install-python-source:

From Source Using Pip
~~~~~~~~~~~~~~~~~~~~~

Fastscapelib's Python bindings require Python (3.8+), numpy, pybind11_ and
xtensor-python_, which are all available on conda-forge:

.. code-block:: bash

  $ conda install python numpy pybind11 xtensor-python -c conda-forge

After :ref:`downloading the Fastscapelib source <download-fastscapelib>`, you
can build and install the Python package using ``pip``. Run the following
commands from the source root directory:

.. code-block:: bash

   $ cd python
   $ pip install .

.. _cmake: https://cmake.org/
.. _conda: https://conda.io/docs/
.. _mamba: https://mamba.readthedocs.io/en/latest/
.. _pybind11: https://github.com/pybind/pybind11
.. _xtensor: https://xtensor.readthedocs.io
.. _xtensor-python: https://xtensor-python.readthedocs.io
