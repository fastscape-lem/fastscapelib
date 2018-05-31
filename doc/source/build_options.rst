Build and configuration
=======================

Build options
-------------

``fastscapelib`` build supports the following options (all disabled by
default). See below for more explanations.

- ``BUILD_TESTS``: enables the ``run_test`` target.
- ``DOWNLOAD_GTEST``: downloads ``gtest`` and builds it locally
  instead of using a binary installation.
- ``GTEST_SRC_DIR``: indicates where to find the ``gtest`` sources
  instead of downloading them.
- ``BUILD_PYTHON_MODULE``: enables building fastscapelib as a Python
  extension.

Build and run tests
-------------------

Fastscapelib has a test suite based on GTest_. The enabled
``BUILD_TESTS`` adds the target ``run_tests`` which builds and runs
the whole test suite, e.g.,

.. code::

   $ mkdir build
   $ cd build
   $ cmake -DBUILD_TESTS=ON ..
   $ make run_tests

GTest must have been already installed, e.g., using conda:

.. code::

  $ conda install gtest -c conda-forge

Alternatively, GTest may be downloaded automatically by enabling
``DOWNLOAD_GTEST``, or a custom install path may be given by setting
``GTEST_SRC_DIR``. Note that enabling ``DOWNLOAD_GTEST`` or setting
``GTEST_SRC_DIR`` enables ``BUILD_TESTS``.

.. _GTest: https://github.com/google/googletest

Build the Python extension
--------------------------

Fastscapelib has Python bindings, which can be built by enabling
``BUILD_PYTHON_MODULE``, e.g.,

.. code::

   $ mkdir build
   $ cd build
   $ cmake -DBUILD_PYTHON_MODULE=ON ..
   $ make run_tests

Note that the created Python library is not intended to be directly
used within a regular Python installation. The preferred way to build
and install fastscapelib locally is to use pip:

.. code::

   $ cd python
   $ pip install -e .
