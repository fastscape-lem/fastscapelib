.. _build_options:

Build and Configuration
=======================

Build options
-------------

``fastscapelib`` build supports the following options (all disabled by
default). See below for more explanations.

- ``BUILD_TESTS``: enables the ``run_test`` target.
- ``BUILD_BENCHMARK``: enables the ``run_benchmark`` target.
- ``DOWNLOAD_GTEST``: downloads google-test and builds it locally
  instead of using a binary installation.
- ``GTEST_SRC_DIR``: indicates where to find the google-test sources
  instead of downloading them.
- ``BUILD_PYTHON_MODULE``: enables building fastscapelib as a Python
  extension.
- ``DOWNLOAD_XTENSOR``: downloads xtensor development version (master
  branch on github) and uses it to build fastscapelib (useful for
  testing - might be needed for building fastscapelib development
  version).

Build and run tests
-------------------

Fastscapelib has a test suite based on google-test_. The enabled
``BUILD_TESTS`` adds the target ``run_tests`` which builds and runs
the whole test suite, e.g.,

.. code::

   $ mkdir build
   $ cd build
   $ cmake -DBUILD_TESTS=ON ..
   $ make run_tests

Google-test must have been already installed, e.g., using conda:

.. code::

  $ conda install gtest -c conda-forge

Alternatively, google-test may be downloaded automatically by enabling
``DOWNLOAD_GTEST``, or a custom install path may be given by setting
``GTEST_SRC_DIR``. Note that enabling ``DOWNLOAD_GTEST`` or setting
``GTEST_SRC_DIR`` enables ``BUILD_TESTS``.

.. _google-test: https://github.com/google/googletest

Build and run benchmarks
------------------------

Fastscapelib has also a benchmark suite based on
google-benchmark_. Building and running benchmarks is similar to
building and running tests (the ``BUILD_BENCHMARK`` option and
``run_benchmark`` target are used instead). Note that
``BUILD_BENCHMARK`` automatically downloads google-benchmark.

.. _google-benchmark: https://github.com/google/benchmark

Build the Python extension
--------------------------

Fastscapelib has Python bindings, which can be built by enabling
``BUILD_PYTHON_MODULE``, e.g.,

.. code::

   $ mkdir build
   $ cd build
   $ cmake -DBUILD_PYTHON_MODULE=ON ..
   $ make _fastscapelib_py

Note that the created Python library is not intended to be directly
used within a regular Python installation. The preferred way to build
and install fastscapelib locally is to use pip:

.. code::

   $ cd python
   $ pip install -e .
