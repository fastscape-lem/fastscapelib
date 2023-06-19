.. _build_options:

Build and Configuration
=======================

Include Fastscapelib in a CMake Project
---------------------------------------

After :ref:`installing Fastscapelib C++ headers <install-cpp>` (in a default
location), you should be able to use CMake's `find_package
<https://cmake.org/cmake/help/latest/command/find_package.html>`_ to ensure
these will be found when building your CMake project:

.. code:: cmake

   find_package(fastscapelib REQUIRED)

Don't forget to link your target library or application with Fastscapelib:

.. code:: cmake

   target_link_libraries(_your_target INTERFACE fastscapelib)

Build Options
-------------

Fastscapelib provides the following CMake build options (all disabled by
default). See below for more explanations.

.. list-table::
   :widths: 25 75

   * - ``BUILD_TESTS``
     - Enables :ref:`building the C++ tests <run-cpp-tests>`
   * - ``BUILD_BENCHMARK``
     - Enables :ref:`building the micro-benchmarks <run-benchmarks>`
   * - ``DOWNLOAD_GTEST``
     - Downloads google-test and builds it locally instead of using a binary
       installation
   * - ``GTEST_SRC_DIR``
     - Indicates where to find the google-test sources instead of downloading
       them
   * - ``BUILD_PYTHON_MODULE``
     - Enables building fastscapelib as a Python extension (internal use only!
       see :ref:`install-python` for instructions on how to build and install the
       Python library from source)
   * - ``DOWNLOAD_XTENSOR``
     - Downloads xtensor development version (master branch on GitHub) and uses
       it to build fastscapelib (useful for testing)

.. _run-cpp-tests:

Build and Run the C++ Tests
---------------------------

Fastscapelib has a test suite based on google-test_.

You can install google-test, e.g., using conda:

.. code::

  $ conda install gtest -c conda-forge

Alternatively, google-test may be downloaded automatically by enabling
``DOWNLOAD_GTEST``, or a custom install path may be given by setting
``GTEST_SRC_DIR`` (setting ``DOWNLOAD_GTEST=ON`` or
``GTEST_SRC_DIR=/path/to/gtest`` automatically sets ``BUILD_TESTS=ON``).

To build the tests, run the following commands from the source root directory:

.. code::

   $ cmake -S . -B build/tests -DBUILD_TESTS=ON
   $ cmake --build build/tests

Then to run all the tests:

.. code::

   $ ctest -T test --output-on-failure build/tests

.. _google-test: https://github.com/google/googletest

.. _run-python-tests:

Run the Python Tests
--------------------

Running the Python tests requires pytest_. You can install it using, e.g., conda:

.. code::

  $ conda install pytest -c conda-forge


After :ref:`(re)installing the Fastscapelib Python library <install-python>`,
you can run the tests using the following command (from the repository root
directory):

.. code::

   $ pytest -v .

.. _pytest: https://docs.pytest.org/

.. _run-benchmarks:

Build and Run the Benchmarks
----------------------------

Fastscapelib has also a micro-benchmark suite based on
google-benchmark_.

You can install google-benchmark, e.g., using conda:

.. code::

  $ conda install benchmark -c conda-forge

To build the benchmarks, run the following commands from the source root
directory:

.. code::

   $ cmake -S . -B build/benchmarks -DBUILD_BENCHMARK=ON
   $ cmake --build build/benchmarks

Then to run all the benchmarks:

.. code::

   $ build/benchmarks/./benchmark_fastscapelib

.. _google-benchmark: https://github.com/google/benchmark
