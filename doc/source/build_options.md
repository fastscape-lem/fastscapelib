(build-options)=

# Build and Configuration

## Include Fastscapelib in a CMake Project

After {ref}`installing Fastscapelib C++ headers <install-cpp>` (in a default
location), you should be able to use CMake's
[find_package](https://cmake.org/cmake/help/latest/command/find_package.html) to
ensure these will be found when building your CMake project:

```cmake
find_package(fastscapelib REQUIRED)
```

Don't forget to link your target library or application with Fastscapelib:

```cmake
target_link_libraries(your_target INTERFACE fastscapelib)
```

## Build Options

Fastscapelib provides the following CMake build options (all disabled by
default). See below for more explanations.

```{list-table}
:widths: 25 75

* - ``FS_BUILD_TESTS``
  - Enables {ref}`building the C++ tests <run-cpp-tests>`
* - ``FS_DOWNLOAD_GTEST``
  - Downloads google-test and builds it locally instead of using a binary
    installation
* - ``FS_GTEST_SRC_DIR``
  - Indicates where to find the google-test sources instead of downloading
    them
* - ``FS_BUILD_BENCHMARKS``
  - Enables {ref}`building the micro-benchmarks <run-benchmarks>`
* - ``FS_DOWNLOAD_GBENCHMARK``
  - Downloads google-benchmark and builds it locally instead of using a binary
    installation
* - ``FS_GBENCHMARK_SRC_DIR``
  - Indicates where to find the google-benchmark sources instead of downloading
    them
* - ``FS_DOWNLOAD_XTENSOR``
  - Downloads xtensor development version (master branch on GitHub) and uses
    it to build fastscapelib (useful for testing)
```

(run-cpp-tests)=

## Build and Run the C++ Tests

Fastscapelib has a test suite based on [google-test].

You can install google-test, e.g., using conda:

```
$ conda install gtest -c conda-forge
```

Alternatively, google-test may be downloaded automatically by enabling
`FS_DOWNLOAD_GTEST`, or a custom install path may be given by setting
`FS_GTEST_SRC_DIR` (setting `FS_DOWNLOAD_GTEST=ON` or
`FS_GTEST_SRC_DIR=/path/to/gtest` automatically sets `FS_BUILD_TESTS=ON`).

:::{tip}

Download and build google-test from source as part of the Fastscapelib build
process may prevent possible issues with reusing a pre-installed version of
google-test that has been built with different options or flags.

:::

To build the tests, run the following commands from the source root directory:

```
$ cmake -S . -B build/tests -DFS_BUILD_TESTS=ON
$ cmake --build build/tests
```

Then to run all the tests:

```
$ ctest -T test --output-on-failure build/tests
```

(run-python-tests)=

## Run the Python Tests

Running the Python tests requires [pytest]. You can install it using, e.g.,
conda:

```
$ conda install pytest -c conda-forge
```

After {ref}`(re)installing the Fastscapelib Python library <install-python>`,
you can run the tests using the following command (from the repository root
directory):

```
$ pytest -v .
```

(run-benchmarks)=

## Build and Run the Benchmarks

Fastscapelib has also a micro-benchmark suite based on [google-benchmark].

You can install google-benchmark, e.g., using conda:

```
$ conda install benchmark -c conda-forge
```

Alternatively, google-benchmark may be downloaded automatically by enabling
`FS_DOWNLOAD_GBENCHMARK`, or a custom install path may be given by setting
`FS_GBENCHMARK_SRC_DIR` (setting `FS_DOWNLOAD_GBENCHMARK=ON` or
`FS_GBENCHMARK_SRC_DIR=/path/to/gbenchmark` automatically sets
`FS_BUILD_BENCHMARKS=ON`).

To build the benchmarks, run the following commands from the source root
directory:

```
$ cmake -S . -B build/benchmarks -DFS_BUILD_BENCHMARKS=ON
$ cmake --build build/benchmarks
```

Then to run all the benchmarks:

```
$ build/benchmarks/./benchmark_fastscapelib
```

[google-benchmark]: https://github.com/google/benchmark
[google-test]: https://github.com/google/googletest
[pytest]: https://docs.pytest.org/
