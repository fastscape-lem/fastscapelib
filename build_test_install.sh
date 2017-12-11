# Helper script that
#  - reinstall C++ headers in the activated conda env,
#  - build and run the test suite (GTest)
#  - reinstall the Python package in the activated conda env
#
# Run this script from the (root) source directory.

ACTIVE_ENV=$(conda info | grep 'default environment' | awk '{print $4}')

rm -rf build; mkdir build; cd build;
cmake -DCMAKE_INSTALL_PREFIX=$ACTIVE_ENV ..
make install

cmake -DBUILD_TESTS=ON ..
make run_tests

cd ../fastscape-python/
pip uninstall fastscape --yes
pip install .

cd ..
