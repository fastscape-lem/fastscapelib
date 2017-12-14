# Helper script that
#  - reinstall C++ headers in the activated conda env,
#  - build and run the test suite (GTest)
#  - reinstall the Python package in the activated conda env
#
# Run this script from the (root) source directory.
YLW='\033[1;33m'
NC='\033[0m'

ACTIVE_ENV=$(conda info | grep 'default environment' | awk '{print $4}')
echo -e "${YLW}The packages will be installed in ${ACTIVE_ENV}${NC}"

echo -e "${YLW}Install the C++ library${NC}"
rm -rf build; mkdir build; cd build;
cmake -DCMAKE_INSTALL_PREFIX=$ACTIVE_ENV ..
make install

echo -e "${YLW}Build and run tests${NC}"
cmake -DBUILD_TESTS=ON ..
make run_tests

echo -e "${YLW}Install the Python package${NC}"
cd ../fastscape-python/
pip uninstall fastscape --yes
pip install .

cd ..
