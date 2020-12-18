cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

include(VersionUtils)

set_version_str(${INCLUDEDIR})
message(${VERSION_STR})
