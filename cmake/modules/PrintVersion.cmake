cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

include(VersionUtils)

get_version_pieces(VERSION_PIECES)
format_version_pep440(VERSION_STR VERSION_PIECES)

message(${VERSION_STR})
