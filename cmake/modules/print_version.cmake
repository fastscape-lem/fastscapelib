cmake_minimum_required(VERSION 3.1)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})

include(GetGitRevisionAddons)

get_git_version_pieces(GIT_VERSION_PIECES)

format_version_pep440(VERSION_STR GIT_VERSION_PIECES)

message(${VERSION_STR})
