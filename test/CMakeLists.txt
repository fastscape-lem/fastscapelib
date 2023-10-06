cmake_minimum_required(VERSION 3.15)

# Dependencies
# ============

if(FS_DOWNLOAD_GTEST OR FS_GTEST_SRC_DIR)
    include(FetchContent)

    if(FS_DOWNLOAD_GTEST)
        message(STATUS "Downloading googletest v1.13.0")

        FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest
            GIT_TAG tags/v1.13.0
            GIT_SHALLOW TRUE)
    else()
        message(STATUS "Build googletest from local directory")

        FetchContent_Declare(googletest
            SOURCE_DIR FS_GTEST_SRC_DIR
            )
    endif()

    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" FORCE)

    FetchContent_MakeAvailable(googletest)

    # do not install googletest by default when running cmake --install
    set_property(DIRECTORY ${googletest_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
else()
    find_package(GTest REQUIRED)
endif()

include(GoogleTest)

# Compiler configuration
# ======================

# -- Set compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT WIN32))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wreorder")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wunused-variable")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj")
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(NOT WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wreorder")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Wunused-variable")
    else()   # We are using clang-cl
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /bigobj")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST:NO")
        add_definitions(-D_CRT_SECURE_NO_WARNINGS)
        add_definitions(-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
    endif()
endif()

# Build
# =====

set(FASTSCAPELIB_TEST_SRC
    test_basin_graph.cpp
    test_diffusion_adi.cpp
    test_flow_graph.cpp
    test_flow_operator.cpp
    test_flow_router.cpp
    test_flow_snapshot.cpp
    test_single_flow_router.cpp
    test_multi_flow_router.cpp
    test_pflood.cpp
    test_profile_grid.cpp
    test_raster_grid.cpp
    test_raster_rook_grid.cpp
    test_raster_queen_grid.cpp
    test_raster_bishop_grid.cpp
    test_structured_grid.cpp
    test_sink_resolver.cpp
    test_spl.cpp
    test_trimesh.cpp
)

enable_testing()

# -- build a target for each test
foreach(filename IN LISTS FASTSCAPELIB_TEST_SRC)
    string(REPLACE ".cpp" "" targetname ${filename})
    add_executable(${targetname} ${filename})
    target_link_libraries(${targetname} PRIVATE fastscapelib GTest::gtest_main)
    target_compile_features(${targetname} PRIVATE cxx_std_17)
    gtest_discover_tests(${targetname})
endforeach()
