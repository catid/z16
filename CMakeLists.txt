cmake_minimum_required(VERSION 3.5)
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
set (CMAKE_CXX_STANDARD 11)
project(z16)


################################################################################
# Config

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# In debug mode, add -DDEBUG
add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
    # Warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -march=native -fstack-protector")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
endif()


################################################################################
# Dependencies

set(ZSTD_INCLUDE_FILES
    zstd_lib/zstd.h
)

file(GLOB ZSTD_SOURCE_FILES
    zstd_lib/common/*.c
    zstd_lib/common/*.h
    zstd_lib/compress/*.c
    zstd_lib/compress/*.h
    zstd_lib/decompress/*.c
    zstd_lib/decompress/*.h
    zstd_lib/dictBuilder/*.c
    zstd_lib/dictBuilder/*.h
)

# Zstd library
add_library(zstd ${ZSTD_SOURCE_FILES})
target_include_directories(zstd PUBLIC
    zstd_lib
)


################################################################################
# Targets

enable_testing()


# Library: z16

add_library(z16 STATIC
    include/z16.h
    src/z16.cpp
    src/Codec.hpp
    src/Codec.cpp
    src/Tools.hpp
    src/Tools.cpp
)
target_include_directories(z16 PUBLIC
    include
)
target_link_libraries(z16
    PUBLIC
        zstd
)

install(TARGETS z16 DESTINATION lib)
install(FILES include/z16.h DESTINATION include)


# App: z16_tool

add_executable(z16_tool
    app/z16_tool.cpp
    app/argagg.hpp
    app/lodepng.cpp
    app/lodepng.h
    app/teenypath.cpp
    app/teenypath.h
)
target_link_libraries(z16_tool
    PUBLIC
        z16
)

set_target_properties(z16_tool PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

install(TARGETS z16_tool DESTINATION bin)
