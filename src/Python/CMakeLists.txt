cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(pythonwrapper)

#Not tested
find_package(pybind11 REQUIRED)
pybind11_add_module(pythonwrapper pythonwrapper.cpp)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # using Clang
  target_compile_options(pythonwrapper PRIVATE -march=native PRIVATE -Wall PRIVATE -Wextra PRIVATE -flto)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  target_compile_options(pythonwrapper PRIVATE -march=native PRIVATE -Wall PRIVATE -Wextra PRIVATE -flto)
# elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
  # change architecture settings depending target
  target_compile_options(pythonwrapper PRIVATE /arch:AVX2 PRIVATE /W4 PRIVATE /GL)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # using Clang
    target_compile_options(pythonwrapper PRIVATE -Ofast)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    target_compile_options(pythonwrapper PRIVATE -Ofast)
  # elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # using Intel C++ 
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
    # change architecture settings depending target
    target_compile_options(pythonwrapper PRIVATE /Ot)
  endif()
endif()