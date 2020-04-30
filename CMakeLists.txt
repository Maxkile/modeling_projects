cmake_minimum_required(VERSION 3.5)

project(meshbuilder LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(OpenMP REQUIRED)

set(INCLUDE_DIRECTORIES include)
set(SRC_DIRECTORIES src)

include_directories(meshbuilder ${INCLUDE_DIRECTORIES})

set(HEADERS_ONLY_INCLUDE ${INCLUDE_DIRECTORIES}/solver.hpp ${INCLUDE_DIRECTORIES}/decomposition.hpp ${INCLUDE_DIRECTORIES}/toposBuild.hpp)

set(HEADERS_COMPILE ${INCLUDE_DIRECTORIES}/FixedSizeMeshContainer.hpp  ${INCLUDE_DIRECTORIES}/IO.hpp ${INCLUDE_DIRECTORIES}/platformDependencies.hpp
    ${INCLUDE_DIRECTORIES}/Sparse.hpp ${INCLUDE_DIRECTORIES}/stdafx.hpp
    ${INCLUDE_DIRECTORIES}/VariableSizeMeshContainer.hpp ${INCLUDE_DIRECTORIES}/vmo.hpp ${INCLUDE_DIRECTORIES}/vtkGenerator.hpp)

set(SOURCES ${SRC_DIRECTORIES}/solver.cpp ${SRC_DIRECTORIES}/decomposition.cpp ${SRC_DIRECTORIES}/prod.cpp ${SRC_DIRECTORIES}/toposBuild.cpp
    ${SRC_DIRECTORIES}/vmo.cpp ${SRC_DIRECTORIES}/IO.cpp)


add_executable(meshbuilder ${HEADERS_COMPILE} ${SOURCES})

target_link_libraries(meshbuilder "${OpenMP_CXX_FLAGS}")





