cmake_minimum_required(VERSION 3.0)
project(SqLatticePercolation)

set(CMAKE_CXX_STANDARD 11)

#set (CMAKE_C_COMPILER               /usr/bin/gcc)
#set (CMAKE_CXX_COMPILER             /home/shahnoor/software/pgi/linux86-64-llvm/2018/bin/pgc++)
#set (CMAKE_MAKE_PROGRAM             /usr/bin/make)
SET( CMAKE_CXX_FLAGS  "-pthread -fopenmp")
#SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

set(SOURCE_FILES
        src/main.cpp
        src/types.h
        src/exception/exceptions.h
        src/index/index.cpp
        src/index/index.h
        src/lattice/bond.cpp
        src/lattice/bond.h
        src/lattice/bond_type.h
        src/lattice/lattice.cpp
        src/lattice/lattice.h
        src/lattice/site.cpp
        src/lattice/site.h
        src/percolation/cluster.cpp
        src/percolation/cluster.h
        src/percolation/percolation.cpp
        src/percolation/percolation.h
        src/percolation/percolation_site_ballistic_deps_v2.cpp
        src/util/printer.h
        src/util/time_tracking.cpp
        src/util/time_tracking.h
        src/util/printer.cpp
        src/percolation/percolation_site_v9.cpp
        src/tests/test_percolation.cpp
        src/tests/test_percolation.h)

add_executable(SqLatticePercolation ${SOURCE_FILES})
