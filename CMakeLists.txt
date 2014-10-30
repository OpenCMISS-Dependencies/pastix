cmake_minimum_required(VERSION 3.0)
project(PASTIX VERSION 5.2.2.16 LANGUAGES C Fortran)

# no starpu yet
# no cuda support yet
#     not caring about DCUDA_SM_VERSION (in sopalin)
# no advanced mpi options (MPI_THREAD_MULTIPLE / MPI_Datatype)
# no semaphore barrier
# no hwloc
# unclear on murge options
# no scotch/metis support

if (OPENCMISS_DEPENDENCIES_LIBRARIES)
    SET(CMAKE_PREFIX_PATH ${OPENCMISS_DEPENDENCIES_CONFIGS_DIR} ${CMAKE_PREFIX_PATH})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OPENCMISS_DEPENDENCIES_LIBRARIES})
else()
    SET(CMAKE_PREFIX_PATH ../cmake/release ../cmake/debug ${CMAKE_PREFIX_PATH})
    SET(OPENCMISS_DEPENDENCIES_CONFIGS_DIR ${CMAKE_CURRENT_BINARY_DIR})
    SET(BUILD_TESTING ON)
    SET(USE_THREADS ON)
    SET(USE_CUDA OFF)
    SET(USE_METIS OFF)
    SET(USE_PTSCOTCH OFF)
endif()
include(CMakePackageConfigHelpers)

find_package(BLAS CONFIG REQUIRED)
find_package(MPI REQUIRED)

# Create the include directory
SET(PASTIX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
SET(PASTIX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PASTIX_INCLUDE_DIR})

SET(BUILD_DOUBLE_PREC ON)
SET(USE_STARPU OFF)
SET(USE_MAGMABLAS OFF) # BLAS for GPU?!

# Header file names
SET(VERSIONINT _int)
SET(VERSIONPRC _simple)
SET(VERSIONFLT _real)
SET(VERSIONTYPE ${VERSIONINT}${VERSIONPRC}${VERSIONFLT})

# it seems that sparsekit.f90 is not used
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -x f95-cpp-input")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
add_definitions(-DINTSIZE32)
if (BUILD_DOUBLE_PREC)
    add_definitions(-DPREC_DOUBLE)
endif()
if (NOT USE_THREADS)
    add_definitions(-DFORCE_NOSMP)
endif()
if (NOT USE_CUDA)
    add_definitions(-DFORCE_NO_CUDA)
endif()

if (USE_PTSCOTCH)
    find_package(PTSCOTCH CONFIG REQUIRED)
    add_definitions(-DWITH_SCOTCH -DDISTRIBUTED)
elseif (USE_SCOTCH)
    find_package(SCOTCH CONFIG REQUIRED)
    add_definitions(-DWITH_SCOTCH)
endif()
if (USE_METIS)
    find_package(METIS 4 CONFIG REQUIRED)
    add_definitions(-DMETIS)
endif()

add_subdirectory(src)