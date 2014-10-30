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
# building only parallel version


if (OPENCMISS_DEPENDENCIES_LIBRARIES)
    SET(CMAKE_PREFIX_PATH ${OPENCMISS_DEPENDENCIES_CONFIGS_DIR} ${CMAKE_PREFIX_PATH})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OPENCMISS_DEPENDENCIES_LIBRARIES})
else()
    SET(CMAKE_PREFIX_PATH ../cmake/release ../cmake/debug ${CMAKE_PREFIX_PATH})
    SET(OPENCMISS_DEPENDENCIES_CONFIGS_DIR ${CMAKE_CURRENT_BINARY_DIR})
    SET(BUILD_TESTING ON)
    SET(USE_THREADS ON)
    SET(USE_CUDA OFF)
    SET(USE_METIS ON)
    SET(USE_PTSCOTCH ON)
endif()
include(CMakePackageConfigHelpers)

find_package(BLAS CONFIG REQUIRED)
find_package(MPI REQUIRED)

# build double version (s,d,c,z)
SET(BUILD_TYPE d)
# integer versions (int32, long, int64)
SET(INT_TYPE int32)
SET(USE_STARPU OFF)
SET(USE_MAGMABLAS OFF) # BLAS for GPU?!

# Ensure the include directory exists
SET(PASTIX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
SET(PASTIX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PASTIX_INCLUDE_DIR})

SET(HOSTARCH i686_pc_linux)

# it seems that sparsekit.f90 is not used
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -x f95-cpp-input")

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_Fortan_COMPILER ${MPI_Fortran_COMPILER})
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
SET(COMMON_DEFINITIONS VERSION=opencmiss X_ARCH${HOSTARCH} DOF_CONSTANT)
########################################################################
# Check for variable types and precision
########################################################################
SET(TYPE_DEFINITIONS )
if (INT_TYPE STREQUAL long)
    SET(VERSIONINT _long)
    LIST(APPEND TYPE_DEFINITIONS FORCE_LONG INTSIZELONG)
    message(STATUS "Building with long integers..")
elseif(INT_TYPE STREQUAL int64)
    SET(VERSIONINT _int64)
    LIST(APPEND TYPE_DEFINITIONS INTSIZE64)
    message(STATUS "Building with int64 integers..")
elseif(INT_TYPE STREQUAL int32)
    SET(VERSIONINT _int32)
    LIST(APPEND TYPE_DEFINITIONS INTSIZE32)
    message(STATUS "Building with int32 integers..")
else()
    message(FATAL_ERROR "Unknown integer type: ${INT_TYPE}")
endif()
SET(VERSIONPRC _simple)
if(BUILD_TYPE STREQUAL s OR BUILD_TYPE STREQUAL d)
    message(STATUS "Building double precision..")
    LIST(APPEND TYPE_DEFINITIONS PREC_DOUBLE)
    SET(VERSIONPRC _double)
endif()
SET(VERSIONFLT _real)
if(BUILD_TYPE STREQUAL c OR BUILD_TYPE STREQUAL z)
    message(STATUS "Building complex arithmetics..")
    LIST(APPEND TYPE_DEFINITIONS TYPE_COMPLEX)
    SET(VERSIONFLT _complex)
endif()
SET(VERSIONTYPE ${VERSIONINT}${VERSIONPRC}${VERSIONFLT})
########################################################################

if (NOT USE_THREADS)
    add_definitions(-DFORCE_NOSMP)
endif()
if (NOT USE_CUDA)
    add_definitions(-DFORCE_NO_CUDA)
endif()

if (USE_PTSCOTCH)
    find_package(PTSCOTCH CONFIG REQUIRED)
    add_definitions(-DWITH_SCOTCH -DDISTRIBUTED)
    message(STATUS "Using PT-Scotch ${PTSCOTCH_VERSION} orderings")
elseif (USE_SCOTCH)
    find_package(SCOTCH CONFIG REQUIRED)
    add_definitions(-DWITH_SCOTCH)
    message(STATUS "Using Scotch ${SCOTCH_VERSION} orderings")
endif()
if (USE_METIS)
    find_package(METIS 4 CONFIG REQUIRED)
    add_definitions(-DMETIS)
    message(STATUS "Using Metis ${METIS_VERSION} orderings")
endif()

add_subdirectory(src)