cmake_minimum_required(VERSION 3.0)
project(PASTIX VERSION 5.2.2.16 LANGUAGES C Fortran CXX)

# no starpu yet
# no cuda support yet
#     not caring about DCUDA_SM_VERSION (in sopalin)
# no advanced mpi options (MPI_THREAD_MULTIPLE / MPI_Datatype)
# no semaphore barrier
# no hwloc
# unclear on murge options
# building only parallel version

# Hints:
# src/* builds libpastix
# src/utils creates the different headers for each precision
# src/matrix_drivers uses those headers for the matrix_drivers library
# So far: unclear how pastix is going to be used, so we go with this version now.

option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
set(BUILD_PRECISION sdcz CACHE STRING "${PROJECT_NAME} - Build precisions (sdcz)")
set(INT_TYPE int32 CACHE STRING "${PROJECT_NAME} - Integer version (int32, long, int64)")
option(USE_THREADS "Use PTHREAD library" ON)
option(USE_METIS "Use METIS orderings" ON)
option(USE_PTSCOTCH "Use PT-Scotch orderings (will include Scotch)" ON)
option(USE_SCOTCH "Use Scotch orderings" OFF)
# Not yet implemented to be ON
SET(USE_CUDA OFF)
SET(USE_STARPU OFF)
SET(USE_MAGMABLAS OFF) # BLAS for GPU?!

include(CMakePackageConfigHelpers)

find_package(BLAS CONFIG REQUIRED)
find_package(MPI REQUIRED)

if (BUILD_TESTS)
    enable_testing()
endif()

# it seems that sparsekit.f90 is not used
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -x f95-cpp-input")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_Fortan_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
SET(COMMON_DEFINITIONS 
    VERSION=opencmiss
    X_ARCH${HOSTARCH}
    DOF_CONSTANT)
########################################################################
# Check for variable types and precision
########################################################################
if (INT_TYPE STREQUAL long)
    SET(VERSIONINT _long)
    LIST(APPEND COMMON_DEFINITIONS FORCE_LONG INTSIZELONG)
    message(STATUS "Building with long integers..")
elseif(INT_TYPE STREQUAL int64)
    SET(VERSIONINT _int64)
    LIST(APPEND COMMON_DEFINITIONS INTSIZE64)
    message(STATUS "Building with int64 integers..")
elseif(INT_TYPE STREQUAL int32)
    SET(VERSIONINT _int32)
    LIST(APPEND COMMON_DEFINITIONS INTSIZE32)
    message(STATUS "Building with int32 integers..")
else()
    message(FATAL_ERROR "Unknown integer type: ${INT_TYPE}")
endif()
########################################################################

# Ensure the include directory exists
SET(PASTIX_INCDIR_SUFFIX include/pastix/${INT_TYPE})
SET(PASTIX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${PASTIX_INCDIR_SUFFIX})
SET(PASTIX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PASTIX_INCLUDE_DIR})

macro(SETINCDIR precision)
    SET(PASTIX_INCLUDE_DIR_PREC ${PASTIX_INCLUDE_DIR}/${precision})
    SET(PASTIX_INCDIR_SUFFIX_PREC ${PASTIX_INCDIR_SUFFIX}/${precision})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PASTIX_INCLUDE_DIR_PREC})
endmacro()
macro(SETTYPEDEFS precision)
    SET(TYPEDEFS WITH_TYPE_PREFIX)
    if (${precision} MATCHES [dz])
        LIST(APPEND TYPEDEFS PREC_DOUBLE)
    endif()
    if (${precision} MATCHES [cz])
        LIST(APPEND TYPEDEFS TYPE_COMPLEX)     
    endif()
    #message(STATUS "Call to SETTYPEDEFS with prec '${precision}' gave ${TYPEDEFS}")
endmacro()

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
    # Fixed up the pastix.c file to have it compile against Metis 5.
    # Just a couple of option renames had to be done, see pastix.c:1650
    find_package(METIS 5 CONFIG REQUIRED)
    add_definitions(-DMETIS)
    message(STATUS "Using Metis ${METIS_VERSION} orderings")
endif()

foreach(prec s d c z)
    if(BUILD_PRECISION MATCHES [${prec}])
        # We need a neutral precision for notype-includes
        if(NOT NEUTRAL_PRECISION)
            SET(NEUTRAL_PRECISION ${prec})
        endif()
        message(STATUS "Building pastix precision '${prec}'")
    endif()
endforeach()
if(NOT NEUTRAL_PRECISION)
    message(FATAL_ERROR "No build precision set. Please set BUILD_PRECISION accordingly")
endif()

add_subdirectory(src)

install(DIRECTORY include/ DESTINATION include)
install(EXPORT pastix-config DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/pastix-config-version.cmake
    VERSION ${METIS_VERSION} COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-version.cmake DESTINATION lib/cmake)