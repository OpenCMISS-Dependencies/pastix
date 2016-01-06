cmake_minimum_required(VERSION 3.3)
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

# Instructions for PASTIX builds on windows:
# http://forge-bd1.inria.fr/forum/message.php?msg_id=148918&group_id=186
# PThreads from ftp://sourceware.org/pub/pthreads-win32/pthreads-w32-2-9-1-release.zip

option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
option(PACKAGE_CONFIG_DIR "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)" "lib/cmake")
set(BUILD_PRECISION sdcz CACHE STRING "${PROJECT_NAME} - Build precisions (sdcz)")
set(INT_TYPE int32 CACHE STRING "${PROJECT_NAME} - Integer version (int32, long, int64)")
option(USE_THREADS "Use PTHREAD library" ON)
option(USE_METIS "Use METIS orderings" ON)
option(USE_PTSCOTCH "Use PT-Scotch orderings (will include Scotch)" ON)
option(USE_SCOTCH "Use Scotch orderings" OFF)
option(WITH_OPENMP "${PROJECT_NAME} - Use OpenMP" ON)
if (BUILD_SHARED_LIBS)
    # This will cause all the OBJECT targets to have -fPIC as well.
    SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Not yet implemented to be ON
SET(USE_CUDA OFF)
SET(USE_STARPU OFF)
SET(USE_MAGMABLAS OFF) # BLAS for GPU?!

include(CMakePackageConfigHelpers)

find_package(BLAS ${BLAS_VERSION} REQUIRED)
find_package(MPI REQUIRED)

if (WITH_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  else()
    SET(WITH_OPENMP FALSE)
    message(WARNING "OpenMP requested but not found. Skipping.")
  endif()
endif()
if (USE_THREADS)
    find_package(Threads QUIET)
    if (Threads_FOUND)
        if (MINGW)
            add_definitions(-D_MINGW_)
        endif()
    else()
        message(WARNING "Threading was requested (USE_THREADS=${USE_THREADS}), but package could not be found. Disabling.")
        SET(USE_THREADS NO)
    endif()
endif()

if (BUILD_TESTS)
    enable_testing()
endif()

#include(CheckFortranCompilerFlag)
# This bit surely does not work as some files dont have Fortran free form
# See https://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html
#CHECK_Fortran_COMPILER_FLAG(-ffree-form HAVE_Fortran_FREEFORM_FLAG)
#if (HAVE_Fortran_FREEFORM_FLAG)
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")
#endif()
#CHECK_Fortran_COMPILER_FLAG(-x HAVE_Fortran_X_FLAG)
#if (HAVE_Fortran_X_FLAG)
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -x")
#endif()
#CHECK_Fortran_COMPILER_FLAG(f95-cpp-input HAVE_Fortran_f95_cpp_FLAG)
#if (HAVE_Fortran_f95_cpp_FLAG)
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} f95-cpp-input")
#endif()

# Detect host architecture
# Yet undetermined:
# Compaq: alpha_compaq_osf1
# SUN: sun
# AIX: power_ibm_aix
if (APPLE)
    SET(HOSTARCH i686_mac)
elseif(UNIX)
    SET(HOSTARCH i686_pc_linux)
endif()
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

# Windows hacks
if (WIN32)
    cmake_host_system_information(RESULT NUMPROC QUERY NUMBER_OF_LOGICAL_CORES)
    configure_file(windows_fix.h.in ${CMAKE_CURRENT_BINARY_DIR}/windows_fix.h @ONLY)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/windows_fix.h DESTINATION include)#/pastix
endif()

SET(PASTIX_INCDIR_SUFFIX include/pastix/${INT_TYPE})
SET(PASTIX_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${PASTIX_INCDIR_SUFFIX})
SET(PASTIX_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(MAKE_DIRECTORY ${PASTIX_INCLUDE_DIR})

macro(SETINCDIR precision)
    SET(PASTIX_INCLUDE_DIR_PREC ${PASTIX_INCLUDE_DIR}/${precision})
    SET(PASTIX_INCDIR_SUFFIX_PREC ${PASTIX_INCDIR_SUFFIX}/${precision})
    file(MAKE_DIRECTORY ${PASTIX_INCLUDE_DIR_PREC})
endmacro()
macro(SETTYPEDEFS precision)
    SET(TYPEDEFS WITH_TYPE_PREFIX)
    #Somehow "${precision} MATCHES [dz]" wont work for "d" at the first cmake run; after the second it works.
    if (${precision} STREQUAL d OR ${precision} STREQUAL z) 
        LIST(APPEND TYPEDEFS PREC_DOUBLE)
    endif()
    if (${precision} MATCHES [cz])
        LIST(APPEND TYPEDEFS TYPE_COMPLEX)     
    endif()
    #message(STATUS "Call to SETTYPEDEFS with prec '${precision}' gave ${TYPEDEFS} (${precision} MATCHES [dz])")
endmacro()

if (NOT USE_THREADS)
    add_definitions(-DFORCE_NOSMP)
endif()
if (NOT USE_CUDA)
    add_definitions(-DFORCE_NO_CUDA)
endif()

set(GRAPH_LIBS )
if (USE_PTSCOTCH)
    find_package(PTSCOTCH ${PTSCOTCH_VERSION} REQUIRED)
    add_definitions(-DWITH_SCOTCH -DDISTRIBUTED)
    message(STATUS "Using PT-Scotch ${PTSCOTCH_VERSION} orderings")
    list(APPEND GRAPH_LIBS ptscotch)
elseif (USE_SCOTCH)
    find_package(SCOTCH ${SCOTCH_VERSION} REQUIRED)
    add_definitions(-DWITH_SCOTCH)
    message(STATUS "Using Scotch ${SCOTCH_VERSION} orderings")
    list(APPEND GRAPH_LIBS scotch)
endif()

if (USE_METIS)
    # Fixed up the pastix.c file to have it compile against Metis 5.
    # Just a couple of option renames had to be done, see pastix.c:1650
    find_package(METIS ${METIS_VERSION} REQUIRED)
    add_definitions(-DMETIS)
    message(STATUS "Using Metis ${METIS_VERSION} orderings")
    list(APPEND GRAPH_LIBS metis)
endif()

foreach(prec d s c z)
    if(BUILD_PRECISION MATCHES [${prec}])
        # We need a neutral precision for notype-includes
        if(NOT NEUTRAL_PRECISION)
            SET(NEUTRAL_PRECISION ${prec})
        endif()
        message(STATUS "Building pastix precision '${prec}'")
        
        # Check which header files need to be generated
        SETINCDIR(${prec})
        if(NOT EXISTS ${PASTIX_INCLUDE_DIR_PREC}/pastix.h)
            SET(NEED_HEADERS_${prec} YES)
            message(STATUS "Generating headers for precision '${prec}'")
        endif()
    endif()
endforeach()
if(NOT NEUTRAL_PRECISION)
    message(FATAL_ERROR "No build precision set. Please set BUILD_PRECISION accordingly")
endif()

add_subdirectory(src)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION include)
install(EXPORT pastix-config DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/pastix-config-version.cmake COMPATIBILITY AnyNewerVersion)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake 
    "include(CMakeFindDependencyMacro)\r\n"
    	"set(${PROJECT_NAME}_IMPORT_PREFIX \${_IMPORT_PREFIX})\r\n"
)
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "find_dependency(BLAS ${BLAS_VERSION})\r\n")
if (USE_METIS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "find_dependency(METIS ${METIS_VERSION})\r\n")
endif()
if (USE_SCOTCH)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "find_dependency(SCOTCH ${SCOTCH_VERSION})\r\n")
endif()
if (USE_PTSCOTCH)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "find_dependency(PTSCOTCH ${PTSCOTCH_VERSION})\r\n")
endif()
if (USE_THREADS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "find_dependency(Threads)\r\n")
endif()
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})")
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/pastix-config-dependencies.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR})