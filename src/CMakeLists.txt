add_subdirectory(utils/src)

macro(IMPORT WHAT)
    STRING(TOUPPER "${WHAT}_SRC" VARNAME)
    SET(DIR ${WHAT}/src)
    add_subdirectory(${DIR})
    # Get variables there
    get_directory_property(REGULAR_FILES
        DIRECTORY ${DIR}
        DEFINITION ${VARNAME})
    APPEND_PREFIX(${DIR} REGULAR_FILES)
    # Get multiarch-files and set -DWITH_TYPE_PREFIX for them
    get_directory_property(MULTIARCH_FILES
        DIRECTORY ${DIR}
        DEFINITION ${VARNAME}_MA)
    APPEND_PREFIX(${DIR} MULTIARCH_FILES)
    
    add_library(${WHAT} OBJECT ${REGULAR_FILES} ${MULTIARCH_FILES})
    #message(STATUS "target_compile_definitions(${WHAT} PRIVATE ${COMMON_DEFINITIONS} CHOL_SOPALIN)")
    target_compile_definitions(${WHAT} 
        PRIVATE ${COMMON_DEFINITIONS} CHOL_SOPALIN)
    target_include_directories(${WHAT} 
        PRIVATE 
        symbol/src common/src order/src sopalin/src perf/src
        fax/src blend/src perf/src/num_recipes kass/src
        PUBLIC
        ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
        $<$<BOOL:${USE_PTSCOTCH}>:
            $<TARGET_PROPERTY:ptscotch,INTERFACE_INCLUDE_DIRECTORIES>>
        $<$<BOOL:${USE_PTSCOTCH}>:
            $<TARGET_PROPERTY:scotch,INTERFACE_INCLUDE_DIRECTORIES>>
        $<$<BOOL:${USE_SCOTCH}>:
            $<TARGET_PROPERTY:scotch,INTERFACE_INCLUDE_DIRECTORIES>>
         $<$<BOOL:${USE_METIS}>:
             $<TARGET_PROPERTY:metis,INTERFACE_INCLUDE_DIRECTORIES>>
    )
    LIST(APPEND ALL_OBJS "$<TARGET_OBJECTS:${WHAT}>")
    if (${WHAT} STREQUAL sopalin)
        get_directory_property(MULTIFAC_FILES
            DIRECTORY ${DIR}
            DEFINITION SRC_FAC)
        APPEND_PREFIX(${DIR} MULTIFAC_FILES)
        # Special case for different factorization codes (per precision)
        # GE: SOPALIN_LU
        # SY: no extra definition, using NOEXTRADEF_SY to have a loop step
        # HE: HERMITIAN
        foreach(extradef SOPALIN_LU NOEXTRADEF_SY HERMITIAN)
            SET(lib sopalin_${extradef})
            add_library(${lib} OBJECT ${MULTIFAC_FILES})
            #message(STATUS "target_compile_definitions(${lib} PRIVATE ${extradef} ${COMMON_DEFINITIONS})")
            target_compile_definitions(${lib} PRIVATE
                ${extradef} ${COMMON_DEFINITIONS})
            target_include_directories(${lib} PRIVATE 
                order/src common/src blend/src symbol/src
                kass/src fax/src perf/src
                ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
            )
            # Make sure the type specific headers are used
            #target_include_directories(${lib} BEFORE PRIVATE ${PASTIX_INCLUDE_DIR_PREC})
            LIST(APPEND ALL_OBJS "$<TARGET_OBJECTS:${lib}>") 
        endforeach()
    endif()
    
    foreach(prec s d c z)
        SETTYPEDEFS(${prec})
        SETINCDIR(${prec})
        list(LENGTH MULTIARCH_FILES _NUM_MA_FILES)
        if(_NUM_MA_FILES GREATER 0 AND BUILD_PRECISION MATCHES [${prec}])
            SET(lib ${prec}${WHAT})
            add_library(${lib} OBJECT ${MULTIARCH_FILES})
            target_compile_definitions(${lib}
                PRIVATE ${COMMON_DEFINITIONS} ${TYPEDEFS} CHOL_SOPALIN)
            #message(STATUS "target_compile_definitions(${lib} PRIVATE ${COMMON_DEFINITIONS} ${TYPEDEFS} CHOL_SOPALIN)")
            target_include_directories(${lib} PRIVATE 
                order/src common/src blend/src symbol/src
                kass/src fax/src perf/src
                ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
                $<$<BOOL:${USE_PTSCOTCH}>:
                    $<TARGET_PROPERTY:ptscotch,INTERFACE_INCLUDE_DIRECTORIES>>
                $<$<BOOL:${USE_PTSCOTCH}>:
                    $<TARGET_PROPERTY:scotch,INTERFACE_INCLUDE_DIRECTORIES>>
                $<$<BOOL:${USE_SCOTCH}>:
                    $<TARGET_PROPERTY:scotch,INTERFACE_INCLUDE_DIRECTORIES>>
                $<$<BOOL:${USE_METIS}>:
                    $<TARGET_PROPERTY:metis,INTERFACE_INCLUDE_DIRECTORIES>>
            )
            
            # Make sure the type specific headers are used
            #target_include_directories(${lib} BEFORE PRIVATE ${PASTIX_INCLUDE_DIR_PREC})
            LIST(APPEND ALL_OBJS "$<TARGET_OBJECTS:${lib}>")
        endif()
        if (${WHAT} STREQUAL sopalin
            AND BUILD_PRECISION MATCHES [${prec}])
            get_directory_property(MULTIFAC_FILES
                DIRECTORY ${DIR}
                DEFINITION SRC_FAC)
            APPEND_PREFIX(${DIR} MULTIFAC_FILES)
            # Special case for different factorization codes (per precision)
            # GE: SOPALIN_LU
            # SY: no extra definition, using NOEXTRADEF_SY to have a loop step
            # HE: HERMITIAN
            foreach(extradef SOPALIN_LU NOEXTRADEF_SY HERMITIAN)
                SET(lib ${prec}sopalin_${extradef})
                add_library(${lib} OBJECT ${MULTIFAC_FILES})
                target_compile_definitions(${lib} PRIVATE
                    ${extradef} ${COMMON_DEFINITIONS} ${TYPEDEFS})
                #message(STATUS "target_compile_definitions(${lib} PRIVATE ${extradef} ${COMMON_DEFINITIONS} ${TYPEDEFS} CHOL_SOPALIN)")
                target_include_directories(${lib} PRIVATE 
                    order/src common/src blend/src symbol/src
                    kass/src fax/src perf/src
                    ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
                )
                # Make sure the type specific headers are used
                #target_include_directories(${lib} BEFORE PRIVATE ${PASTIX_INCLUDE_DIR_PREC})
                LIST(APPEND ALL_OBJS "$<TARGET_OBJECTS:${lib}>") 
            endforeach()
        endif()
    endforeach()
endmacro()

macro(APPEND_PREFIX PREFIX FILELIST)
    SET(outlist )
    FOREACH(root ${${FILELIST}})
      LIST(APPEND outlist ${PREFIX}/${root})
    ENDFOREACH(root)
    SET(${FILELIST} ${outlist})
endmacro()

SET(ALL_OBJS )
IMPORT(common)
IMPORT(symbol)
IMPORT(order)
IMPORT(sopalin)
IMPORT(blend)
IMPORT(fax)
IMPORT(kass)
IMPORT(sparse-matrix)
    
add_library(pastix ${ALL_OBJS})
set_target_properties(pastix PROPERTIES OUTPUT_NAME pastix-${PASTIX_VERSION})
target_compile_definitions(pastix PRIVATE
    ${COMMON_DEFINITIONS}
    ${TYPE_DEFINITIONS}
)
#message(STATUS "target_compile_definitions(pastix PRIVATE ${COMMON_DEFINITIONS} ${TYPE_DEFINITIONS})")    
if (USE_THREADS)
    #if (MINGW)
    #    target_link_libraries(pastix PRIVATE pthread)#${CMAKE_THREAD_LIBS_INIT}
    #else()
        target_link_libraries(pastix PRIVATE Threads::Threads)#${CMAKE_THREAD_LIBS_INIT}
    #endif()
endif()
target_link_libraries(pastix PRIVATE blas PUBLIC ${GRAPH_LIBS})
if (NOT CMAKE_C_COMPILER STREQUAL MPI_C_COMPILER)
    target_link_libraries(pastix PUBLIC $<BUILD_INTERFACE:mpi>)
endif()
    
target_include_directories(pastix 
    PRIVATE common/src symbol/src order/src blend/src sopalin/src
    perf/src perf/src/num_recipes fax/src kass/src
    matrix_drivers/src
)

SETINCDIR(${NEUTRAL_PRECISION})
install(TARGETS pastix EXPORT pastix-config
    DESTINATION lib
    INCLUDES DESTINATION ${PASTIX_INCDIR_SUFFIX_PREC})

add_subdirectory(perf/src)
add_subdirectory(matrix_drivers/src)

if (BUILD_TESTS)
    add_subdirectory(example/src)
endif()