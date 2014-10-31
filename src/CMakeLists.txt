add_subdirectory(utils/src)

macro(IMPORT DIR VARNAME)
    add_subdirectory(${DIR})
    get_directory_property(REGULAR_FILES
        DIRECTORY ${DIR}
        DEFINITION ${VARNAME})
    APPEND_PREFIX(${DIR} REGULAR_FILES)
    # Get multiarch-files and set -DWITH_TYPE_PREFIX for them
    get_directory_property(MULTIARCH_FILES
        DIRECTORY ${DIR}
        DEFINITION ${VARNAME}_MA)
    APPEND_PREFIX(${DIR} MULTIARCH_FILES)
    set_source_files_properties(${MULTIARCH_FILES}
        PROPERTIES COMPILE_DEFINITIONS WITH_TYPE_PREFIX)
    # join into one variable for ease of handling
    SET(${VARNAME} ${REGULAR_FILES} ${MULTIARCH_FILES})
endmacro()
macro(APPEND_PREFIX PREFIX FILELIST)
    SET(outlist )
    FOREACH(root ${${FILELIST}})
      LIST(APPEND outlist ${PREFIX}/${root})
    ENDFOREACH(root)
    SET(${FILELIST} ${outlist})
endmacro()

IMPORT(common/src COMMON_SRC)
IMPORT(symbol/src SYMBOL_SRC)
IMPORT(order/src ORDER_SRC)
IMPORT(perf/src PERF_SRC)
IMPORT(sopalin/src SOPALIN_SRC)
IMPORT(blend/src BLEND_SRC)
IMPORT(fax/src FAX_SRC)
IMPORT(kass/src KASS_SRC)
IMPORT(sparse-matrix/src SPM_SRC)
#SET(WRAPPER wrapper/petsc/src_mat_impl_aij_mpi_pastix_pastix.c)
    
add_library(pastix STATIC ${COMMON_SRC} ${SYMBOL_SRC} ${ORDER_SRC}
     ${SOPALIN_SRC} ${PERF_SRC} ${BLEND_SRC} ${FAX_SRC} ${KASS_SRC}
     ${SPM_SRC})
target_compile_definitions(pastix PRIVATE
    ${COMMON_DEFINITIONS} ${TYPE_DEFINITIONS}
)
target_link_libraries(pastix PUBLIC
    $<$<BOOL:${USE_PTSCOTCH}>:ptscotch>
    $<$<BOOL:${USE_SCOTCH}>:scotch>
    $<$<BOOL:${USE_METIS}>:metis>
    blas
    ${CMAKE_Fortran_IMPLICIT_LIBRARIES})
target_include_directories(pastix 
    PUBLIC ${PASTIX_INCLUDE_DIR}
    PRIVATE common/src symbol/src order/src blend/src sopalin/src
    perf/src perf/src/num_recipes fax/src kass/src
    matrix_drivers/src
    ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
    )

add_subdirectory(matrix_drivers/src)

if (BUILD_TESTING)
    add_subdirectory(example/src)
endif()