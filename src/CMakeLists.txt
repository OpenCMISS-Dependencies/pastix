add_subdirectory(utils/src)

macro(IMPORT DIR SRC)
    add_subdirectory(${DIR})
    get_directory_property(${SRC}
        DIRECTORY ${DIR}
        DEFINITION ${SRC})
    SET(outlist )
    FOREACH(root ${${SRC}})
      LIST(APPEND outlist ${DIR}/${root})
    ENDFOREACH(root)
    SET(${SRC} ${outlist})
endmacro()

IMPORT(common/src COMMON_SRC)
IMPORT(symbol/src SYMBOL_SRC)
IMPORT(order/src ORDER_SRC)
IMPORT(perf/src PERF_SRC)
IMPORT(sopalin/src SOPALIN_SRC)
IMPORT(blend/src BLEND_SRC)
IMPORT(fax/src FAX_SRC)
IMPORT(kass/src KASS_SRC)
# still need: wrapper sparse-matrix
    
add_library(pastix STATIC ${COMMON_SRC} ${SYMBOL_SRC} ${ORDER_SRC}
     ${SOPALIN_SRC} ${PERF_SRC} ${BLEND_SRC} ${FAX_SRC} ${KASS_SRC})
target_link_libraries(pastix PUBLIC
    $<$<BOOL:${USE_PTSCOTCH}>:ptscotch>
    $<$<BOOL:${USE_SCOTCH}>:scotch>
    blas
    ${CMAKE_Fortran_IMPLICIT_LIBRARIES})
target_include_directories(pastix PRIVATE
    common/src symbol/src order/src blend/src sopalin/src
    perf/src perf/src/num_recipes fax/src kass/src
    matrix_drivers/src
    ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})

add_subdirectory(matrix_drivers/src)