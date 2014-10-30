add_executable(genheader genheader.c)
target_include_directories(genheader PRIVATE
    ../../common/src)

SET(PASTIXH ${PASTIX_INCLUDE_DIR}/pastix.h)
SET(PASTIXH_VERSION ${PASTIX_INCLUDE_DIR}/pastix${VERSIONTYPE}.h)
SET(PASTIXFH ${PASTIX_INCLUDE_DIR}/pastix_fortran.h)
SET(PASTIXFH_VERSION ${PASTIX_INCLUDE_DIR}/pastix_fortran${VERSIONTYPE}.h)
SET(MURGEINC ${PASTIX_INCLUDE_DIR}/murge.inc)
SET(MURGEINC_VERSION ${PASTIX_INCLUDE_DIR}/murge${VERSIONTYPE}.inc)

add_custom_target(headers
    # C
    COMMAND genheader ${PASTIXH_VERSION} ${PASTIXFH_VERSION} ${MURGEINC_VERSION}
        ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
        ${PASTIX_SRC}/murge/include/murge.h
        ${PASTIX_SRC}/murge/scripts/genfort.pl
        ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh C
    COMMAND ${CMAKE_COMMAND} 
        -DFILE=${PASTIX_SRC}/sopalin/src/pastix.h
        -DTARGET=${PASTIXH_VERSION}
        -P ${PASTIX_SOURCE_DIR}/Append.cmake 
    COMMAND ${CMAKE_COMMAND} 
        -DFILE=${PASTIX_SRC}/common/src/api.h
        -DTARGET=${PASTIXH_VERSION}
        -P ${PASTIX_SOURCE_DIR}/Append.cmake
    COMMAND ${CMAKE_COMMAND} -E copy ${PASTIXH_VERSION} ${PASTIXH}
    # Fortran
    COMMAND genheader ${PASTIXH_VERSION} ${PASTIXFH_VERSION} ${MURGEINC_VERSION}
        ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
        ${PASTIX_SRC}/murge/include/murge.h
        ${PASTIX_SRC}/murge/scripts/genfort.pl
        ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh Fortran
    COMMAND ${CMAKE_COMMAND} 
        -DFILE=${PASTIX_SRC}/common/src/pastix_fortran.h
        -DTARGET=${PASTIXFH_VERSION}
        -P ${PASTIX_SOURCE_DIR}/Append.cmake 
    COMMAND ${CMAKE_COMMAND} -E copy ${PASTIXFH_VERSION} ${PASTIXFH}
    # Murge
    COMMAND genheader ${PASTIXH_VERSION} ${PASTIXFH_VERSION} ${MURGEINC_VERSION}
        ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
        ${PASTIX_SRC}/murge/include/murge.h
        ${PASTIX_SRC}/murge/scripts/genfort.pl
        ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh Murge
    COMMAND ${CMAKE_COMMAND} 
        -DFILE=${PASTIX_SRC}/common/src/pastix_fortran.inc
        -DTARGET=${MURGEINC_VERSION}
        -P ${PASTIX_SOURCE_DIR}/Append.cmake
    COMMAND ${CMAKE_COMMAND} -E copy ${MURGEINC_VERSION} ${MURGEINC}
)