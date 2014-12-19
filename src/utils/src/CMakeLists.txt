# Need different gen_header binary for each precision
# (dont ask me...)
foreach(prec s d c z)
    if(BUILD_PRECISION MATCHES [${prec}])
        if (NEED_HEADERS_${prec})
            add_executable(${prec}genheader genheader.c)
            SET(TYPEDEFS )
            if (${prec} MATCHES [dz])
                LIST(APPEND TYPEDEFS PREC_DOUBLE)     
            endif()
            if (${prec} MATCHES [cz])
                LIST(APPEND TYPEDEFS TYPE_COMPLEX)     
            endif()
            target_compile_definitions(${prec}genheader
                PRIVATE ${COMMON_DEFINITIONS} ${TYPEDEFS})
            target_include_directories(${prec}genheader PRIVATE
                ../../common/src)
    
            SETINCDIR(${prec})
            SET(PASTIXH ${PASTIX_INCLUDE_DIR_PREC}/pastix.h)
            SET(PASTIXFH ${PASTIX_INCLUDE_DIR_PREC}/pastix_fortran.h)
            SET(MURGEINC ${PASTIX_INCLUDE_DIR_PREC}/murge.inc)
            SET(CSCDUTIL ${PASTIX_INCLUDE_DIR_PREC}/cscd_utils.h)
            add_custom_target(${prec}headers
            # C
            COMMAND ${prec}genheader ${PASTIXH} ${PASTIXFH} ${MURGEINC}
                ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
                ${PASTIX_SRC}/murge/include/murge.h
                ${PASTIX_SRC}/murge/scripts/genfort.pl
                ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh C
            COMMAND ${CMAKE_COMMAND} 
                -DFILE=${PASTIX_SRC}/sopalin/src/pastix.h
                -DTARGET=${PASTIXH}
                -P ${PASTIX_SOURCE_DIR}/Append.cmake 
            COMMAND ${CMAKE_COMMAND} 
                -DFILE=${PASTIX_SRC}/common/src/api.h
                -DTARGET=${PASTIXH}
                -P ${PASTIX_SOURCE_DIR}/Append.cmake
            # Fortran
            COMMAND ${prec}genheader ${PASTIXH} ${PASTIXFH} ${MURGEINC}
                ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
                ${PASTIX_SRC}/murge/include/murge.h
                ${PASTIX_SRC}/murge/scripts/genfort.pl
                ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh Fortran
            COMMAND ${CMAKE_COMMAND} 
                -DFILE=${PASTIX_SRC}/common/src/pastix_fortran.h
                -DTARGET=${PASTIXFH}
                -P ${PASTIX_SOURCE_DIR}/Append.cmake 
            # Murge
            COMMAND ${prec}genheader ${PASTIXH} ${PASTIXFH} ${MURGEINC}
                ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
                ${PASTIX_SRC}/murge/include/murge.h
                ${PASTIX_SRC}/murge/scripts/genfort.pl
                ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh Murge
            COMMAND ${CMAKE_COMMAND} 
                -DFILE=${PASTIX_SRC}/common/src/pastix_fortran.inc
                -DTARGET=${MURGEINC}
                -P ${PASTIX_SOURCE_DIR}/Append.cmake
            # CSCD Utils
            COMMAND ${prec}genheader ${CSCDUTIL} ${PASTIXFH} ${MURGEINC}
                ${PASTIX_SRC}/sopalin/src/murge_pastix.inc
                ${PASTIX_SRC}/murge/include/murge.h
                ${PASTIX_SRC}/murge/scripts/genfort.pl
                ${PASTIX_SRC}/utils/src/insert-fortran-inc.sh C
            COMMAND ${CMAKE_COMMAND} 
                -DFILE=${PASTIX_SRC}/sopalin/src/cscd_utils.h
                -DTARGET=${CSCDUTIL}
                -P ${PASTIX_SOURCE_DIR}/Append.cmake 
            )
            # Ensure the genheader targets are run after another (write race conditions)
            foreach(doneprec ${DONE})
                add_dependencies(${prec}genheader ${doneprec}genheader)
            endforeach()
            LIST(APPEND DONE ${prec})
        endif()
    endif()
endforeach()