if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/api_str_to_int.c)
    message(STATUS "Matrix drivers: Generating api_str_to_int source file")
    #add_custom_command(
    #    OUTPUT GREPOUT
    #    COMMAND grep "[ID]PARM_[A-Z_]*[ ]*=[ ]*[0-9]*" ${CMAKE_CURRENT_SOURCE_DIR}/../../common/src/api.h > GREPOUT
    #)
    execute_process(
        COMMAND grep "[ID]PARM_[A-Z_]*[ ]*=[ ]*[0-9]*" ${CMAKE_CURRENT_SOURCE_DIR}/../../common/src/api.h
        OUTPUT_FILE GREPOUT
    )    
    execute_process(
        COMMAND sed -f ${CMAKE_CURRENT_SOURCE_DIR}/sed.in
        INPUT_FILE GREPOUT
        OUTPUT_VARIABLE STRTOINTVALUES
    )
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/api_str_to_int.c.in
     ${CMAKE_CURRENT_BINARY_DIR}/api_str_to_int.c @ONLY)     
endif()

SET(MD_SRC
    skitf.f		
    iohb.c		
    mmio.c		
    common_drivers.c	
    get_options.c
    ${CMAKE_CURRENT_BINARY_DIR}/api_str_to_int.c
)
SET(MD_SRC_MA
    read_matrix.c
    rsaread.c		
    hbread.c		
    mmread.c		
    mmdread.c		
    petscread.c		
    cccread.c		
    olafread.c		
    chbread.c		
    cscdread.c		
    peerread.c		
    threefilesread.c	
    fdupread.c		
    laplacian.c
)

SET(OBJS $<TARGET_OBJECTS:matrix_driver_common>)
foreach(prec s d c z)
    if(BUILD_PRECISION MATCHES [${prec}])
        SETINCDIR(${prec})
        
        SET(PASTIXH ${PASTIX_INCLUDE_DIR_PREC}/pastix.h)
        SET(PASTIXFH ${PASTIX_INCLUDE_DIR_PREC}/pastix_fortran.h)
        SET(MURGEINC ${PASTIX_INCLUDE_DIR_PREC}/murge.inc)
        SET(CSCDUTIL ${PASTIX_INCLUDE_DIR_PREC}/cscd_utils.h)
        add_custom_command(OUTPUT ${PASTIXH} ${PASTIXFH} ${CSCDUTIL}
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
            # Misc headers
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/read_matrix.h ${PASTIX_INCLUDE_DIR_PREC}
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/get_options.h ${PASTIX_INCLUDE_DIR_PREC}
        )
        if (NOT MINGW)
            # Somehow for MSYS environments something breaks within the rubbish-code to generate headers for murge.
            # As this is not used within OpenCMISS thus far, we simply wont create the headers on windows. 
            add_custom_command(#${prec}headers_murge
                OUTPUT ${MURGEINC}
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
            )
        endif()
        
        add_library(${prec}matrix_driver ${MD_SRC} ${MD_SRC_MA} ${PASTIX_INCLUDE_DIR_PREC}/pastix.h)
        target_include_directories(${prec}matrix_driver
            PRIVATE
            ${PASTIX_INCLUDE_DIR_PREC}
            ${CMAKE_CURRENT_SOURCE_DIR}
            PUBLIC	
            ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
        )
        SETTYPEDEFS(${prec})
        target_compile_definitions(${prec}matrix_driver
            PRIVATE ${COMMON_DEFINITIONS} ${TYPEDEFS})
        
        set_target_properties(${prec}matrix_driver PROPERTIES OUTPUT_NAME pastix-${prec}matrix-driver-${PASTIX_VERSION}) 
        target_link_libraries(${prec}matrix_driver PUBLIC pastix)
        SETINCDIR(${NEUTRAL_PRECISION})
        install(TARGETS ${prec}matrix_driver EXPORT pastix-config
            DESTINATION lib
            COMPONENT ${prec}drivers
            INCLUDES DESTINATION ${PASTIX_INCDIR_SUFFIX_PREC})
    endif()
endforeach()
