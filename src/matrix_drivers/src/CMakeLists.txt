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
 ${CMAKE_CURRENT_SOURCE_DIR}/api_str_to_int.c @ONLY)

SET(MD_SRC
    skitf.f		
    iohb.c		
    mmio.c		
    common_drivers.c	
    get_options.c
    api_str_to_int.c
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
        add_library(${prec}matrix_driver ${MD_SRC} ${MD_SRC_MA})
        #add_library(${prec}matrix_driver OBJECT ${MD_SRC_MA})
        if (NEED_HEADERS_${prec})
            add_dependencies(${prec}matrix_driver ${prec}headers)
        endif()
        target_include_directories(${prec}matrix_driver PRIVATE
            ${PASTIX_INCLUDE_DIR_PREC}
            ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
        SETTYPEDEFS(${prec})
        target_compile_definitions(${prec}matrix_driver
            PRIVATE ${COMMON_DEFINITIONS} ${TYPEDEFS})
           
        #LIST(APPEND OBJS "$<TARGET_OBJECTS:${prec}matrix_driver>")
        configure_file(read_matrix.h ${PASTIX_INCLUDE_DIR_PREC} COPYONLY)
        configure_file(get_options.h ${PASTIX_INCLUDE_DIR_PREC} COPYONLY)
        
        set_target_properties(${prec}matrix_driver PROPERTIES OUTPUT_NAME pastix-${prec}matrix-driver-${PASTIX_VERSION}) 
        target_link_libraries(${prec}matrix_driver PUBLIC pastix)
        SETINCDIR(${NEUTRAL_PRECISION})
        install(TARGETS ${prec}matrix_driver EXPORT pastix-config
            DESTINATION lib
            COMPONENT ${prec}drivers
            INCLUDES DESTINATION ${PASTIX_INCDIR_SUFFIX_PREC})
    endif()
endforeach()