# Need different gen_header binary for each precision
foreach(prec s d c z)
    if(BUILD_PRECISION MATCHES [${prec}])
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

        # Ensure the genheader targets are run after another (write race conditions)
        foreach(doneprec ${DONE})
            add_dependencies(${prec}genheader ${doneprec}genheader)
        endforeach()
        LIST(APPEND DONE ${prec})
    endif()
endforeach()