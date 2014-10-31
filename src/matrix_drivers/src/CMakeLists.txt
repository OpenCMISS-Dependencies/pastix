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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/get_options.h
    ${PASTIX_INCLUDE_DIR}/get_options.h COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/read_matrix.h
    ${PASTIX_INCLUDE_DIR}/read_matrix.h COPYONLY)    

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

add_library(matrix_driver STATIC ${MD_SRC} ${MD_SRC_MA})
add_dependencies(matrix_driver headers)
target_link_libraries(matrix_driver PUBLIC pastix)
target_include_directories(matrix_driver PRIVATE
    ${PASTIX_INCLUDE_DIR}
    ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})