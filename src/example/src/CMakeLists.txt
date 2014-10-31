# Generate specific headers
# Only BUILD_TYPE=d done so far
execute_process(COMMAND sed 
        -f ${CMAKE_CURRENT_SOURCE_DIR}/sed_${BUILD_TYPE}.in
        ${CMAKE_CURRENT_SOURCE_DIR}/utils.h
    OUTPUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${BUILD_TYPE}utils.h)

    
add_executable(simple simple.c)
target_link_libraries(simple pastix matrix_driver)