#message("Reading content from file ${FILE}")
FILE(READ ${FILE} CONTENT)
#message("Appending content to file ${TARGET}")
FILE(APPEND ${TARGET} "${CONTENT}")