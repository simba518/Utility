SET(LIBRARY_NAME math)
AUX_SOURCE_DIRECTORY(. math_src)

ADD_LIBRARY(${LIBRARY_NAME} 
  ${math_src}
)

# install 
FILE(GLOB header ${header} *.h)
INSTALL(FILES ${header} DESTINATION /usr/local/include/utility)
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION /usr/local/lib/)
