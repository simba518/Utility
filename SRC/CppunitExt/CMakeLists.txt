SET(LIBRARY_NAME cppunitext)
AUX_SOURCE_DIRECTORY(. basetestor_src)

ADD_LIBRARY(${LIBRARY_NAME} 
  ${basetestor_src}
)

# install 
FILE(GLOB header ${header} *.h)
INSTALL(FILES ${header} DESTINATION /usr/local/include/utility)
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION /usr/local/lib/)
