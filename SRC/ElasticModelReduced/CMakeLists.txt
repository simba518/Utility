SET(LIBRARY_NAME utility_elastic_reduced)
AUX_SOURCE_DIRECTORY(. elastic_src)
AUX_SOURCE_DIRECTORY(./cubacode cuba_src)

ADD_LIBRARY(${LIBRARY_NAME} 
  ${elastic_src}
  ${cuba_src}
)

# install 
FILE(GLOB header ${header} *.h)
INSTALL(FILES ${header} DESTINATION ${HEADFILE_DEST})
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_DEST})
