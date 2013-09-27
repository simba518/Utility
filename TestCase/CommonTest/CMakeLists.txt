SET(APPNAME CommonTest)

AUX_SOURCE_DIRECTORY(. CommonTest_src)

ADD_EXECUTABLE(${APPNAME} ${CommonTest_src})

LINK_PROJECT(${APPNAME}
  utility_common
)

TARGET_LINK_LIBRARIES( ${APPNAME} 
  
  ${Boost_LIBRARIES}
  )
