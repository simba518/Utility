SET(APPNAME tempt)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 
  ${LIBDIR}/libutility_common.a
  ${Boost_LIBRARIES}
  )
