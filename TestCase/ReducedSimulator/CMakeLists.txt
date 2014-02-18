SET(APPNAME ReducedSimulator)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_common.a
  ${LIBDIR}/libutility_math.a
  ${LIBDIR}/libutility_elastic.a
  ${LIBDIR}/libutility_elastic_reduced.a
  ${LIBDIR}/libutility_collision.a
  /usr/lib/x86_64-linux-gnu/libgomp.so.1
  ${Boost_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )
