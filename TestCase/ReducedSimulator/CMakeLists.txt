SET(APPNAME ReducedSimulator)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_common.a
  ${LIBDIR}/libutility_math.a
  ${LIBDIR}/libutility_elastic.a
  ${LIBDIR}/libutility_elastic_reduced.a
  ${LIBDIR}/libutility_collision.a
  ${Boost_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )
