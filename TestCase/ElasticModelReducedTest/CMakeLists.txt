SET(APPNAME ElasticModelReducedTest)

AUX_SOURCE_DIRECTORY(. ElasticModel_src)

ADD_EXECUTABLE(${APPNAME} ${ElasticModel_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 
  ${LIBDIR}/libutility_elastic.a
  ${LIBDIR}/libutility_elastic_reduced.a
  ${LIBDIR}/libutility_common.a
  ${Boost_LIBRARIES}
  )
