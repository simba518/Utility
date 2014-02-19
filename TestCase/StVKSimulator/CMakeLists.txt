SET(APPNAME StVKSimulator)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

FIND_PACKAGE(OpenMP REQUIRED)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_common.a
  ${LIBDIR}/libutility_math.a
  ${LIBDIR}/libutility_elastic.a
  ${LIBDIR}/libutility_collision.a
  ${Boost_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )