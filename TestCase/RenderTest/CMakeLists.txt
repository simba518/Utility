SET(APPNAME RenderTest)

AUX_SOURCE_DIRECTORY(. Test_src)

ADD_EXECUTABLE(${APPNAME} ${Test_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_ui.a
  ${LIBDIR}/libutility_common.a
  ${LIBDIR}/libutility_math.a
  ${QGLViewer_LIBRARIES}
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${Boost_LIBRARIES}
  )
