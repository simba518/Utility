SET(APPNAME RenderTest)

AUX_SOURCE_DIRECTORY(. Test_src)

ADD_EXECUTABLE(${APPNAME} ${Test_src})

LINK_PROJECT(${APPNAME}

  utility_common
  utility_math
  utility_ui
)

TARGET_LINK_LIBRARIES( ${APPNAME} 

  /home/simba/Workspace/Utility/Lib/Release/libutility_ui.a
  ${QGLViewer_LIBRARIES}
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${Boost_LIBRARIES}
  )
