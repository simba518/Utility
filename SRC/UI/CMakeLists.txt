SET(LIBRARY_NAME utility_ui)

SET(NONE_QT_CPP

  Render/SelfRenderRect.cpp
  Render/FileDialog.cpp
  Manipulation/SelectionGroup.cpp
)

SET(UI_SOURCES

  Manipulation/SelectCtrl.cpp
  Manipulation/DragCtrl.cpp
  Animation/AniCtrl.cpp
  Render/QGLViewerExt.cpp
)

SET(UI_HEADERS 

  Manipulation/SelectCtrl.h
  Manipulation/DragCtrl.h
  Animation/AniCtrl.h
  Animation/AniSliderBar.h
  Render/QGLViewerExt.h
)

SET(UI_FORMS 

)

QT4_WRAP_CPP(UI_HEADERS_MOC ${UI_HEADERS})
QT4_WRAP_UI(UI_FORMS_HEADERS ${UI_FORMS})

ADD_LIBRARY(${LIBRARY_NAME}

  ${NONE_QT_CPP}
  ${UI_SOURCES}
  ${UI_HEADERS_MOC}
  ${UI_FORMS_HEADERS}
)

# install
FILE(GLOB header ${header} ./Animation/*.h ./Manipulation/*.h ./Render/*.h CACHE INTERNAL)
INSTALL(FILES ${header} DESTINATION ${HEADFILE_DEST})
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_DEST})
