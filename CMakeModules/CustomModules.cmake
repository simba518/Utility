MACRO( INCLUDE_FOR_PROJECT PROJNAME )
  FOREACH( varname ${ARGN} )
    IF( ${varname}_INCLUDE_DIR )
      INCLUDE_DIRECTORIES( "${${varname}_INCLUDE_DIR}" )
    ELSE( ${varname}_INCLUDE_DIR )
      INCLUDE_DIRECTORIES( "${varname}" )
    ENDIF( ${varname}_INCLUDE_DIR )
  ENDFOREACH( varname )
ENDMACRO( INCLUDE_FOR_PROJECT PROJNAME )




MACRO( LINK_PROJECT PROJNAME )
  FOREACH( varname ${ARGN} )
    IF( ${varname}_LIBRARY )
      TARGET_LINK_LIBRARIES( ${PROJNAME} "${${varname}_LIBRARY}" )
    ELSE( ${varname}_LIBRARY )
	  IF( varname )
        TARGET_LINK_LIBRARIES( ${PROJNAME} "${${varname}}" )
	  ELSE( varname )
		MESSAGE("error: can't find library: " ${varname})
      ENDIF( varname )
    ENDIF( ${varname}_LIBRARY )
  ENDFOREACH( varname )
ENDMACRO( LINK_PROJECT PROJNAME )





MACRO(FIND_INTERNAL_PACKAGE INTERNAL_LIB_PATH)

  MESSAGE("finding libraries at " ${INTERNAL_LIB_PATH})
  SET(${varname}_LIBRARY)

  FOREACH(varname ${ARGN})
	FIND_LIBRARY( ${varname}_LIBRARY 
	  NAMES lib${varname}.a
	  PATHS
	  ${INTERNAL_LIB_PATH}
	  )

	IF( !${${varname}_LIBRARY} )
	  FIND_LIBRARY( ${varname}_LIBRARY 
		NAMES lib${varname}.so
		PATHS
		${INTERNAL_LIB_PATH}
		)
	ENDIF( !${${varname}_LIBRARY} )

	IF( !${${varname}_LIBRARY} )
	  MESSAGE("error: unfound library" ${varname})
	ELSE( !${${varname}_LIBRARY} )
	  MESSAGE( "find " ${${varname}_LIBRARY} )
	ENDIF( !${${varname}_LIBRARY} )

  ENDFOREACH(varname)

ENDMACRO(FIND_INTERNAL_PACKAGE INTERNAL_LIB_PATH)