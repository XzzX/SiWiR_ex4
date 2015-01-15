CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( hello )

SET( SRCS

   hello.cpp
   # more *.cpp files follow here
)

OPTION( BIG_HELLO   "enables/disables big letter mode (off by default)"                                          OFF )
OPTION( ENABLE_ANSI "enables/disables addition of '-ansi -pedantic' flags for gnu c++ compiler (off by default)" OFF )

IF( BIG_HELLO )
   ADD_DEFINITIONS( -DBIG_HELLO )
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCC )
   IF( ENABLE_ANSI )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pedantic" )
   ENDIF()
ENDIF()

ADD_EXECUTABLE( hello ${SRCS} )
