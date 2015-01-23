set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( SiWiR_ex4 CXX )

SET( SRCS

   heat.cpp
   # more *.cpp files follow here
)

OPTION( ARCHITECTURE_OPTIMIZATION   "enables/disables platform specific optimizations (on by default)" ON )
OPTION( ENABLE_WERROR   "enables/disables treat warnings as errors (on by default)" ON )

IF( CMAKE_COMPILER_IS_GNUCXX )
	ADD_DEFINITIONS( -Wall -Wextra -Wshadow -std=c++0x )
	IF( CMAKE_COMPILER_IS_GNUCC )
		ADD_DEFINITIONS( -march=native )
	ENDIF ()	
	IF( ENABLE_WERROR )
		ADD_DEFINITIONS( -Werror )
	ENDIF ()
ENDIF()

ADD_EXECUTABLE( heat ${SRCS} )
