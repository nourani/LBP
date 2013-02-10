# LBP

cmake_minimum_required(VERSION 2.8)

# project name
PROJECT(NOURANI)

SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}/proj/acfr/nourani")
#SET(CMAKE_CXX_COMPILER "g++-3.4")
#SET(CMAKE_CXX_LINKER "g++-3.4")

INCLUDE_DIRECTORIES(
	./
	/usr/local/include
	/opt/local/include
)

LINK_DIRECTORIES(
	/usr/local/lib
	/opt/local/lib
)

ADD_DEFINITIONS( -ggdb -O3 -D_x86 -D_OSX )



ADD_LIBRARY(LBP
	STATIC
	LBP.cpp
)
TARGET_LINK_LIBRARIES(LBP
	
)


ADD_EXECUTABLE(LBPTest
	LBPTest.cpp
)
TARGET_LINK_LIBRARIES(LBPTest
	LBP opencv_core opencv_imgproc opencv_highgui fftw3 m
)

###############################################################################
## 								Install
# Binaries and libs
#INSTALL(TARGETS CopyFiles
#	RUNTIME DESTINATION bin
#	LIBRARY DESTINATION lib
#	ARCHIVE DESTINATION lib
#)
# Header files
#INSTALL(FILES 
#	DESTINATION include/apps
#)
