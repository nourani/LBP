# LBP
cmake_minimum_required(VERSION 2.8)

# project name
project(LBP)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )

# find dependencies 
find_path(FFTW_INCLUDE_DIR fftw3.h  ${FFTW_INCLUDE_DIRS})
find_library(FFTW_LIBRARY fftw3 ${FFTW_LIBRARY_DIRS})
find_package(OpenCV REQUIRED)

add_library(LBP
    SHARED
	LBP.cpp
)
target_link_libraries(LBP
	${OpenCV_LIBS} ${FFTW_LIBRARY} m
)

add_executable(LBPTest
	LBPTest.cpp
)
target_link_libraries(LBPTest
	LBP
)

add_executable(LBPMain
	LBPMain.cpp
)
target_link_libraries(LBPMain
	LBP
)

option(BUILD_PYTHON "Build SWIG/Python module" OFF)
if(BUILD_PYTHON)
    add_subdirectory(python)
endif(BUILD_PYTHON)


###############################################################################
## 								Install
# Binaries and libs
install(TARGETS LBP LBPMain
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
#	ARCHIVE DESTINATION lib
)
# Header files
install(FILES LBP.hpp
	DESTINATION include 
)
