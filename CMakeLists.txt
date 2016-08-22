# LBP
cmake_minimum_required(VERSION 2.8)

# project name
project(LBP)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )

# find dependencies 
find_path(FFTW_INCLUDE_DIR fftw3.h  ${FFTW_INCLUDE_DIRS})
find_library(FFTW_LIBRARY fftw3 ${FFTW_LIBRARY_DIRS})
find_package(OpenCV REQUIRED)


INCLUDE_DIRECTORIES(
        ./
        /usr/local/include
        /opt/local/include
)

LINK_DIRECTORIES(
	/usr/local/cuda/lib/
	/opt/local/lib
)

add_library(LBP
    SHARED
	LBP.cpp 
)
target_link_libraries(LBP
	${OpenCV_LIBS} ${FFTW_LIBRARY}  m
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


option(BUILD_LBPGPU "Build GPU module" OFF)
if(BUILD_LBPGPU)
	find_package(CUDA REQUIRED)

	IF(APPLE)
		SET(CUDA_HOST_COMPILER /usr/bin/clang CACHE FILEPATH "Setting clang as the CUDA compiler" FORCE)
		SET(CUDA_NVCC_FLAGS "-Xcompiler -stdlib=libstdc++; -Xlinker -stdlib=libstdc++; -arch=sm_20" CACHE STRING "Setting NVCC compiler flags" FORCE)
	ENDIF(APPLE)

	CUDA_ADD_EXECUTABLE( LBPGPU
		LBPGPU.cu
	)
	TARGET_LINK_LIBRARIES(LBPGPU
		${CUDA_LIBRARIES} nppc nppi X11
	)

endif(BUILD_LBPGPU)



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
if(BUILD_LBPGPU)
install(TARGETS LBPGPU
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
#	ARCHIVE DESTINATION lib
)
endif(BUILD_LBPGPU)
# Header files
install(FILES LBP.hpp
	DESTINATION include 
)
