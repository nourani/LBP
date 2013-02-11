# Build Python extension module

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_SWIG_FLAGS "")

include_directories(${CMAKE_SOURCE_DIR})

set_source_files_properties(lbp.i PROPERTIES CPLUSPLUS ON)
swig_add_module(lbp python lbp.i ../LBP.cpp)
swig_link_libraries(lbp 
    ${PYTHON_LIBRARIES} 
    opencv_core opencv_features2d
    fftw3
)
