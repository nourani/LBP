# Build Python extension module

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Python headers
find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

# Python executable


# Python site packages
# http://stackoverflow.com/a/14796658
execute_process(
    COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" 
    OUTPUT_VARIABLE PYTHON_DIST_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Find NumPy by asking python where it is
find_package(PythonInterp 2.7 REQUIRED)
exec_program ("${PYTHON_EXECUTABLE}"
    ARGS "-c 'import numpy; print numpy.get_include()'"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    RETURN_VALUE NUMPY_NOT_FOUND)
set (NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE PATH "Numpy include path")
include_directories(${NUMPY_INCLUDE_DIR})

set(CMAKE_SWIG_FLAGS "")

include_directories(${CMAKE_SOURCE_DIR})

set_source_files_properties(lbp.i PROPERTIES CPLUSPLUS ON)
swig_add_module(lbp python lbp.i ../LBP.cpp)
swig_link_libraries(lbp 
    ${PYTHON_LIBRARIES} 
    ${OpenCV_LIBS} 
    ${FFTW_LIBRARY}
)



install(
    TARGETS _lbp DESTINATION ${PYTHON_DIST_PACKAGES}
)

install(
    FILES ${CMAKE_BINARY_DIR}/python/lbp.py DESTINATION ${PYTHON_DIST_PACKAGES}
)
