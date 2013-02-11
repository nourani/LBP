%module lbp

%{ 
    #define SWIG_FILE_WITH_INIT
    #include "LBP.hpp"
%}

%include "numpy.i"
%include exception.i

%init %{
    import_array();
%}

/*
 * Convert numpy array to opencv
 */
%typemap(typecheck, fragment="NumPy_Object_to_Array,NumPy_Array_Requirements") 
Mat 
{
   $1 = is_array($input) ? 1 : 0;
}
%typemap(in, fragment="NumPy_Object_to_Array,NumPy_Array_Requirements") 
Mat 
{ 
    PyArrayObject* ary = obj_to_array_no_conversion($input, NPY_FLOAT64);
    if(!ary || !require_contiguous(ary) || !require_dimensions(ary, 2)) 
        SWIG_fail;

    unsigned char *data = (unsigned char *)array_data(ary); 
    npy_intp* dims = array_dimensions(ary);
    Mat mat(dims[0], dims[1], CV_64FC1, data);
    $1 = mat;
}


/*
 * Convert std::vector<double> to numpy arrays 
 */
%typemap(out) vector<double> {
    npy_intp dims[1];
    dims[0] = $1.size(); 
    
    /* Make copy of the data */
    double *data = (double *)malloc(dims[0] * sizeof(double));
    memcpy(data, &$1[0], dims[0] * sizeof(double));   
    PyArrayObject* py_result = (PyArrayObject* )PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data);
    py_result->flags |= NPY_OWNDATA;
    $result = (PyObject*)py_result;
}


%include "LBP.hpp"
