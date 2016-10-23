#include <Python.h>
#include <stdio.h>
#include <math.h>


double _c_qvc(double T, double omega)
{
    const double kb1 = 1.43877516016792017517e+00;
    const double K_BOLTZMANN_RY = 6.33363068565903064295e-06;
    double x, expx;

    if (T<1E-9 || omega<1E-9) {
        return 0.0;
    }
    x = omega * kb1 / T; 
    expx = exp(-x);    // exponential term
    if (expx>1E-3) {         // compute normally
        return x*x*K_BOLTZMANN_RY*expx/pow(expx-1.0,2);
    }
    else                   // Taylor series
        return K_BOLTZMANN_RY*expx*x/pow(x-0.5*x*x+0.16666666666666667*x*x*x+0.04166666666666666667*x*x*x*x,2);
}


static PyObject* grunc_c_qvc(PyObject* self, PyObject* args)
{
    double T,omega;

    if (!PyArg_ParseTuple(args, "dd", &T, &omega))
        return NULL;

    return Py_BuildValue("d", _c_qvc(T,omega));
}


//Method definition object for this extension, these argumens mean:
//ml_name: The name of the method
//ml_meth: Function pointer to the method implementation
//ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
//ml_doc:  Contents of this method's docstring
static PyMethodDef grunc_methods[] = { 
    {   
        "c_qvc",
        grunc_c_qvc,
        METH_VARARGS,
        "Calculate the c_qv terms in a C extension."
    },  
    {NULL, NULL, 0, NULL}   // Sentinel
};

//Module definition
//The arguments of this structure tell Python what to call your extension,
//what it's methods are and where to look for it's method definitions
static struct PyModuleDef grunc_definition = { 
    PyModuleDef_HEAD_INIT,
    "grunc",
    "A Python module from C code.",
    -1, 
    grunc_methods
};

//Module initialization
//Python calls this function when importing your extension. It is important
//that this function is named PyInit_[[your_module_name]] exactly, and matches
//the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_grunc(void)
{
    Py_Initialize();

    return PyModule_Create(&grunc_definition);
}
