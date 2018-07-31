/*
 * Follow the paper A&A 367, 809-825 (2001)
 * The multifrequency emission of Mrk 501 from radio to TeV gamma-rays
 * The emission and absorption coefficient calculation using the approximation (Appendix A)
 * The equation (3) has some difficult in high energy, donot use it.
 */
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "Synchrotron.h"
#include "PhysicalConstant.h"


double U_B(double B){
    return B * B / (8 * M_PI);
}

double Nu_B(double B){
    return e * B / (2 * M_PI * m_e * c);
}

double f_(double x){
    int N = 10000;
    double dx = (70*x - x) / N;
    double sum = 0;
    for(int i=0; i<N; ++i){
        double xi = x + i*dx;
        double kv = gsl_sf_bessel_Knu(5./3., xi);
        sum += kv * dx;
    }
    return x * sum;
}


static PyObject* F(PyObject* self, PyObject* args){
    PyArrayObject* x;
    PyArg_ParseTuple(args, "O!", &PyArray_Type, &x);
    
    PyArrayObject* in_array = x;
    PyObject* out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = f_(*in_dataptr);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}

//Follow the paper ApJ 167, 26 (2006) Models for nonthermal photon spectra
double W_(double x, double a, double b){
    printf("%f\t%f\t%f\n", 1/2.+b-a, 1+2*b, x);
    double result = gsl_sf_hyperg_U(a, b, x);//???????????????????????????
    printf("%f\n", result);
    return exp(-x/2.) * pow(x, b + 1/2.) * gsl_sf_hyperg_U(1/2.+b-a, 1+2*b, x);
}//Why it doesn't work????????????????????????????????????????????????????

double R_(double x){
    return M_PI / 2. * x * (W_(x, 0, 4./3.) * W_(x, 0, 1./3.) - W_(x, 1./2., 5./6.) * W_(x, -1./2., 5./6.));
}

static PyObject* R(PyObject* self, PyObject* args){
    PyArrayObject* x;
    PyArg_ParseTuple(args, "O!", &PyArray_Type, &x);
    
    PyArrayObject* in_array = x;
    PyObject* out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = R_(*in_dataptr);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


double c1 = 0.78;
double c2 = 0.25;
double c3 = 2.175;

double t(double nu, double gamma, double B){
    double nu_B = Nu_B(B);
    return nu / (3 * pow(gamma, 2) * nu_B);
}

double g_(double nu, double t, double n, double K, double nu_B){
    double a = c2 + (n-1)/2;
    double x = c3*t;
    double fx = gsl_sf_gamma_inc(a, x);

    return K * pow((nu / (3 * nu_B)), -n/2) * pow(t, c2 + (n-1)/2) * pow(c3*t, (1-n)/2 - c2) * fx;
}


//Power-Law
double j_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tmax = t(nu, gamma_max, B);
    
    return 9 * sigma_T * c * u_B * c1 / (24 * pow(M_PI, 2) * nu_B) * sqrt(nu/nu_B)\
             * (g_(nu, tmax, n, K, nu_B) - g_(nu, tmin, n, K, nu_B));
}

double J_PowerLaw(double nu, double B, double* pars){
    double K = pars[0];
    double n = pars[1];
    double gamma_min = pars[2];
    double gamma_max = pars[3];
    return j_PowerLaw(nu, B, K, n, gamma_min, gamma_max);
}


//Break Power-Law
double j_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tb = t(nu, gamma_b, B);
    double tmax = t(nu, gamma_max, B);
    
    return 9 * sigma_T * c * u_B * c1 / (24 * pow(M_PI, 2) * nu_B) * sqrt(nu/nu_B)\
             * (g_(nu, tb, n1, K1, nu_B) - g_(nu, tmin, n1, K1, nu_B) + g_(nu, tmax, n2, K2, nu_B) - g_(nu, tb, n2, K2, nu_B));
}

double J_BrokenPowerLaw(double nu, double B, double* pars){
    double K1 = pars[0];
    double n1 = pars[1];
    double n2 = pars[2];
    double gamma_min = pars[3];
    double gamma_break = pars[4];
    double gamma_max = pars[5];
    double K2 = K1 * pow(gamma_break, n2-n1);
    return j_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_break, gamma_max);
}

double J_Sychrotron(double nu, double B, int SpectrumType, double* pars){
    if(SpectrumType == 1){
        double K = pars[0];
        double n = pars[1];
        double gamma_min = pars[2];
        double gamma_max = pars[3];
        return j_PowerLaw(nu, B, K, n, gamma_min, gamma_max);
    }
    if(SpectrumType == 2){
        double K1 = pars[0];
        double n1 = pars[1];
        double n2 = pars[2];
        double gamma_min = pars[3];
        double gamma_break = pars[4];
        double gamma_max = pars[5];
        double K2 = K1 * pow(gamma_break, n2-n1);
        return j_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_break, gamma_max);
    }
}

//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
static PyObject* J(PyObject* self, PyObject* args){
    PyArrayObject* nu;
    double B;
    int SpectrumType;
    PyTupleObject* SpectrumPars;
    PyArg_ParseTuple(args, "O!diO!", &PyArray_Type, &nu, &B, &SpectrumType, &PyTuple_Type, &SpectrumPars);
    
    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ++ipar)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* in_array = nu;
    PyObject* out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = J_Sychrotron(*in_dataptr, B, SpectrumType, pars);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


double h_(double nu, double t, double n, double K, double nu_B){
    double a = c2 + n/2;
    double x = c3*t;
    double fx = gsl_sf_gamma_inc(a, x);
    return K * (n+2) * pow(3, n/2) * pow(nu/nu_B, -n/2) * pow(t, c2+n/2) * pow(c3*t, -c2-n/2) * fx;
}


//Power-Law
double k_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tmax = t(nu, gamma_max, B);
    
    return 3 * sqrt(3) * sigma_T * c * u_B * c1 / (16 * pow(M_PI, 2) * m_e * pow(nu, 2) * nu_B)\
             * (h_(nu, tmax, n, K, nu_B) - h_(nu, tmin, n, K, nu_B));
}

double K_PowerLaw(double nu, double B, double* pars){
    double K = pars[0];
    double n = pars[1];
    double gamma_min = pars[2];
    double gamma_max = pars[3];
    return k_PowerLaw(nu, B, K, n, gamma_min, gamma_max);
}


//Broken Power-Law
double k_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tb = t(nu, gamma_b, B);
    double tmax = t(nu, gamma_max, B);
    
    return 3 * sqrt(3) * sigma_T * c * u_B * c1 / (16 * pow(M_PI, 2) * m_e * pow(nu, 2) * nu_B)\
             * (h_(nu, tb, n1, K1, nu_B) - h_(nu, tmin, n1, K1, nu_B) + h_(nu, tmax, n2, K2, nu_B) - h_(nu, tb, n2, K2, nu_B));
}

double K_BrokenPowerLaw(double nu, double B, double* pars){
    double K1 = pars[0];
    double n1 = pars[1];
    double n2 = pars[2];
    double gamma_min = pars[3];
    double gamma_break = pars[4];
    double gamma_max = pars[5];
    double K2 = K1 * pow(gamma_break, n2-n1);
    return k_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_break, gamma_max);
}

double K_Sychrotron(double nu, double B, int SpectrumType, double* pars){
    if(SpectrumType == 1){
        double K = pars[0];
        double n = pars[1];
        double gamma_min = pars[2];
        double gamma_max = pars[3];
        return k_PowerLaw(nu, B, K, n, gamma_min, gamma_max);
    }
    if(SpectrumType == 2){
        double K1 = pars[0];
        double n1 = pars[1];
        double n2 = pars[2];
        double gamma_min = pars[3];
        double gamma_break = pars[4];
        double gamma_max = pars[5];
        double K2 = K1 * pow(gamma_break, n2-n1);
        return k_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_break, gamma_max);
    }
}

//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
static PyObject* K(PyObject* self, PyObject* args){
    PyArrayObject* nu;
    double B;
    int SpectrumType;
    PyTupleObject* SpectrumPars;
    PyArg_ParseTuple(args, "O!diO!", &PyArray_Type, &nu, &B, &SpectrumType, &PyTuple_Type, &SpectrumPars);

    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ++ipar)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* in_array = nu;
    PyObject* out_array;
    out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = K_Sychrotron(*in_dataptr, B, SpectrumType, pars);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


static PyMethodDef Synchrotron[] = 
{
    {"F", (PyCFunction)F, METH_VARARGS, 'F(x)'},
    {"R", (PyCFunction)R, METH_VARARGS, 'R(x)'},
    {"J", (PyCFunction)J, METH_VARARGS, "calculate the emission coefficient of Synchrotron"},
    {"K", (PyCFunction)K, METH_VARARGS, "calculate the absorption coefficient of Synchrotron"},
    {NULL, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "Synchrotron",
        NULL,
        -1,
        Synchrotron,
        NULL,
        NULL,
        NULL,
        NULL
};

PyMODINIT_FUNC
PyInit_Synchrotron(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
