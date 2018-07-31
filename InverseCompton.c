/*
 * Follow the paper A&A 367, 809-825 (2001)
 * THe multifrequency emission of Mrk 501 from radio to TeV gamma-rays
 * 
 * gamma * m_e * c**2 is the energy of electron;
 * Es * m_e * c**2 is the energy of Synchrotron photon;
 * Ec * m_e * c**2 is the energy of InverseCompton photon.
 */ 

//#include "Synchrotron.h"
#include "PhysicalConstant.h"
#include "ElectronDistribution.h"
#include "InverseCompton.h"
#include <gsl/gsl_spline.h>

double Kappa(double Ec, double gamma, double Es){
    return Ec / (4 * Es * gamma * (gamma - Ec));
}

double CrossSection(double Ec, double gamma, double Es){
    double k = Kappa(Ec, gamma, Es);
    return 2 * M_PI * r_e * r_e * c / (gamma * gamma * Es)\
    * (2 * k * log(k) + (1+2*k)*(1-k) + (4*Es*gamma*k)*(4*Es*gamma*k) / (2 * (1+4*Es*gamma*k)) * (1-k));
}

double N_Es(double Es, double j, double k, double R){
    return (3./4) * (4*M_PI) / (h*c*Es) * (j/k) * (1 - exp(-k*R));
}

double Emin(double gamma, double Es){
    return Es;
}

double Emax(double gamma, double Es){
    return gamma * (4 * Es * gamma) / (1 + 4 * Es * gamma);
}

double Interp_log10(double xi, double* x, double* y, int N){
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline, x, y, N);
    double yi = gsl_spline_eval(spline, xi, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return yi;
}

void GetArray_log10(PyObject* List, double* Array){
    Py_ssize_t NList = PyList_Size(List);
    for(int i=0; i<NList; ++i){
        Array[i] = log10(PyFloat_AsDouble(PyList_GetItem(List, i)));
        if(Array[i] < -99)
            Array[i] = -99;
    }
}

double Q(double Ec, double B, double R, int SpectrumType, double* pars, 
        PyListObject* nuList, PyListObject* jList, PyListObject* kList){
    Py_ssize_t NList = PyList_Size(nuList);
    double nu_min = PyFloat_AsDouble(PyList_GetItem(nuList, 0));
    double nu_max = PyFloat_AsDouble(PyList_GetItem(nuList, NList-1));

    int Nbins_Es = 300;
    double Es_min = Nu2E(nu_min);
    double Es_max = Nu2E(nu_max);
    double dlog_Es = (log10(Es_max) - log10(Es_min)) / Nbins_Es;

    int Nbins_gamma = 200;
    double gamma_min;
    double gamma_max;
    if(SpectrumType == 1){
        gamma_min = pars[2];
        gamma_max = pars[3];
    }
    if(SpectrumType == 2){
        gamma_min = pars[3];
        gamma_max = pars[5];
    }
    double dlog_gamma = (log10(gamma_max) - log10(gamma_min)) / Nbins_gamma;

    double q = 0;
    for(int i=0; i<Nbins_Es; ++i){
        double Es = pow(10, log10(Es_min) + dlog_Es * i);
        double nu = E2Nu(Es);

        //double j = J_Sychrotron(nu, B, SpectrumType, pars);
        //double k = K_Sychrotron(nu, B, SpectrumType, pars);
        double nuArray[NList];
        double jArray[NList];
        double kArray[NList];
        GetArray_log10(nuList, nuArray);
        GetArray_log10(jList, jArray);
        GetArray_log10(kList, kArray);
        double j = pow(10, Interp_log10(log10(nu), nuArray, jArray, NList));
        double k = pow(10, Interp_log10(log10(nu), nuArray, kArray, NList));
        double n_Es = N_Es(Es, j, k, R);

        double int_gamma = 0;
        for(int j=0; j<Nbins_gamma; ++j){
            double gamma = pow(10, log10(gamma_min) + dlog_gamma * i);
            if((Ec < Es) || (Ec > Emax(gamma, Es))) continue;
            double n_gamma = N_gamma(gamma, SpectrumType, pars);
            double C = CrossSection(Ec, gamma, Es);
            int_gamma += C * n_gamma * gamma * dlog_gamma;
        }
        q += int_gamma * n_Es * Es * dlog_Es;
    }
    return q;
}

double Nu2E(double nu){
    double E = (h * nu) / (m_e * c * c);
    return E;
}

double E2Nu(double E){
    double nu = (E * m_e * c * c) / h;
    return nu;
}

double j_InverseCompton(double nu, double B, double R, int SpectrumType, double* pars, PyListObject* nuList, PyListObject* jList, PyListObject* kList){
    double Ec = Nu2E(nu);
    double q = Q(Ec, B, R, SpectrumType, pars, nuList, jList, kList);
    return h / (4 * M_PI) * Ec * q;
}

static PyObject* J(PyObject* self, PyObject* args){
    PyArrayObject* nu;
    double B;
    double R;
    int SpectrumType;
    PyTupleObject* SpectrumPars;
    PyListObject* nuList;
    PyListObject* jList;
    PyListObject* kList;
    PyArg_ParseTuple(args, "O!ddiO!O!O!O!", &PyArray_Type, &nu, &B, &R, &SpectrumType, &PyTuple_Type, &SpectrumPars, &PyList_Type, &nuList, &PyList_Type, &jList, &PyList_Type, &kList);

    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ++ipar)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    double nu_min = PyFloat_AsDouble(PyArray_Min(nu, 0, NULL));
    double nu_max = PyFloat_AsDouble(PyArray_Max(nu, 0, NULL));
    PyArrayObject* in_array = nu;
    PyObject* out_array;
    out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = j_InverseCompton(*in_dataptr, B, R, SpectrumType, pars, nuList, jList, kList);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


static PyMethodDef InverseCompton[] = 
{
    {"J", (PyCFunction)J, METH_VARARGS, "calculate the emission coefficient of InverseCompton"},
    {NULL, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "InverseCompton",
        NULL,
        -1,
        InverseCompton,
        NULL,
        NULL,
        NULL,
        NULL
};

PyMODINIT_FUNC
PyInit_InverseCompton(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
