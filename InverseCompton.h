#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

double Nu2E(double nu);
double E2Nu(double E);
double j_InverseCompton(double nu, double B, double R, int SpectrumType, double* pars, PyListObject* nuList, PyListObject* jList, PyListObject* kList);
static PyObject* J(PyObject* self, PyObject* args);