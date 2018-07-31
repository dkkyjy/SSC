#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

double j_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double J_PowerLaw(double nu, double B, double* pars);
double j_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double J_BrokenPowerLaw(double nu, double B, double* pars);
double J_Sychrotron(double nu, double B, int SpectrumType, double* pars);
static PyObject* J(PyObject* self, PyObject* args);

double k_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double K_PowerLaw(double nu, double B, double* pars);
double k_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double K_BrokenPowerLaw(double nu, double B, double* pars);
double K_Sychrotron(double nu, double B, int SpectrumType, double* pars);
static PyObject* K(PyObject* self, PyObject* args);
