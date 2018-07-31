#include <math.h>
#include "ElectronDistribution.h"
double N_PowerLaw(double gamma, double K, double n, double gamma_min, double gamma_max){
    if((gamma > gamma_min) && (gamma < gamma_max))
        return K * pow(gamma, -n);
    else
        return 0;
}

double N_BrokePowerLaw(double gamma, double K, double n1, double n2, double gamma_min, double gamma_break, double gamma_max){
    if((gamma > gamma_min) && (gamma < gamma_break))
        return K * pow(gamma, -n1);
    else{ 
        if ((gamma > gamma_break) && (gamma < gamma_max))
            return K * pow(gamma_break, n2-n1) * pow(gamma, -n2);
        else
            return 0;
    }
}

double N_gamma(double gamma, int SpectrumType, double* pars){
    if(SpectrumType == 1){
        double K = pars[0];
        double n = pars[1];
        double gamma_min = pars[2];
        double gamma_max = pars[3];
        double n_gamma = N_PowerLaw(gamma, K, n, gamma_min, gamma_max);
        return n_gamma;
    }
    if(SpectrumType == 2){
        double K = pars[0];
        double n1 = pars[1];
        double n2 = pars[2];
        double gamma_min = pars[3];
        double gamma_break = pars[4];
        double gamma_max = pars[5];
        double n_gamma = N_BrokePowerLaw(gamma, K, n1, n2, gamma_min, gamma_break, gamma_max);
        return n_gamma;
    }
}
