// NTUtilities: constants and functions of use
#ifndef NTUTILITIES_HH
#define NTUTILITIES_HH


#include <cmath>
#include <vectors>


typedef std::vector<double> doubles;

// thermal energy (25 C)
const double Etherm = 0.025692579120652998*1e-6;

// bounds of fast energies
const double Emin = 1e-7, Emax = 20.; 

// kernel of a Maxwell-Boltzmann thermal distribution of energies
double MaxwellBoltzmannKernel(double E) {
    return std::sqrt(E)*std::exp(-E/Etherm);
}

// assumed kernel of fast flux
double fast_kernel(double E) {
    return 1./E;
}

// trapezoidal integration of (x,y) data
double trap(const doubles &x, const doubles &y) {
    size_t n = x.size();
    if (n != y.size()) {
        throw std::invalid_argument("Invalid arguments to function trap(). " 
        "Vectors must have the same length.");
    }
    double s = 0;
    for (size_t i = 1; i < n; ++i) {
        s += (y[i] + y[i-1])*(x[i] - x[i-1])/2.;
    }
    return s;
}

double interp(double x, double x1, double y1, double x2, double y2) {
    return y1 + (x - x1)*(y2 - y1)/(x2 - x1);
}


#endif // NTUTILITIES_HH