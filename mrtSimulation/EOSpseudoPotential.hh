#ifndef EOS_PSEUDO_POTENTIAL_HH
#define EOS_PSEUDO_POTENTIAL_HH

#include "EOSpseudoPotential.h"
#include <cmath>

namespace plb {
namespace interparticlePotential {

template <typename T>
cePseudoPotential<T>::cePseudoPotential(T G_) : G(G_) {};

template <typename T>
T cePseudoPotential<T>::compute(T rho, T temp) const
{
    T Temperature = temp * Tc;
    T mid1 = 4. * rho - 2. * rho * rho;
    T mid2 = 1. + mid1 / (1. - rho) / (1. - rho) / (1. - rho);
    T mid3 = R * Temperature * mid2 - a * rho - 1. / 3.;
    T psiSquare = abs(6. * rho * mid3 / G);
    return sqrt(psiSquare);
};

template <typename T>
cePseudoPotential<T> *cePseudoPotential<T>::clone() const
{
    return new cePseudoPotential<T>(*this);
}

};  // namespace interparticlePotential

};  // namespace plb

#endif