#ifndef EOS_PSEUDO_POTENTIAL_H
#define EOS_PSEUDO_POTENTIAL_H

namespace plb {
namespace interparticlePotential {

template <typename T>
class eosPseudoPotentialBase {
public:
    // compute pseudo-potential
    T compute(T rho, T temp) const = 0;

    eosPseudoPotentialBase<T> *clone() const = 0;
};

template <typename T>
class cePseudoPotential : public eosPseudoPotentialBase<T> {
public:
    cePseudoPotential(T G_);
    T compute(T rho, T temp) const;
    cePseudoPotential<T> *clone() const;

private:
    const T G;
    const T R = 1., b = 4., a = 1.;
    const T Tc = 0.3773 * a / b / R;
};
};  // namespace interparticlePotential
};  // namespace plb

#endif  // !EOS_PSEUDO_POTENTIAL
