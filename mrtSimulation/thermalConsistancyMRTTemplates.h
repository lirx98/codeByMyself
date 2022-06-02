#include "core/globalDefs.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/mrtLattices.h"

namespace plb {

template <typename T, class Descriptor>
struct TCmrtTemplatesImpl;

template <typename T, template <typename U> class Descriptor>
struct TCmrtTemplates {
    static T mrtCollisionWithForce(
        Cell<T, Descriptor> &cell, const T &rhoBar, const Array<T, Descriptor<T>::d> &u,
        const T &omega, T amplitude)
    {
        Array<T, Descriptor<T>::d> force;
        force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
        return TCmrtTemplatesImpl<T, typename Descriptor<T>::SecondBaseDescriptor>::
            mrtCollisionWithForce(cell.getRawPopulations(), rhoBar, u, omega, force, amplitude);
    }
};

template <typename T, class Descriptor>
struct TCmrtTemplatesImpl {
    static T equilibrium(plint iPop, T rhoBar, Array<T, Descriptor::d> const &j, const T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        T equ = T();
        for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
            equ += Descriptor::M[iPop][jPop]
                   * dynamicsTemplatesImpl<T, Descriptor>::bgk_ma2_equilibrium(
                       jPop, rhoBar, invRho, j, jSqr);
        }

        return equ;
    }

    static void computeInvM_S(T invM_S[Descriptor::q][Descriptor::q], const T &omega)
    {
        Array<T, Descriptor::q> s;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            s[iPop] = Descriptor::S[iPop];
        }
        for (plint iA = 0; iA < Descriptor::shearIndexes; ++iA) {
            plint iPop = Descriptor::shearViscIndexes[iA];
            s[iPop] = omega;
        }

        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                invM_S[iPop][jPop] = Descriptor::invM[iPop][jPop] * s[jPop];
            }
        }
    }

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibriumMoments(
        Array<T, Descriptor::q> &momentsEq, T rhoBar, Array<T, Descriptor::d> const &j,
        const T jSqr)
    {
        T invRho = Descriptor::invRho(rhoBar);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            momentsEq[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                momentsEq[iPop] += Descriptor::M[iPop][jPop]
                                   * dynamicsTemplatesImpl<T, Descriptor>::bgk_ma2_equilibrium(
                                       jPop, rhoBar, invRho, j, jSqr);
            }
        }
    }

    static void computeMoments(Array<T, Descriptor::q> &moments, const Array<T, Descriptor::q> &f)
    {
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            moments[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                moments[iPop] += Descriptor::M[iPop][jPop] * f[jPop];
            }
        }
    }

    static T mrtCollision(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &j,
        const T &omega)
    {
        Array<T, Descriptor::q> momentsEq;
        Array<T, Descriptor::q> moments;

        computeMoments(moments, f);
        T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);
        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);

        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);

        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }

        return jSqr;
    }

    static T mrtCollision(Array<T, Descriptor::q> &f, const T &omega)
    {
        Array<T, Descriptor::q> momentsEq;
        Array<T, Descriptor::q> moments;

        computeMoments(moments, f);
        T rhoBar = moments[0];
        Array<T, Descriptor::d> j;
        for (plint iA = 0; iA < Descriptor::jIndexes; ++iA) {
            plint iPop = Descriptor::momentumIndexes[iA];

            j[iA] = moments[iPop];
        }
        T jSqr = VectorTemplateImpl<T, Descriptor::d>::normSqr(j);
        computeEquilibriumMoments(momentsEq, rhoBar, j, jSqr);

        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);

        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }

        return jSqr;
    }

    static void addForceModify(Array<T, Descriptor::q> &forceMoments)
    {
        
        forceMoments[Descriptor::forceModifyIndex] += Descriptor::forceModify[forceModifyIndex];
    }

    static void addGuoForce(
        Array<T, Descriptor::q> &f, const Array<T, Descriptor::d> &force,
        Array<T, Descriptor::d> const &u, const T &omega, T amplitude)
    {
        Array<T, Descriptor::q> forcing;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T c_u = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                c_u += Descriptor::c[iPop][iD] * u[iD];
            }
            c_u *= Descriptor::invCs2 * Descriptor::invCs2;
            T forceTerm = T();
            for (int iD = 0; iD < Descriptor::d; ++iD) {
                forceTerm += (((T)Descriptor::c[iPop][iD] - u[iD]) * Descriptor::invCs2
                              + c_u * (T)Descriptor::c[iPop][iD])
                             * force[iD];
            }
            forceTerm *= Descriptor::t[iPop];
            forceTerm *= amplitude;
            forcing[iPop] = forceTerm;
        }

        Array<T, Descriptor::q> forceMoments;
        computeMoments(forceMoments, forcing);
        // add force modification
        addForceModify(forceMoments);

        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);

        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += -invM_S[iPop][jPop] * forceMoments[jPop];
            }
            collisionTerm *= (T)0.5;
            collisionTerm += forcing[iPop];
            f[iPop] += collisionTerm;
        }
    }


    
    static T mrtCollisionWithForce(
        Array<T, Descriptor::q> &f, const T &rhoBar, const Array<T, Descriptor::d> &u,
        const T &omega, const Array<T, Descriptor::d> &force, T amplitude)
    {
        Array<T, Descriptor::d> j = Descriptor::fullRho(rhoBar) * u;
        T jSqr = mrtCollision(f, rhoBar, j, omega);
        addGuoForce(f, force, u, omega, amplitude);
        addForceModify(f);
        return jSqr;
    }
};
}  // namespace plb