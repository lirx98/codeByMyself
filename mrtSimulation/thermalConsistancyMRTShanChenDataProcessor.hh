#ifndef TC_MRT_SHANCHEN_DATA_PROCESSOR_HH
#define TC_MRT_SHANCHEN_DATA_PROCESSOR_HH

#include "thermalConsistancyMRTLattice3D.h"
#include "thermalConsistancyMRTLattice3D.hh"
#include "EOSpseudoPotential.h"
#include "EOSpseudoPotential.hh"

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiPhysics/interparticlePotential.h"
#include "ShanChenParams.h"
#include "ShanChenParams.hh"
#include "thermalConsistancyMRTShanChenDataProcessor.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>::
    TCMRTShanChenSingleComponentProcessor3D(
        T G_, interparticlePotential::cePseudoPotential<T> *Psi_,
        ShanChenParamsBase<T> parameters_) :
    G(G_), Psi(Psi_), parameters(parameters_) {};

//此函数还需进一步更改，将forceModify项加入mrt格子Descirptor中,并将力项对动量的影响消除
template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
void TCMRTShanChenSingleComponentProcessor3D<T,nsDescriptor,aDescriptor>:: process(
    Box3D domain, BlockLattice3D<T, nsDescriptor> &nsLattice,
    BlockLattice3D<T, aDescriptor> &adLattice)
{  // Short-hand notation for the lattice descriptor
    typedef nsDescriptor<T> nsD;
    typedef aDescriptor<T> aD;
    // Handle to external scalars
    enum {
        densityOffset = nsD::ExternalField::densityBeginsAt,
        momentumOffset = nsD::ExternalField::momentumBeginsAt,
        sourceOffset = aD::ExternalField::scalarBeginsAt
    };

    T rhoAverage = 0.75 * parameters.getLiquidRho() + 0.25 * parameters.getGasRho();
    plint nx = domain.getNx() + 2;  // Include a one-cell boundary
    plint ny = domain.getNy() + 2;  // Include a one-cell boundary
    plint nz = domain.getNz() + 2;  // Include a one-cell boundary
    plint offsetX = domain.x0 - 1;
    plint offsetY = domain.y0 - 1;
    plint offsetZ = domain.z0 - 1;
    ScalarField3D<T> psiField(nx, ny, nz);
    ScalarField3D<T> divergence(nx, ny, nz);
    TensorField3D<T, 3> velocity(nx, ny, nz);
    // Compute density and momentum on every site and store result in external scalars;
    //   furthermore, evaluate the interaction potential Psi and store it into a ScalarField.
    //   Envelope cells are included, because they are needed to compute the interaction
    //   potential in the following. Note that the value of the momentum is stored temporarily
    //   only, as it is corrected later on to include corrections due to the interaction
    //   potential.
    for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
        for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
            for (plint iZ = domain.z0 - 1; iZ <= domain.z1 + 1; ++iZ) {
                // Get "intelligent" value of density through cell object, to account
                //   for the fact that the density value can be user-defined, for example
                //   on boundaries.
                Cell<T, nsDescriptor> &nsCell = nsLattice.get(iX, iY, iZ);
                Cell<T, aDescriptor> &adCell = adLattice.get(iX, iY, iZ);
                T rho = nsCell.computeDensity();
                T temp = adCell.computeDensity();
                // Evaluate potential function psi.
                psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ) = Psi->compute(rho, temp);
                // Store potential into the modification term
                nsD::forceModify[nsD::forceModifyIndex] =
                    psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ);

                // Store density into the corresponding external scalar.
                *nsCell.getExternal(densityOffset) = rho;
                // Compute momentum through direct access to particle populations, and store
                //   result in corresponding external scalars. Note that Cell::computeVelocity
                //   cannot be used, because it returns the velocity of the external scalars,
                //   not the velocity computed from the particle populations.
                Array<T, nsDescriptor<T>::d> j;
                momentTemplates<T, nsDescriptor>::get_j(nsCell, j);
                for (int iD = 0; iD < nsDescriptor<T>::d; ++iD) {
                    *(nsCell.getExternal(momentumOffset) + iD) = j[iD];
                }
            }
        }
    }

    // Compute the interparticle forces, and store they by means of a
    //   velocity correction in the external velocity field.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, nsD::d> rhoContribution;
                rhoContribution.resetToZero();
                // Compute the term \sum_i ( t_i psi(x+c_i,t) c_i )
                for (plint iPop = 0; iPop < nsD::q; ++iPop) {
                    plint nextX = iX + nsD::c[iPop][0];
                    plint nextY = iY + nsD::c[iPop][1];
                    plint nextZ = iZ + nsD::c[iPop][2];
                    T psi = psiField.get(nextX - offsetX, nextY - offsetY, nextZ - offsetZ);
                    for (int iD = 0; iD < nsD::d; ++iD) {
                        rhoContribution[iD] += nsD::t[iPop] * psi * nsD::c[iPop][iD];
                    }
                }

                // Computation and storage of the final momentum, including tho momentum
                //   difference due to interaction potential and the external force.
                Cell<T, nsDescriptor> &nsCell =
                    nsLattice.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                T *momentum = nsCell.getExternal(momentumOffset);
                for (int iD = 0; iD < nsD::d; ++iD) {
                    // Initialize force contribution with force from external fields if there
                    //   is any, or with zero otherwise.
                    T forceContribution = getExternalForceComponent(nsCell, iD);
                    forceContribution -=
                        gravityContibution[iD] * 0.0002 * (nsCell.computeDensity() - rhoAverage);
                    // Add interaction term.
                    T psi = psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                    forceContribution -= G * psi * rhoContribution[iD];
                    // Include into total momentum.
                    momentum[iD] += (T)1 / nsCell.getDynamics().getOmega() * forceContribution;
                    velocity.get(iX, iY, iZ)[iD] =
                        (momentum[iD] - nsCell.getDynamics().getOmega() * forceContribution / 2.)
                        / nsCell.computeDensity();
                }
            }
        }
    }

    // compute divergence of velocity and calculate temperature caused by latent heat
    Dot3D offset = computeRelativeDisplacement(divergence, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                plint iX2 = iX + offset.x;
                plint iY2 = iY + offset.y;
                plint iZ2 = iZ + offset.y;
                divergence.get(iX - offsetX, iY - offsetY, iZ - offsetZ) =
                    fdDataField::bulkXderiv(velocity, iX2, iY2, iZ2, 0)
                    + fdDataField::bulkXderiv(velocity, iX2, iY2, iZ2, 1)
                    + fdDataField::bulkXderiv(velocity, iX2, iY2, iZ2, 2);
            }
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                Cell<T, nsDescriptor> nsCell =
                    nsLattice.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                Cell<T, aDescriptor> adCell =
                    adLattice.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                T *rho = nsCell.getExternal(densityOffset);
                T temp = adCell.computeDensity();
                T Cv = parameters.getCv();
                T midst1 = 1. + *rho + *rho * *rho - *rho * *rho * *rho;
                T midst2 = midst1 / (1. - *rho) / (1. - *rho) / (1. - *rho);
                *adCell.getExternal(sourceOffset) =
                    temp * (1. - midst2 / *rho / Cv) * divergence.get(iX, iY, iZ);
            }
        }
    }
};

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>* TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>::clone() const
{
    return new TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>(*this);
};

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
void TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
};

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
TCMRTShanChenSingleComponentProcessor3D<
    T, nsDescriptor, aDescriptor>::~TCMRTShanChenSingleComponentProcessor3D()
{
    delete Psi;
}
};  // namespace plb

#endif  // !TC_MRT_SHANCHEN_DATA_PROCESSOR_HH
