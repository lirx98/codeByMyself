#include <cstdlib>
#include <iostream>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor
#define ADESCRITPOR  descriptors::AdvectionDiffusionWithSourceD3Q7Descriptor

class Params {
public:
    Params(
        T lx_, T ly_, T lz_, plint resolution_, T liquidRho_, T gasRho_, T wallRho_, T uMax_,
        Array<T, 3> center_, T radius_, Array<T, 3> gravityOrientation_) :
        lx(lx_),
        ly(ly_),
        lz(lz_),
        liquidRho(liquidRho_),
        gasRho(gasRho_),
        wallRho(wallRho_),
        resolution(resolution_),
        uMax(uMax_),
        center(center_),
        radius(radius_),
        gravityOrientation(gravityOrientation_) {};

    plint getNx() const
    {
        PLB_ASSERT(plint(lx * resolution) == lx * resolution);
        return lx * resolution;
    };
    plint getNy() const
    {
        PLB_ASSERT(plint(ly * resolution) == ly * resolution);
        return ly * resolution;
    };
    plint getNz() const
    {
        PLB_ASSERT(plint(lz * resolution) == lz * resolution);
        return lz * resolution;
    }
    T getRadius() const
    {
        return radius * (T)resolution;
    };
    Array<T, 3> getCenter() const
    {
        T cx = center[0] * resolution;
        T cy = center[1] * resolution;
        T cz = center[2] * resolution;
        Array<T, 3> centerInLattice(cx, cy, cz);
        return centerInLattice;
    }
    T getLiquidRho() const
    {
        return liquidRho;
    };
    T getGasRho() const
    {
        return gasRho;
    };
    T getWallRho() const
    {
        return wallRho;
    };
    T getDeltaX() const
    {
        return 1. / resolution;
    }
    T getLatticeU() const
    {
        return uMax;
    }
    T getDeltaT() const
    {
        return getLatticeU() / (T)resolution;
    }
    T getCv() const
    {
        return 1.;
    };
    Array<T, 3> getGravityOrientation() const
    {
        return gravityOrientation;
    };

private:
    T liquidRho, gasRho, wallRho, lx, ly, lz, uMax, radius;
    plint resolution;
    Array<T, 3> center, gravityOrientation;
};

class PsiCSEOS {
public:
    PsiCSEOS(T G_) : G(G_) {};
    T compute(T rho, T temp) const
    {
        T Temperature = temp * Tc;
        T mid1 = 4. * rho - 2. * rho * rho;
        T mid2 = 1. + mid1 / (1. - rho) / (1. - rho) / (1. - rho);
        T mid3 = R * Temperature * mid2 - a * rho - 1. / 3.;
        T psiSquare = abs(6. * rho * mid3 / G);
        return sqrt(psiSquare);
    };
    PsiCSEOS *clone() const
    {
        return new PsiCSEOS(*this);
    };

private:
    T G;
    const T R = 1., b = 4., a = 1.;
    const T Tc = 0.3773 * a / b / R;
};

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>


class NShanChenSingleComponentProcessor3D :
    public BoxProcessingFunctional3D_LL<T, nsDescriptor, T, aDescriptor> {
public:
    NShanChenSingleComponentProcessor3D(T G_, PsiCSEOS *Psi_, Params parameters_) :
        G(G_), Psi(Psi_), parameters(parameters_) {};
    virtual void process(
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
                    psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ) =
                        Psi->compute(rho, temp);
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
        Array<T, 3> gravityContibution(parameters.getGravityOrientation());
        //Array<T, 3> gravityContibution(0., 0., 0.);
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
                        forceContribution -= gravityContibution[iD] * 0.0002
                                             * (nsCell.computeDensity() - rhoAverage);
                        // Add interaction term.
                        T psi = psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                        forceContribution -= G * psi * rhoContribution[iD];
                        // Include into total momentum.
                        momentum[iD] += (T)1 / nsCell.getDynamics().getOmega() * forceContribution;
                        velocity.get(iX, iY, iZ)[iD] =
                            (momentum[iD]
                             - nsCell.getDynamics().getOmega() * forceContribution / 2.)
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
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    };
    virtual NShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor> *clone() const
    {
        return new NShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor>(*this);
    };

private:
    T G;
    PsiCSEOS *Psi;
    Params parameters;
};

T positionRho(plint iZ,Params parameters)
{
    return iZ < 0.75 * parameters.getNz() ? parameters.getLiquidRho() : parameters.getGasRho();
};

class initialDensityAndZeroVelocity {
public:
    initialDensityAndZeroVelocity(Params parameters_) : parameters(parameters_) {};
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T, 3> &u) const
    {
        rho = positionRho(iZ, parameters);
        u[0] = 0.;
        u[1] = 0.;
        u[2] = 0.;
    };

private:
    Params parameters;
};

class initialTemperatureOnBottomWall {
public:
    initialTemperatureOnBottomWall(Params parameters_) : parameters(parameters_) {};
    T operator()(plint iX, plint iY, plint iZ) const
    {
        T cx = parameters.getNx() / 2.;
        T cy = parameters.getNy() / 2.;
        T cz = 0.;
        T distance2 = (iX - cx) * (iX - cx) + (iY - cy) * (iY - cy);

        return distance2 > 100 ? 0.85 : 1.1;
    }

private:
    Params parameters;
};

void programSetup(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiBlockLattice3D<T, ADESCRITPOR> &adLattice,
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRITPOR>&adBoundaryCondition,
    Params &parameters)
{
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNx();
    T rhoWall = parameters.getWallRho();
    T rhoAverge = 0.75 * parameters.getLiquidRho() + 0.25 * parameters.getGasRho();
    Box3D heatedWall(0, nx - 1, 0, ny - 1, 0, 0);
    Box3D topWall(0, nx - 1, 0, ny - 1, nz - 1, nz - 1);

    //initialization of nsLattice and boundary
    initializeAtEquilibrium(
        nsLattice, nsLattice.getBoundingBox(), initialDensityAndZeroVelocity(parameters));
    setExternalVector(
        nsLattice, nsLattice.getBoundingBox(), NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
        Array<T, 3>(0., 0., 0.));

    defineDynamics(nsLattice, heatedWall, new BounceBack<T, NSDESCRIPTOR>(parameters.getWallRho()));
    defineDynamics(nsLattice, topWall, new BounceBack<T, NSDESCRIPTOR>(parameters.getWallRho()));

    //initialization of adLattice and boundary
    initializeAtEquilibrium(adLattice, adLattice.getBoundingBox(), 0.85, Array<T, 3>(0., 0., 0.));


    adBoundaryCondition.addTemperatureBoundary2N(heatedWall, adLattice);
    adBoundaryCondition.addTemperatureBoundary2P(topWall, adLattice);
    setBoundaryDensity(adLattice, heatedWall, initialTemperatureOnBottomWall(parameters));
    setBoundaryDensity(adLattice, topWall, 0.85);
    setExternalScalar(adLattice, adLattice.getBoundingBox(), ADESCRITPOR<T>::ExternalField::scalarBeginsAt, 0.);

    plint processLevel = 1;
    integrateProcessingFunctional(
        new NShanChenSingleComponentProcessor3D<T, NSDESCRIPTOR, ADESCRITPOR>(
            -1. / 3., new PsiCSEOS(-1. / 3.), parameters),
        nsLattice.getBoundingBox(), nsLattice, adLattice, processLevel);

    nsLattice.initialize();
    adLattice.initialize();
}

void writePLTs(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiBlockLattice3D<T, ADESCRITPOR> &adLattice,
    int iT, Params &parameters)
{
    enum {
        momentOffset = NSDESCRIPTOR<T>::ExternalField::momentumBeginsAt,
        sourceOffset = ADESCRITPOR<T>::ExternalField::scalarBeginsAt
    };
    stringstream valString;
    valString << setw(8) << setfill('0') << iT;
    string outName1 = "./tmp/rho" + valString.str() + ".plt";
    const char *outName;
    outName = outName1.c_str();
    plb_ofstream ofile(outName);
    ofile << "variables = x, y, z, rho, temperature" << endl;
    ofile << "zone i = " << parameters.getNx() << ", j = " << parameters.getNy()
          << ", k = " << parameters.getNz() << ", f = point" << endl;
    for (int iZ = 0; iZ < parameters.getNz(); iZ++) {
        for (int iY = 0; iY < parameters.getNy(); iY++) {
            for (int iX = 0; iX < parameters.getNx(); iX++) {
                ofile << fixed << setprecision(6) << iX << "\t" << iY << "\t" << iZ << "\t"
                      << nsLattice.get(iX, iY, iZ).computeDensity()
                      << "\t"
                      //<< *nsLattice.get(iX, iY, iZ).getExternal(momentOffset + 0) << "\t"
                      //<< *nsLattice.get(iX, iY, iZ).getExternal(momentOffset + 1) << "\t"
                      //<< *nsLattice.get(iX, iY, iZ).getExternal(momentOffset + 2) << "\t"
                      << adLattice.get(iX, iY, iZ).computeDensity() << "\t" << endl;
                      //<< *adLattice.get(iX, iY, iZ).getExternal(sourceOffset) << endl;
            }
        }
    }
}

void writeVTK(
    MultiBlockLattice3D<T, NSDESCRIPTOR>& nsLattice, MultiBlockLattice3D<T, ADESCRITPOR>& adLattice,
    plint iT, Params parameters)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iT, 6), 1.);
    vtkOut.writeData<3, float>(*computeVelocity(nsLattice), "velocity", 1.);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    T lx = 1., ly = 1., lz = 1., gasRho = 0.032, liquidRho = 0.279, wallRho = 0.12, radius = 0.3,
      uMax = 0.1;
    Array<T, 3> center(0.5, 0.5, 0.1);
    Array<T, 3> gravityOrientation(0., 0., 0.);
    plint resolution = 100;
    Params parameters(
        lx, ly, lz, resolution, liquidRho, gasRho, wallRho, uMax, center, radius,
        gravityOrientation);
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice(
        nx, ny, nz, new ExternalMomentBGKdynamics<T, NSDESCRIPTOR>(1.));
    MultiBlockLattice3D<T, ADESCRITPOR> adLattice(
        nx, ny, nz, new AdvectionDiffusionWithSourceBGKdynamics<T, ADESCRITPOR>(1.));
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, ADESCRITPOR> *adBoundaryCondition =
        createLocalAdvectionDiffusionBoundaryCondition3D<T, ADESCRITPOR>();

    nsLattice.periodicity().toggleAll(true);
    adLattice.periodicity().toggleAll(true);
    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);

    programSetup(nsLattice, adLattice, *adBoundaryCondition, parameters);

    for (int iT = 0; iT < 1001; iT++) {
        if (iT % 200 == 0) {
            writePLTs(nsLattice, adLattice, iT, parameters);
            writeVTK(nsLattice, adLattice, iT, parameters);
        }
        nsLattice.collideAndStream();
        adLattice.collideAndStream();
        pcout << iT << endl;
    }
}