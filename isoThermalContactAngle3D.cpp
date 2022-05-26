#include <cstdlib>
#include <iostream>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor

class Params {
public:
    Params(
        T lx_, T ly_, T lz_, plint resolution_, T liquidRho_, T gasRho_, T wallRho_, T uMax_,
        Array<T, 3> center_, T radius_) :
        lx(lx_),
        ly(ly_),
        lz(lz_),
        liquidRho(liquidRho_),
        gasRho(gasRho_),
        wallRho(wallRho_),
        resolution(resolution_),
        uMax(uMax_),
        center(center_),
        radius(radius_) {};

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

private:
    T liquidRho, gasRho, wallRho, lx, ly, lz, uMax, radius;
    plint resolution;
    Array<T, 3> center;
};

class PsiCSEOS {
public:
    PsiCSEOS(T G_, T temp_) : G(G_), temp(temp_) {};
    T compute(T rho) const
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
    T G, temp;
    const T R = 1., b = 4., a = 1.;
    const T Tc = 0.3773 * a / b / R;
};

template <typename T, template <typename U> class Descriptor>
class NShanChenSingleComponentProcessor3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    NShanChenSingleComponentProcessor3D(T G_, PsiCSEOS *Psi_) : G(G_), Psi(Psi_) {};
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
    {  // Short-hand notation for the lattice descriptor
        typedef Descriptor<T> D;
        // Handle to external scalars
        enum {
            densityOffset = D::ExternalField::densityBeginsAt,
            momentumOffset = D::ExternalField::momentumBeginsAt
        };

        plint nx = domain.getNx() + 2;  // Include a one-cell boundary
        plint ny = domain.getNy() + 2;  // Include a one-cell boundary
        plint nz = domain.getNz() + 2;  // Include a one-cell boundary
        plint offsetX = domain.x0 - 1;
        plint offsetY = domain.y0 - 1;
        plint offsetZ = domain.z0 - 1;
        ScalarField3D<T> psiField(nx, ny, nz);

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
                    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                    T rho = cell.computeDensity();
                    // Evaluate potential function psi.
                    psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ) = Psi->compute(rho);
                    // Store density into the corresponding external scalar.
                    *cell.getExternal(densityOffset) = rho;
                    // Compute momentum through direct access to particle populations, and store
                    //   result in corresponding external scalars. Note that Cell::computeVelocity
                    //   cannot be used, because it returns the velocity of the external scalars,
                    //   not the velocity computed from the particle populations.
                    Array<T, Descriptor<T>::d> j;
                    momentTemplates<T, Descriptor>::get_j(cell, j);
                    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                        *(cell.getExternal(momentumOffset) + iD) = j[iD];
                    }
                }
            }
        }

        // Compute the interparticle forces, and store they by means of a
        //   velocity correction in the external velocity field.
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    Array<T, D::d> rhoContribution;
                    rhoContribution.resetToZero();
                    // Compute the term \sum_i ( t_i psi(x+c_i,t) c_i )
                    for (plint iPop = 0; iPop < D::q; ++iPop) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];
                        T psi = psiField.get(nextX - offsetX, nextY - offsetY, nextZ - offsetZ);
                        for (int iD = 0; iD < D::d; ++iD) {
                            rhoContribution[iD] += D::t[iPop] * psi * D::c[iPop][iD];
                        }
                    }

                    // Computation and storage of the final momentum, including tho momentum
                    //   difference due to interaction potential and the external force.
                    Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                    T *momentum = cell.getExternal(momentumOffset);
                    for (int iD = 0; iD < D::d; ++iD) {
                        // Initialize force contribution with force from external fields if there
                        //   is any, or with zero otherwise.
                        T forceContribution = getExternalForceComponent(cell, iD);
                        // Add interaction term.
                        T psi = psiField.get(iX - offsetX, iY - offsetY, iZ - offsetZ);
                        forceContribution -= G * psi * rhoContribution[iD];
                        // Include into total momentum.
                        momentum[iD] += (T)1 / cell.getDynamics().getOmega() * forceContribution;
                    }
                }
            }
        }
    };
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    };
    virtual NShanChenSingleComponentProcessor3D<T, Descriptor> *clone() const
    {
        return new NShanChenSingleComponentProcessor3D<T, Descriptor>(*this);
    };

private:
    T G;
    PsiCSEOS *Psi;
};

T contactIniRho(plint iX, plint iY, plint iZ, Params parameters)
{
    T cx = parameters.getCenter()[0];
    T cy = parameters.getCenter()[1];
    T cz = parameters.getCenter()[2];
    T distance2 =
        ((T)iX - cx) * ((T)iX - cx) + ((T)iY - cy) * ((T)iY - cy) + ((T)iZ - cz) * ((T)iZ - cz);
    if (distance2 > parameters.getRadius() * parameters.getRadius()) {
        return parameters.getGasRho();
    } else {
        return parameters.getLiquidRho();
    }
};

class initialDensityAndZeroVelocity {
public:
    initialDensityAndZeroVelocity(Params parameters_) : parameters(parameters_) {};
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T, 3> &u) const
    {
        rho = contactIniRho(iX, iY, iZ, parameters);
        u[0] = 0.;
        u[1] = 0.;
        u[2] = 0.;
    };

private:
    Params parameters;
};

void programSetup(MultiBlockLattice3D<T, DESCRIPTOR> &lattice, Params &parameters)
{
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNx();
    T rhoWall = parameters.getWallRho();

    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(), initialDensityAndZeroVelocity(parameters));

    defineDynamics(
        lattice, Box3D(0, nx - 1, 0, ny - 1, 0, 0),
        new BounceBack<T, DESCRIPTOR>(parameters.getWallRho()));
    defineDynamics(
        lattice, Box3D(0, nx - 1, 0, ny - 1, nz - 1, nz - 1),
        new BounceBack<T, DESCRIPTOR>(parameters.getWallRho()));

    plint processLevel = 1;

    integrateProcessingFunctional(
        new NShanChenSingleComponentProcessor3D<T, DESCRIPTOR>(
            -1. / 3., new PsiCSEOS(-1. / 3., 0.85)),
        lattice.getBoundingBox(), lattice, processLevel);

    lattice.initialize();
}

void writePLTs(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,Params &parameters)
{
    stringstream valString;
    valString << setw(8) << setfill('0') << parameters.getWallRho();
    string outName1 = "./tmp/rho" + valString.str() + ".plt";
    const char *outName;
    outName = outName1.c_str();
    plb_ofstream ofile(outName);
    ofile << "variables = x, y, z, rho" << endl;
    ofile << "zone i = " << parameters.getNx() << ", j = " << parameters.getNy()
          << ", k = " << parameters.getNz() << ", f = point" << endl;
    for (int iZ = 0; iZ < parameters.getNz(); iZ++) {
        for (int iY = 0; iY < parameters.getNy(); iY++) {
            for (int iX = 0; iX < parameters.getNx(); iX++) {
                ofile << fixed << setprecision(6) << iX << "\t" << iY << "\t" << iZ << "\t"
                      << lattice.get(iX, iY, iZ).computeDensity() << endl;
            }
        }
    }
}

void writeVKTs(MultiBlockLattice3D<T, DESCRIPTOR> &lattice, Params &parameters, int iT)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iT, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "rho", (T)1);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    for (int ind = 0; ind < 10; ind++) {
        T wallRho = 0.1143 + (T)ind * 0.00823;
        T lx = 1., ly = 1., lz = 1., gasRho = 0.032, liquidRho = 0.279, radius = 0.3, uMax = 0.1;
        Array<T, 3> center(0.5, 0.5, 0.3);
        plint resolution = 100;
        Params parameters(lx, ly, lz, resolution, liquidRho, gasRho, wallRho, uMax, center, radius);
        plint nx = parameters.getNx();
        plint ny = parameters.getNy();
        plint nz = parameters.getNz();
        MultiBlockLattice3D<T, DESCRIPTOR> lattice(
            nx, ny, nz, new ExternalMomentBGKdynamics<T, DESCRIPTOR>(1.));
        lattice.periodicity().toggleAll(true);
        lattice.toggleInternalStatistics(false);
        programSetup(lattice, parameters);

        pcout << "begin to simulate, rho_w = " << wallRho << endl;
        //global::timer("mainLoop").start();

        for (int iT = 0; iT < 501; iT++) {
            lattice.collideAndStream();
            // pcout << iT << endl;

            /*if (iT % 200 == 0) {
                writePLTs(lattice, iT, parameters);
            }*/
        }
        writePLTs(lattice, parameters);
        // T processTiming = global::timer("mainLoop").stop();
        // pcout << "Calculating domain is " << nx << " x " << ny << " x " << nz
        //      << ", and time for main loop is " << processTiming << endl; 
    }
    
}