#ifndef MRT_LATTICES_H
#define MRT_LATTICES_H

#include <vector>

#include "core/globalDefs.h"
#include "latticeBoltzmann/externalFields.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "multiPhysics/shanChenLattices2D.h"
#include "multiPhysics/shanChenLattices3D.h"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"

namespace plb {
namespace descriptors {

template <typename T>
struct TCMRTD3Q19DescriptorBase : public D3Q19DescriptorBase<T> {
    typedef D3Q19DescriptorBase<T> BaseDescriptor;
    typedef TCMRTD3Q19DescriptorBase<T> SecondBaseDescriptor;
    enum { numPop = BaseDescriptor::q };

    static const T M[BaseDescriptor::q][BaseDescriptor::q];  // Matrix of base change between f and
                                                             // moments : moments=M.f
    static const T invM[BaseDescriptor::q]
                       [BaseDescriptor::q];  // inverse of base change matrix : f=invM.moments
    static const T S[BaseDescriptor::q];     // relaxation times
    enum { jIndexes = 3 };
    static const int momentumIndexes[jIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { shearIndexes = 5 };
    static const int
        shearViscIndexes[shearIndexes];  // relevant indexes of r. t. for shear viscosity
    static const int bulkViscIndex = 1;  // relevant index of r. t. for bulk viscosity
    enum { qIndexes = 3 };
    static const int
        qViscIndexes[qIndexes];  // relevant indexes of r. t. for q indices (s4,s6,s8 of the paper)
    static const int epsilonIndex =
        2;  // relevant index of r. t. for epsilon (s2 in the original paper)
    static const T forceModify[BaseDescriptor::q];
    static const int forceModifyIndex = 1;
};

template <typename T>
struct TCForcedMRTShanChenD3Q19Descriptor :
    public MRTD3Q19DescriptorBase<T>,
    public ForcedShanChenExternalBase3D {
    static const char name[];
};
};  // namespace descriptors
};
#endif
