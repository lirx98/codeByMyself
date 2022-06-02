#ifndef TC_MRT_SHANCHEN_DATA_PROCESSOR_H
#define TC_MRT_SHANCHEN_DATA_PROCESSOR_H

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



namespace plb {

template <
    typename T, template <typename U1> class nsDescriptor, template <typename U2> class aDescriptor>
class TCMRTShanChenSingleComponentProcessor3D :
    public BoxProcessingFunctional3D_LL<T, nsDescriptor, T, aDescriptor> {
public:
    TCMRTShanChenSingleComponentProcessor3D(
        T G_, interparticlePotential::cePseudoPotential<T> *Psi_,
        ShanChenParamsBase<T> parameters_);
    ~TCMRTShanChenSingleComponentProcessor3D();

    virtual void process(
        Box3D domain, BlockLattice3D<T, nsDescriptor> &nsLattice,
        BlockLattice3D<T, aDescriptor> &adLattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual TCMRTShanChenSingleComponentProcessor3D<T, nsDescriptor, aDescriptor> *clone() const;

private:
    T G;
    interparticlePotential::cePseudoPotential<T> *Psi;
    ShanChenParamsBase<T> parameters;
};
};


#endif  // !TC_MRT_SHANCHEN_DATA_PROCESSOR_H
