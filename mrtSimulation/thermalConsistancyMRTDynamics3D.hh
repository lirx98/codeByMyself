#include <algorithm>
#include <limits>

#include "core/latticeStatistics.h"
#include "thermalConsistancyMRTDynamics3D.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "thermalConsistancyMRTTemplates.h"  //修改为自定义的templates

namespace plb {
/* *************** Class TCGuoExternalForceMRTdynamics ***********************************************
 */
template <typename T, template <typename U> class Descriptor>
int TCGuoExternalForceMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, TCGuoExternalForceMRTdynamics<T, Descriptor> >(
        "MRT_ExternalForce_Guo_TC");

template <typename T, template <typename U> class Descriptor>
TCGuoExternalForceMRTdynamics<T, Descriptor>::TCGuoExternalForceMRTdynamics(T omega_) :
    ExternalForceDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
TCGuoExternalForceMRTdynamics<T, Descriptor>::TCGuoExternalForceMRTdynamics(
    HierarchicUnserializer &unserializer) :
    ExternalForceDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
TCGuoExternalForceMRTdynamics<T, Descriptor> *TCGuoExternalForceMRTdynamics<T, Descriptor>::clone()
    const
{
    return new TCGuoExternalForceMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int TCGuoExternalForceMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void TCGuoExternalForceMRTdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    //此处应重新定义一个struct，并重新定义内部的mrtCollisionWithForce函数。
    typedef TCmrtTemplates<T, Descriptor> mrtTemp;

    T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);

    Array<T, Descriptor<T>::d> u;
    this->computeVelocity(cell, u);

    T jSqr = mrtTemp::mrtCollisionWithForce(cell, rhoBar, u, this->getOmega(), (T)1);

    if (cell.takesStatistics()) {
        T invRho = Descriptor<T>::invRho(rhoBar);
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
T TCGuoExternalForceMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

}