/** \file
 * This object is a MRT LB dynamics with modification 
 *for thermaldynamic consistancy in 3D case. 
 */
#include "basicDynamics/externalForceDynamics.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "core/globalDefs.h"
#include "complexDynamics/mrtDynamics.h"

namespace plb {
/// Implementation of the MRT collision step with thermal consitancy
template <typename T, template <typename U> class Descriptor>
class TCGuoExternalForceMRTdynamics : public ExternalForceDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    TCGuoExternalForceMRTdynamics(T omega_);
    TCGuoExternalForceMRTdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual TCGuoExternalForceMRTdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    static int id;
};
};  // namespace plb