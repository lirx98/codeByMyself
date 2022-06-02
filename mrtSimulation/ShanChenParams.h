#ifndef SHANCHEN_PARAMS_H
#define SHANCHEN_PARAMS_H

#include "core/array.h"

namespace plb {

template <typename T>
class ShanChenParamsBase {
public:
    ShanChenParamsBase(
        T lx_, T ly_, T lz_, plint resolution_, T liquidRho_, T gasRho_, T wallRho_, T uMax_,
        T Cv_, Array<T, 3> gravityOrientation_);

    plint getNx() const;
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

    T getLiquidRho() const;
    T getGasRho() const;
    T getWallRho() const;

    T getDeltaX() const;

    T getLatticeU() const;

    T getDeltaT() const;
    T getCv() const
    {
        return 1.;
    };
    Array<T, 3> getGravityOrientation() const;

private:
    T liquidRho, gasRho, wallRho, lx, ly, lz, uMax, Cv;
    plint resolution;
    Array<T, 3> gravityOrientation;
};

};      // namespace plb
#endif  // !SHANCHEN_PARAMS_H
