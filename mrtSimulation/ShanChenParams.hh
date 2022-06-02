#ifndef SHANCHEN_PARAMS_HH
#define SHANCHEN_PARAMS_HH

#include "ShanChenParams.h"
#include "core/array.h"

namespace plb {

template <typename T>
inline ShanChenParamsBase<T>::ShanChenParamsBase(
    T lx_, T ly_, T lz_, plint resolution_, T liquidRho_, T gasRho_, T wallRho_, T uMax_,T Cv_,
    Array<T, 3> gravityOrientation_) :
    lx(lx_),
    ly(ly_),
    lz(lz_),
    resolution(resolution_),
    liquidRho(liquidRho_),
    gasRho(gasRho_),
    wallRho(wallRho_),
    uMax(uMax_),
    Cv(Cv_),
    gravityOrientation(gravityOrientation_) {};

template <typename T>
plint ShanChenParamsBase<T>::getNx() const
{
    PLB_ASSERT(plint(lx * resolution) == lx * resolution);
    return lx * resolution;
};

template <typename T>
plint ShanChenParamsBase<T>::getNy() const
{
    PLB_ASSERT(plint(ly * resolution) == ly * resolution);
    return ly * resolution;
};

template <typename T>
plint ShanChenParamsBase<T>::getNz() const
{
    PLB_ASSERT(plint(lz * resolution) == l * resolution);
    return lz * resolution;
};

template <typename T>
T ShanChenParamsBase<T>::getLiquidRho() const
{
    return liquidRho;
};

template <typename T>
T ShanChenParamsBase<T>::getGasRho() const
{
    return gasRho;
};

template <typename T>
T ShanChenParamsBase<T>::getWallRho() const
{
    return wallRho;
};

template <typename T>
T ShanChenParamsBase<T>::getDeltaX() const
{
    return 1. / resolution;
};

template <typename T>
T ShanChenParamsBase<T>::getLatticeU() const
{
    return uMax;
}

template <typename T>
T ShanChenParamsBase<T>::getDeltaT() const
{
    return getLatticeU() / (T)resolution;
}
template <typename T>
T ShanChenParamsBase<T>::getCv() const
{
    return Cv;
}

template <typename T>
Array<T, 3> ShanChenParamsBase<T>::getGravityOrientation() const
{
    return gravityOrientation;
}

};
#endif  // !SHANCHEN_PARAMS_HH
