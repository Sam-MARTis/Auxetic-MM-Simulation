#pragma once

#include <cmath>
#include <stdint.h>

#include "constants.hpp"
#include "vec.hpp"


inline double randValueMine(const double min, const double max)
{
    return min + (max-min)*(1.0*rand()/(1.0*RAND_MAX));
}


// The legacy hand-rolled 2D vector type now delegates to the generic Vec struct.
using Vec2D = Vec<double, 2>;
