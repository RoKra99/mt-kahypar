#pragma once
#include "mt-kahypar/utils/floating_point_comparisons.h"

namespace mt_kahypar::math {

inline Volume fast_power(Volume base, HyperedgeWeight exp) {
    
    Volume result = 1;
    //const Volume int_max = std::numeric_limits<Volume>::max();
    while (true) {
        if (exp & 1) {
            //assert(result < int_max / base);
            result *= base;
        }
        exp >>= 1;
        if (!exp)
            break;
        //assert(base < int_max / base);
        base *= base;
    }
    return result;
}
}