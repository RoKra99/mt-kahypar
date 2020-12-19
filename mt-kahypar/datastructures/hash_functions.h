#pragma once
#include "mt-kahypar/datastructures/xxhash.hpp"

namespace mt_kahypar {
namespace ds {

// credit to KaHIP repositiory
template <typename T>
struct xxhash {
    using hash_type = xxh::hash32_t;

    // explicit xxhash(uint32_t seed = 0)
    //     : _seed(seed)         {}

    xxhash() = default;

    // inline void reset(uint32_t seed) {
    //     _seed = seed;
    // }

    inline hash_type operator()(T x) const {
        return xxh::xxhash<32>({x});
    }

// private:
//     uint32_t _seed;
};

}
}