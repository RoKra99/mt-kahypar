#pragma once

namespace mt_kahypar::math {

inline HyperedgeWeight fast_power2(HyperedgeWeight base, HyperedgeWeight exponent) {
    auto fast_power_acc = [&](HyperedgeWeight base2, HyperedgeWeight exponent2, HyperedgeWeight acc_i, HyperedgeWeight acc_j, auto& ref) mutable {
        if (exponent2 == 0) return 1;
        else if (exponent2 == 1) {
            ASSERT(acc_i <= std::numeric_limits<HyperedgeWeight>::max() / acc_j);
            return acc_i * acc_j;
        }
        else if (exponent2 % 2 == 0) {
            ASSERT(acc_i <= std::numeric_limits<HyperedgeWeight>::max() / acc_i);
            return ref(base2, exponent2 / 2, acc_i * acc_i, acc_j, ref);
        }
        else {
            ASSERT(acc_i <= std::numeric_limits<HyperedgeWeight>::max() / acc_j);
            return ref(base2, exponent2 - 1, acc_i, acc_i * acc_j, ref);
        }
    };
    return fast_power_acc(base, exponent, base, 1, fast_power_acc);
}

inline HyperedgeWeight fast_power(HyperedgeWeight base, HyperedgeWeight exp) {
    HyperedgeWeight result = 1;
    const HyperedgeWeight int_max = std::numeric_limits<HyperedgeWeight>::max();
    while (true) {
        if (exp & 1) {
            assert(result < int_max / base);
            result *= base;
        }
        exp >>= 1;
        if (!exp)
            break;
        assert(base < int_max / base);
        base *= base;
    }
    return result;
}
}