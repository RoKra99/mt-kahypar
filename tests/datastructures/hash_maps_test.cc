#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/hash_maps.h"
#include "mt-kahypar/datastructures/hash_functions.h"
#include "mt-kahypar/utils/timer.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(AHashMap, InitializesProperly) {
    HashMap<size_t, size_t, xxhash<size_t>> hm(100);
    ASSERT_TRUE(hm.empty());
}

TEST(AHashMap, InsertsElementsProperly) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.insert(i, i);
    }
    size_t count = 0;
    for (auto& e : hm.entries()) {
        ++count;
        ASSERT_TRUE(e.first == e.second);
    }
    ASSERT_TRUE(!hm.empty());
    ASSERT_TRUE(hm.size() == INSERTS);
    ASSERT_TRUE(count == INSERTS);
    for (size_t i = 0; i < INSERTS; ++i) {
        auto e = hm.get(i);
        ASSERT_TRUE(e.first == i && e.second == i);
    }
}

TEST(AHashMap, ErasesElementsProperly) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.insert(i, i);
    }
    ASSERT(hm.size() == INSERTS);
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.erase(i);
        for (auto& e : hm.entries()) {
            ASSERT_TRUE(e.first != i);
        }
        ASSERT_TRUE(hm.size() == INSERTS - i - 1);
    }
}

}
}