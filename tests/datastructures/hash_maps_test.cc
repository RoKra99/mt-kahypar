#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/hash_maps.h"
#include "mt-kahypar/datastructures/hash_functions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/randomize.h"

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

TEST(AHashMap, HoldsTheSameElementsAfterResizing) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(8);
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

TEST(AHashMap, ContainsOnlyTheInsertedElements) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    ASSERT_TRUE(hm.empty());
    for (size_t i = 0; i < INSERTS; ++i) {
        ASSERT_FALSE(hm.contains(i));
    }
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.insert(i, i);
    }
    for (size_t i = 0; i < INSERTS; ++i) {
        ASSERT_TRUE(hm.contains(i));
    }
    for (size_t i = INSERTS; i < 2 * INSERTS; ++i) {
        ASSERT_FALSE(hm.contains(i));
    }
}

TEST(AHashMap, DecrementsValuesProperlyWithoutRemovingThem) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    for (size_t i = INSERTS; i < INSERTS*2; ++i) {
        hm.insert(i,i);
        hm.decrement(i);
    }
    ASSERT_TRUE(hm.size() == INSERTS);
    for (size_t i = INSERTS; i < INSERTS*2; ++i) {
        ASSERT_TRUE(hm.get(i).second == i-1);
    }
}


TEST(AHashMap, RemovesElementsDecrementedToZero) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.insert(i, 1);
    }
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.decrement(i);
        ASSERT_TRUE(hm.size() == INSERTS - i - 1);
        ASSERT_TRUE(hm.get(i).first != i);
        for (auto e : hm.entries()) {
            ASSERT_TRUE(e.first > i);
        }
    }
}

TEST(AHashMap, IncrementsValuesProperlyWithoutRemovingThem) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    for (size_t i = INSERTS; i < INSERTS*2; ++i) {
        hm.insert(i,i);
        hm.increment(i);
    }
    ASSERT_TRUE(hm.size() == INSERTS);
    for (size_t i = INSERTS; i < INSERTS*2; ++i) {
        ASSERT_TRUE(hm.get(i).second == i+1);
    }
}

TEST(AHashMap, DoesNotIncrementNotExistingEntries) {
    const size_t INSERTS = 10;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(100);
    ASSERT_FALSE(hm.increment(1));
    ASSERT_TRUE(hm.empty());
    for (size_t i = INSERTS; i < INSERTS*2; ++i) {
        hm.insert(i,i);
    }
    for (size_t i = 0; i < INSERTS; ++i) {
        ASSERT_FALSE(hm.increment(i));
    }
    for (auto& e : hm.entries()) {
        ASSERT_EQ(e.first, e.second);
    }
}


TEST(AHashMap, ActsLikeAnUnorderedMapWithInsertAndErase) {
    const size_t INSERTS = 50000;
    HashMap<size_t, size_t, xxhash<uint32_t>> hm(1024);
    size_t hm_result = 0;
    for (size_t i = 0; i < INSERTS; ++i) {
        hm.insert(i, i);
        for (auto& e : hm.entries()) {
            hm_result += e.second;
        }
    }

    for (size_t i = 0; i < INSERTS; ++i) {
        hm.erase(i);
        for (auto& e : hm.entries()) {
            hm_result += e.second;
        }
    }

    std::unordered_map<size_t, size_t> map;
    size_t map_result = 0;
    for (size_t i = 0; i < INSERTS; ++i) {
        map.insert(std::make_pair(i,i));
        for (auto& e : map) {
            map_result += e.second;
        }
    }

    for (size_t i = 0; i < INSERTS; ++i) {
        map.erase(i);
        for (auto& e : map) {
            map_result += e.second;
        }
    }

    ASSERT_EQ(hm_result, map_result);
}

}
}