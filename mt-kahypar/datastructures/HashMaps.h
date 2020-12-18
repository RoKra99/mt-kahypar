#pragma once
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace ds {

template<typename Key, typename Value>
class TestHashMap {
    public:
    TestHashMap() = default;

    bool increment(Key id) {
        if (map.count(id) > 0) {
            ++map[id];
            return true;
        }
        return false;
    }

    bool decrement(Key id) {
        if (map.count(id) > 0) {
            if (map[id] == 1) {
                map.erase(id);
            } else {
                --map[id];
            }
            return true;
        }
        return false;
    }

    void insert(Key id, Value v) {
        map.insert(id, v);
    }


    private:
     std::unordered_map<Key, Value> map;
};
}
}