#pragma once
#include "mt-kahypar/definitions.h"
#include "robin_hood.h"

namespace mt_kahypar {
namespace ds {


template<typename Map>
class HashMapIterator {
private:
    using Element = typename Map::Element;

public:
    HashMapIterator(const Map& map, const size_t offset = 0) : _map(map), _offset(offset) {}

    const Element& operator*() {
        return _map._entries[_map._positions[_offset]];
    }

    HashMapIterator& operator++ () {
        ++_offset;
        return *this;
    }

    bool operator== (const HashMapIterator& it) {
        return _offset == it._offset;
    }

    bool operator!=(const HashMapIterator& it) {
        return !(*this == it);
    }
private:
    const Map& _map;
    size_t _offset;
};

template<typename Key, typename Value>
class RobinHoodMap {
private:
    robin_hood::unordered_map<Key, Value> map;

public:
    using Element = std::pair<Key, Value>;
    using Iterator = typename robin_hood::unordered_map<Key, Value>::iterator;

    RobinHoodMap(size_t size) {
        map.reserve(size);
    }
    // ! returns true if the element is already in the hashmap and increments its value.
    // ! returns false if not
    bool increment(const Key key) {
        auto it = map.find(key);
        if (it != map.end()) {
            it->second += 1;
            return true;
        }
        return false;
    }

    // ! returns true if the element is already in the hashmap and decrements its value.
    // ! erases the Element if its value becomes zero
    // ! returns false if not
    bool decrement(const Key key) {
        auto it = map.find(key);
        if (it != map.end()) {
            it->second -= 1;
            if (it->second == 0) {
                map.erase(it);
            }
            return true;
        }
        return false;
    }

    void insert(Key k, Value v) {
        map.insert({ k,v });
    }

    Iterator begin() {
        return map.begin();
    }
    Iterator end() {
        return map.end();
    }

    void clear() {
        map.clear();
    }
};

template<typename Key, typename Value, typename Hash>
class HashMap {
private:
    using MyType = HashMap<Key, Value, Hash>;

public:
    using Element = std::pair<Key, Value>;
    using Iterator = HashMapIterator<MyType>;

public:
    HashMap(const size_t size) :
        _empty_element(std::numeric_limits<Key>::max()),
        _size(std::max(2, static_cast<int>(pow(2, floor(log2(size)))))),
        _hash(),
        _entries(_size * 2, std::make_pair(_empty_element, Value())) {
        _positions.reserve(_size);
    }

    // ! amount of occupied spaces
    size_t size() const {
        return _positions.size();
    }

    // ! returns true if there are no entries
    bool empty() const {
        return _positions.size() == 0;
    }

    // ! returns true if the key is in the hashmap already
    bool contains(const Key key) const {
        return _entries[getPosition(key)].first == key;
    }

    // ! inserts the key value pair into the hashmap, might trigger a resize
    void insert(const Key& key, const Value& v) {
        insertAlgorithm(key, v);
    }

    // ! inserts the key value pair into the hashmap, might trigger a resize
    void insert(const Element& e) {
        insertAlgorithm(e.first, e.second);
    }

    // ! swaps the contents of two HashMaps
    void swap(MyType& other) {
        std::swap(_size, other._size);
        _entries.swap(other._entries);
        _positions.swap(other._positions);
        std::swap(_empty_element, other._empty_element);
    }

    // ! returns true if the element is already in the hashmap and increments its value.
    // ! returns false if not
    bool increment(const Key key) {
        size_t pos = getPosition(key);
        if (_entries[pos].first == key) {
            ++_entries[pos].second;
            return true;
        }
        return false;
    }

    // ! returns true if the element is already in the hashmap and decrements its value.
    // ! erases the Element if its value becomes zero
    // ! returns false if not
    bool decrement(const Key key) {
        const size_t pos = getPosition(key);
        if (_entries[pos].first == key) {
            const Value v = --_entries[pos].second;
            if (v == 0) {
                eraseAtPosition(pos);
            }
            return true;
        }
        return false;
    }

    // ! removes the Element with the given key from the Hashmap
    void erase(const Key key) {
        const size_t pos = getPosition(key);
        if (_entries[pos].first == key) {
            eraseAtPosition(pos);
        }
    }

    // ! returns the Key, Value pair with the given key. If it was not inserted before returns a placeholder element.
    Element get(const Key key) const {
        size_t pos = getPosition(key);
        if (_entries[pos].first == key) {
            return _entries[pos];
        }
        return std::make_pair(_empty_element, Value());
    }

    Iterator begin() const {
        return Iterator(*this, 0);
    }

    Iterator end() const {
        return Iterator(*this, size());
    }

    void clear() {
        clearEntries();
        _positions.clear();
    }

private:

    // TODO: This is O(n), improve?
    // ! removes a given position from the list of positions
    void removePosition(size_t pos) {
        for (size_t i = 0; i < _positions.size(); ++i) {
            if (_positions[i] == pos) {
                std::swap(_positions[i], _positions[_positions.size() - 1]);
                _positions.pop_back();
                return;
            }
        }
        ASSERT(false, "Trying to remove a position which is not occupied");
    }

    // ! erases the element at the given position
    // ! expects a valid entry at that position
    void eraseAtPosition(size_t pos) {
        ASSERT(_entries[pos].first != _empty_element, "trying to erase an empty element");
        ASSERT(pos < _entries.size() - 1, "trying to erase a element past the table at position");
        // erasing element
        _entries[pos].first = _empty_element;
        size_t j = pos;
        // fixing the invariant
        while (_entries[++j].first != _empty_element) {
            size_t h = _hash(_entries[j].first);
            size_t s = _size - 1;
            size_t result = h & s;
            if (result <= pos) {
                std::swap(_entries[pos], _entries[j]);
                pos = j;
            }
        }
        removePosition(pos);
    }

    // ! actual implementation of the insertion algorithm
    void insertAlgorithm(const Key key, const Value v) {
        // if half full, resize
        if (_positions.size() >= _size) {
            resize();
        }
        size_t pos = getPosition(key);
        _entries[pos] = std::make_pair(key, v);
        _positions.push_back(pos);
    }

    // ! resizes the hashmap by creating a new one and inserting each element again
    void resize() {
        _size *= 2;
        MyType other(_size);
        ASSERT(other._entries.size() > _entries.size(), "resize doesn't increase the size of the table");
        for (const size_t pos : _positions) {
            other.insert(_entries[pos]);
        }
        swap(other);
    }

    // finds a free position
    size_t getPosition(const Key key) const {
        size_t pos = _hash(key) & (_size - 1);
        // there should always be empty_elements in the table
        while (!(_entries[pos].first == _empty_element || _entries[pos].first == key)) {
            ++pos;
            ASSERT(pos < _entries.size(), "get Position runs out of the table");
        }
        return pos;
    }

    void clearEntries() {
        for(const size_t p : _positions) {
            _entries[p].first = _empty_element;
        }
    }

    friend Iterator;
    Key _empty_element;
    size_t _size;
    Hash _hash;
    std::vector<Element> _entries;
    std::vector<size_t> _positions;
};
}
}