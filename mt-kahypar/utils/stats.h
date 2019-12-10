/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <mutex>
#include <string>
#include <unordered_map>

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace utils {
class Stats {
  static constexpr bool debug = false;

 public:
  enum class Type : uint8_t {
    BOOLEAN = 0,
    INT32 = 1,
    INT64 = 2,
    FLOAT = 3,
    DOUBLE = 4
  };

  class Stat {
   public:
    explicit Stat(const bool value) :
      _type(Type::BOOLEAN),
      _value_1(value),
      _value_2(0),
      _value_3(0),
      _value_4(0.0),
      _value_5(0.0) { }

    explicit Stat(const int32_t value) :
      _type(Type::INT32),
      _value_1(false),
      _value_2(value),
      _value_3(0),
      _value_4(0.0),
      _value_5(0.0) { }

    explicit Stat(const int64_t value) :
      _type(Type::INT64),
      _value_1(false),
      _value_2(0),
      _value_3(value),
      _value_4(0.0),
      _value_5(0.0) { }

    explicit Stat(const float value) :
      _type(Type::FLOAT),
      _value_1(false),
      _value_2(0),
      _value_3(0),
      _value_4(value),
      _value_5(0.0) { }

    explicit Stat(const double value) :
      _type(Type::DOUBLE),
      _value_1(false),
      _value_2(0),
      _value_3(0),
      _value_4(0.0),
      _value_5(value) { }

    template <typename T>
    void update(const T delta) { }

    void update(const bool value) {
      _value_1 = value;
    }

    void update(const int32_t delta) {
      _value_2 += delta;
    }

    void update(const int64_t delta) {
      _value_3 += delta;
    }

    void update(const float delta) {
      _value_4 += delta;
    }

    void update(const double delta) {
      _value_5 += delta;
    }

    friend std::ostream & operator<< (std::ostream& str, const Stat& stat);

   private:
    Type _type;
    bool _value_1;
    int32_t _value_2;
    int64_t _value_3;
    float _value_4;
    double _value_5;
  };

 public:
  Stats(const Stats&) = delete;
  Stats & operator= (const Stats &) = delete;

  Stats(Stats&&) = delete;
  Stats & operator= (Stats &&) = delete;

  static Stats& instance() {
    static Stats instance;
    return instance;
  }

  void enable() {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    _enable = true;
  }

  void disable() {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    _enable = false;
  }
  template <typename T>
  void add_stat(const std::string& key, const T value) {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    if (_enable) {
      if (_stats.find(key) == _stats.end()) {
        _stats.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(key),
          std::forward_as_tuple(value));
      }
    }
  }

  template <typename T>
  void update_stat(const std::string& key, const T delta) {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    if (_enable) {
      if (_stats.find(key) == _stats.end()) {
        _stats.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(key),
          std::forward_as_tuple(delta));
      } else {
        _stats.at(key).update(delta);
      }
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const Stats& stats);

 private:
  explicit Stats() :
    _stat_mutex(),
    _stats(),
    _enable(true) { }

  std::mutex _stat_mutex;
  std::unordered_map<std::string, Stat> _stats;
  bool _enable;
};

std::ostream & operator<< (std::ostream& str, const Stats::Stat& stat) {
  switch (stat._type) {
    case Stats::Type::BOOLEAN:
      str << std::boolalpha << stat._value_1;
      break;
    case Stats::Type::INT32:
      str << stat._value_2;
      break;
    case Stats::Type::INT64:
      str << stat._value_3;
      break;
    case Stats::Type::FLOAT:
      str << stat._value_4;
      break;
    case Stats::Type::DOUBLE:
      str << stat._value_5;
      break;
    default:
      break;  // UNKNOWN TYPE
  }
  return str;
}

std::ostream & operator<< (std::ostream& str, const Stats& stats) {
  std::vector<std::string> keys;
  for (const auto& stat : stats._stats) {
    keys.emplace_back(stat.first);
  }
  std::sort(keys.begin(), keys.end());

  for (const std::string& key : keys) {
    str << " " << key << "=" << stats._stats.at(key);
  }
  return str;
}
}  // namespace utils
}  // namespace mt_kahypar
