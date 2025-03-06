#include "ir/Permutation.hpp"

#include "Definitions.hpp"
#include "ir/operations/Control.hpp"

#include <algorithm>
#include <cassert>

namespace qc {
[[nodiscard]] auto Permutation::apply(const Controls& controls) const
    -> Controls {
  if (empty()) {
    return controls;
  }
  Controls c{};
  for (const auto& control : controls) {
    c.emplace(at(control.qubit), control.type);
  }
  return c;
}
[[nodiscard]] auto Permutation::apply(const Targets& targets) const -> Targets {
  if (empty()) {
    return targets;
  }
  Targets t{};
  for (const auto& target : targets) {
    t.emplace_back(at(target));
  }
  return t;
}

[[nodiscard]] auto Permutation::apply(const Qubit qubit) const -> Qubit {
  if (empty()) {
    return qubit;
  }
  return at(qubit);
}

[[nodiscard]] auto Permutation::maxKey() const -> Qubit {
  if (empty()) {
    return 0;
  }
  return crbegin()->first;
}

[[nodiscard]] auto Permutation::maxValue() const -> Qubit {
  if (empty()) {
    return 0;
  }
  return std::max_element(
             cbegin(), cend(),
             [](const auto& a, const auto& b) { return a.second < b.second; })
      ->second;
}

[[nodiscard]] int Permutation::findPmtIndex(Qubit q) const {
  for(auto &it: *this) {
    if(it.second == q) {
      return static_cast<int>(it.first);
    }
  }
  return -1;
}

[[nodiscard]] std::array<int, MAX_QUBITS_NUMBER> Permutation::generatePreIndex() {
  std::array<int, MAX_QUBITS_NUMBER> res{};
  int prev = -1;
  for(auto &it: *this) {
    assert(it.first >= 0 && it.first < MAX_QUBITS_NUMBER);
    res.at(it.first) = prev;
    prev = static_cast<int>(it.first);
  }
  return res;
}

} // namespace qc
