#pragma once

#include "Definitions.hpp"
#include "operations/Control.hpp"

#include <cstddef>
#include <functional>
#include <array>
#include <map>

#define MAX_QUBITS_NUMBER  32
namespace qc {
class Permutation : public std::map<Qubit, Qubit> {
public:

  [[nodiscard]] auto apply(const Controls& controls) const -> Controls;
  [[nodiscard]] auto apply(const Targets& targets) const -> Targets;
  [[nodiscard]] auto apply(Qubit qubit) const -> Qubit;
  [[nodiscard]] auto maxKey() const -> Qubit;
  [[nodiscard]] auto maxValue() const -> Qubit;

  /**
   * @brief Given a QubitName (named in the input circuit files)'s index
   * and search (and return if successful) its level index in the QMDDs.
   * @return int -1 for errors.
   * @author leeotus
   */
  [[nodiscard]] int findQubitName(Qubit q);

  /**
   * @brief generate "previous index" of this permutation, the index indicates
   * the upper level of each qubit
   * @return std::vector<int> : -1 for none
   * @author leeotus
   */
  [[nodiscard]] std::array<int, MAX_QUBITS_NUMBER> generatePreIndex();

};
}  // namespace qc

// define hash function for Permutation
namespace std {
template <> struct hash<qc::Permutation> {
  std::size_t operator()(const qc::Permutation& p) const {
    std::size_t seed = 0;
    for (const auto& [k, v] : p) {
      qc::hashCombine(seed, k);
      qc::hashCombine(seed, v);
    }
    return seed;
  }
};
} // namespace std
