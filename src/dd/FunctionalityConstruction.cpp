#include "dd/FunctionalityConstruction.hpp"

#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"
#include "dd/DDSifting.hpp"
#include "dd/DDReorder.hpp"

#include <cmath>
#include <cstddef>
#include <stack>

namespace dd {
template <class Config>
MatrixDD buildFunctionality(QuantumComputation* qc, Package<Config>& dd) {
  const auto nq = qc->getNqubits();
  if (nq == 0U) {
    return MatrixDD::one();
  }

  qc::Permutation permutation;
  auto e = dd.createInitialMatrix(qc->ancillary);
  static long int sth = 1000;

  VarOrder *vo = new VarOrder(qc);
  for (const auto& op : *qc) {
    permutation = qc->initialLayout;
    // RESEARCH: 经过dynamic reordering之后op指向的targets和controls内的数值可能需要修改
    auto dd1 = getDD(op.get(), dd, permutation);
    auto tmp = dd.multiply(dd1, e);

    dd.incRef(tmp);
    dd.decRef(e);
    e = tmp;

    if(e.size() > sth) {
      vo->clear();
      // for debug:
      std::cout << "当前的DD大小:" << e.size() << " ";
      std::cout << "超过阈值\r\n";
      DDSiftingAux(e, &dd, qc, vo);
      sth *= 2;
      std::cout << "阈值提升，现在阈值:" << sth << ", ";
      std::cout << "dynamic reordering后的DD大小:" << e.size() << "\r\n";
    }

    dd.garbageCollect();
  }
  // correct permutation if necessary
  changePermutation(e, permutation, qc->outputPermutation, dd);
  e = dd.reduceAncillae(e, qc->ancillary);
  e = dd.reduceGarbage(e, qc->garbage);

  return e;
}

template <class Config>
MatrixDD buildFunctionalityRecursive(QuantumComputation* qc,
                                     Package<Config>& dd) {
  if (qc->getNqubits() == 0U) {
    return MatrixDD::one();
  }

  auto permutation = qc->initialLayout;

  if (qc->size() == 1U) {
    auto e = getDD(qc->front().get(), dd, permutation);
    dd.incRef(e);
    return e;
  }

  std::stack<MatrixDD> s{};
  auto depth = static_cast<std::size_t>(std::ceil(std::log2(qc->size())));
  buildFunctionalityRecursive(qc, depth, 0, s, permutation, dd);
  auto e = s.top();
  s.pop();

  // correct permutation if necessary
  changePermutation(e, permutation, qc->outputPermutation, dd);
  e = dd.reduceAncillae(e, qc->ancillary);
  e = dd.reduceGarbage(e, qc->garbage);

  return e;
}

template <class Config>
bool buildFunctionalityRecursive(QuantumComputation* qc,
                                 std::size_t depth, std::size_t opIdx,
                                 std::stack<MatrixDD>& s,
                                 Permutation& permutation,
                                 Package<Config>& dd) {
  // base case
  if (depth == 1U) {
    auto e = getDD(qc->at(opIdx).get(), dd, permutation);
    ++opIdx;
    if (opIdx == qc->size()) { // only one element was left
      s.push(e);
      dd.incRef(e);
      return false;
    }
    auto f = getDD(qc->at(opIdx).get(), dd, permutation);
    s.push(dd.multiply(f, e)); // ! reverse multiplication
    dd.incRef(s.top());
    return (opIdx != qc->size() - 1U);
  }

  // in case no operations are left after the first recursive call nothing has
  // to be done
  const size_t leftIdx =
      opIdx & ~(static_cast<std::size_t>(1U) << (depth - 1U));
  if (!buildFunctionalityRecursive(qc, depth - 1U, leftIdx, s, permutation,
                                   dd)) {
    return false;
  }

  const size_t rightIdx =
      opIdx | (static_cast<std::size_t>(1U) << (depth - 1U));
  const auto success =
      buildFunctionalityRecursive(qc, depth - 1U, rightIdx, s, permutation, dd);

  // get latest two results from stack and push their product on the stack
  auto e = s.top();
  s.pop();
  auto f = s.top();
  s.pop();
  s.push(dd.multiply(e, f)); // ordering because of stack structure

  // reference counting
  dd.decRef(e);
  dd.decRef(f);
  dd.incRef(s.top());
  dd.garbageCollect();

  return success;
}

template MatrixDD buildFunctionality(qc::QuantumComputation* qc,
                                     Package<DDPackageConfig>& dd);
template MatrixDD
buildFunctionality(qc::QuantumComputation* qc,
                   Package<dd::DensityMatrixSimulatorDDPackageConfig>& dd);
template MatrixDD
buildFunctionality(qc::QuantumComputation* qc,
                   Package<dd::StochasticNoiseSimulatorDDPackageConfig>& dd);

template MatrixDD buildFunctionality(qc::QuantumComputation* qc,
                                     UnitarySimulatorDDPackage& dd);

template MatrixDD buildFunctionalityRecursive(qc::QuantumComputation* qc,
                                              Package<DDPackageConfig>& dd);
template bool buildFunctionalityRecursive(qc::QuantumComputation* qc,
                                          const std::size_t depth,
                                          const std::size_t opIdx,
                                          std::stack<MatrixDD>& s,
                                          qc::Permutation& permutation,
                                          Package<DDPackageConfig>& dd);
template MatrixDD buildFunctionalityRecursive(qc::QuantumComputation* qc,
                                              UnitarySimulatorDDPackage& dd);
template bool buildFunctionalityRecursive(qc::QuantumComputation* qc,
                                          const std::size_t depth,
                                          const std::size_t opIdx,
                                          std::stack<MatrixDD>& s,
                                          qc::Permutation& permutation,
                                          UnitarySimulatorDDPackage& dd);
} // namespace dd
