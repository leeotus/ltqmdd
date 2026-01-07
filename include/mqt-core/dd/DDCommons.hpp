#pragma once
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"

namespace dd {

inline bool single_check_weight(mNode* node) {
  if (node == nullptr) {
    return true;
  }
  auto es = node->e;
  for (int i = 0; i < NEDGE; ++i) {
    Edge<mNode> oe = es[i];
    if (oe.w.approximatelyZero() && !oe.w.exactlyZero()) {
      return false;
    }
  }
  return true;
}

inline bool check_weights(MatrixDD root) {
  mNode* n = root.p;
  std::queue<mNode*> nodesq;
  nodesq.push(n);

  while (!nodesq.empty()) {
    int len = nodesq.size();
    for (int i = 0; i < len; ++i) {
      mNode* p = nodesq.front();
      nodesq.pop();

      // 检测4条出边的weight是否正常
      if (!single_check_weight(p)) {
        std::cout << "weight error!\r\n";
        return false;
      }
      auto es = p->e;
      for (auto oe : es) {
        if (oe.p != nullptr) {
          nodesq.push(oe.p);
        }
      }
    }
  }
  return true;
}

inline bool check_weights(mNode* n) {
  std::queue<mNode*> nodesq;
  nodesq.push(n);

  while (!nodesq.empty()) {
    int len = nodesq.size();
    for (int i = 0; i < len; ++i) {
      mNode* p = nodesq.front();
      nodesq.pop();

      // 检测4条出边的weight是否正常
      if (!single_check_weight(p)) {
        std::cout << "weight error!\r\n";
        return false;
      }
      auto es = p->e;
      for (auto oe : es) {
        if (oe.p != nullptr) {
          nodesq.push(oe.p);
        }
      }
    }
  }
  return true;
}

} // namespace dd
