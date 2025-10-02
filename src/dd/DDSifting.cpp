#include <vector>
#include <algorithm>
#include "dd/DDSifting.hpp"
#include "dd/DDDebug.hpp"
#include "dd/DDLinear.hpp"

namespace dd {

/**
 * @brief 化简之后的单步sifting算法
 */
static void reduced_single_sifting(mNode* node, Package<>* dd, int curPmtIndex,
                                   const Permutation* pmt) {
  if (curPmtIndex == 0) {
    return;
  }

  // for debug:
  for (auto& es : node->e) {
    assert(es.isTerminal() || es.p->ref != 0);
  }

  // 检测该点的四条出边是不是都是skipped
  auto const check = [curPmtIndex](const Edge<mNode>& e) {
    return e.isTerminal() || e.p->v != curPmtIndex - 1;
  };
  if (std::all_of(std::begin(node->e), std::end(node->e), check)) {
    node = dd->mUniqueTable.lookup(node);
    return;
  }
  std::array<std::array<Edge<mNode>, NEDGE>, NEDGE>
      rarEdges{}; // 保存需要重新分配的边
  // 保存需要重新分配的边
  for (size_t i = 0; i < NEDGE; ++i) {
    auto eiw = node->e[i].w;
    if (node->e[i].isTerminal()) {
      if (node->e[i].isOneTerminal()) {
        for (size_t j = 0; j < NEDGE; ++j) {
          rarEdges[i][j] =
              (j == 0 || j == 3) ? (Edge<mNode>::one()) : (Edge<mNode>::zero());
        }
      } else if (node->e[i].isZeroTerminal()) {
        for (size_t j = 0; j < NEDGE; ++j) {
          rarEdges[i][j] = Edge<mNode>::zero();
        }
      }
    } else if (node->e[i].p->v != curPmtIndex - 1) {
      for (size_t j = 0; j < NEDGE; ++j) {
        if (j == 0 || j == 3) {
          // RESEARCH: 分配新的内存放入到rarEdges数组,而不是单纯的复制
          // rarEdges[i][j].w = Complex::one();
          // auto *nodeptr = dd->mMemoryManager.get();
          // assert(nodeptr->ref == 0);
          // nodeptr->e = node->e;
          // rarEdges[i][j].p = nodeptr;

          // RESEARCH: 直接复制进去是否可行?
          rarEdges[i][j] = node->e[i];

          // NOTE: 暂时去掉下面这行
          // rarEdges[i][j].p = node->e[i].p;
        } else {
          rarEdges[i][j] = Edge<mNode>::zero();
        }
        rarEdges[i][j].w = dd->cn.lookup(rarEdges[i][j].w * eiw);
      }
    } else {
      for (size_t j = 0; j < NEDGE; ++j) {
        rarEdges[i][j] = node->e[i].p->e[j];
        rarEdges[i][j].w = dd->cn.lookup(node->e[i].p->e[j].w * eiw);
      }
    }
    node->e[i].w = (!node->e[i].w.exactlyZero())
                       ? dd->cn.lookup(Complex::one())
                       : dd->cn.lookup(Complex::zero());
  }

  // 重新分配边:
  for (size_t i = 0; i < NEDGE; ++i) {
    auto* nodeptr = dd->mMemoryManager.get();
    assert(nodeptr->ref == 0);
    nodeptr->v = static_cast<Qubit>(curPmtIndex - 1);
    nodeptr->flags = 0;

    for (size_t j = 0; j < NEDGE; ++j) {
      nodeptr->e[j] = rarEdges[j][i];
    }

    auto eptr =
        Edge<mNode>::normalize(nodeptr, nodeptr->e, dd->mMemoryManager, dd->cn);
    if (!eptr.isTerminal()) {
      const auto& es = eptr.p->e;
      if ((es[0].p == es[3].p) &&
          (es[0].w.exactlyOne() && es[1].w.exactlyZero() &&
           es[2].w.exactlyZero() && es[3].w.exactlyOne())) {
        auto* ptr = es[0].p;

        if (!node->e[i].isTerminal()) {
          dd->decRef(node->e[i]);
        }
        dd->mMemoryManager.returnEntry(eptr.p);
        node->e[i].p = ptr;
        if (!node->e[i].isTerminal()) {
          dd->incRef(node->e[i]);
        }
        node->e[i].w = eptr.w;
        if (i == NEDGE - 1) {
          // 需要lookup:
          goto lookupNode;
        }
        continue;
      }
    }

    if (eptr.p) {
      auto res = dd->mUniqueTable.searchUp(eptr.p);
      // eptr.p = dd->mUniqueTable.lookup(eptr.p);
      dd->incRef(eptr);
    }

    if (node->e[i].isTerminal()) {
      node->e[i] = eptr;
    } else {
      dd->decRef(node->e[i]);
      node->e[i] = eptr;
    }
  }

lookupNode:
  node = dd->mUniqueTable.lookup(node);

  // for debug:
  for (auto& es : node->e) {
    assert(es.isTerminal() || es.p->ref != 0);
  }
}

// RESEARCH: 目前不知道要怎么处理upper和lower算法
// ERROR!
// static void upper_single_sifting(mNode *node, Package<>*dd, int curPmtIndex, const Permutation *pmt) {
//   if(curPmtIndex == 0) {
//     return;
//   }

//   // DEBUG:
//   for (auto &es : node->e) {
//     assert(es.isTerminal() || es.p->ref != 0);
//   }

//   std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rarEdges{};

//   size_t row = 0;
//   for (size_t i = 0; i < NEDGE; ++i) {
//     auto eiw = node->e[i].w;
//     for (size_t j = 0; j < NEDGE; ++j) {
//       // 先判断这条应该要放在矩阵的哪个位置:
//       row = j ^ i;

//       if (node->e[i].isTerminal()) {
//         rarEdges[row][j] = node->e[i];
//       } else {
//         auto eijw = node->e[i].p->e[j].w;

//         rarEdges[row][j] = node->e[i].p->e[j];
//         if (!eiw.exactlyOne()) {
//           rarEdges[row][j].w = dd->cn.lookup(eiw * eijw);
//         }
//       }
//     }
//     node->e[i].w = dd->cn.lookup(Complex::one());
//   }

//   for (size_t i = 0; i < NEDGE; ++i) {
//     auto* newNode = dd->mMemoryManager.get();
//     assert(newNode->ref == 0);
//     newNode->v = node->v - 1;
//     newNode->flags = 0;

//     for (size_t j = 0; j < NEDGE; ++j) {
//       newNode->e[j] = rarEdges[i][j];
//     }

//     auto newEdge =
//         Edge<mNode>::normalize(newNode, newNode->e, dd->mMemoryManager, dd->cn);
//     newEdge.p = dd->mUniqueTable.lookup(newEdge.p);

//     if (node->e[i].isTerminal()) {
//       node->e[i] = newEdge;
//     } else {
//       dd->decRef(node->e[i]);
//       node->e[i] = newEdge;
//     }

//     if (!node->e[i].isTerminal()) {
//       dd->incRef(node->e[i]);
//     }
//   }
//   node = dd->mUniqueTable.lookup(node);
// }

// void reducedSifting(Qubit qbIndex, Package<>* dd, qc::QuantumComputation* qc,
//                     bool ori) {
//   auto pmtlvl = qbIndex;
//   assert(pmtlvl >= 0 && pmtlvl < qc->getNqubits());
//   if (ori) {
//     pmtlvl = pmtlvl + 1;
//   }

//   if ((pmtlvl == qc->getNqubits() && ori) || (pmtlvl == 0 && !ori)) {
//     return;
//   }

//   auto table = dd->mUniqueTable.getTableColumnAndClear(pmtlvl);
//   for (auto bucket = 0; bucket < table.size(); ++bucket) {
//     auto* node = table[bucket];
//     while (node != nullptr) {
//       auto* next = node->next;
//       if (node->ref != 0 && node->v == pmtlvl) {
//         reduced_single_sifting(node, dd, pmtlvl, &qc->initialLayout);
//       }
//       node = next;
//     }
//   }

//   if (ori) {
//     auto tmp = qc->initialLayout.at(pmtlvl);
//     qc->initialLayout.at(pmtlvl) = qbIndex;
//     qc->initialLayout.at(pmtlvl - 1) = tmp;
//   } else {
//     auto tmp = qc->initialLayout.at(pmtlvl - 1);
//     qc->initialLayout.at(pmtlvl - 1) = qbIndex;
//     qc->initialLayout.at(pmtlvl) = tmp;
//   }
// }

// TODO: upper linear sifting algorithm
void reducedUpper(Qubit qbIndex, Package<>* dd, qc::QuantumComputation* qc,
                  bool ori) {
  auto pmtlvl = qbIndex;
  assert(pmtlvl >= 0 && pmtlvl < qc->getNqubits());
  if (ori) {
    pmtlvl = pmtlvl + 1;
  }

  if ((pmtlvl == qc->getNqubits() && ori) || (pmtlvl == 0 && !ori)) {
    return;
  }

  auto table = dd->mUniqueTable.getTableColumnAndClear(pmtlvl);
  for (auto bucket = 0; bucket < table.size(); ++bucket) {
    auto* node = table[bucket];
    while (node != nullptr) {
      auto* next = node->next;
      if (node->ref != 0 && node->v == pmtlvl) {
        // TODO: 实现单步upper算法步骤
      }
      node = next;
    }
  }
}

/**
 * @brief Linear Sifting算法的入口函数
 */
// template<typename Config=dd::DDPackageConfig>
// void DDSiftingAux(Edge<mNode> root, Package<Config>* dd, QuantumComputation* qc) {
//   VarOrder vo(root, qc);
//   size_t n = qc->getNqubits() - 1;
//   std::vector<bool> freeLevel(n + 1, true);
//   Qubit level{0};

//   OptimalState optimalState{}; // 记录最优位置和采用的方案
//   optimalState.scheme = SCHEME_NONE;
//   for (size_t i = 0; i < n; ++i) {
//     auto minSize = root.size();
//     uint64_t maxActiveLevel{0};

//     for (size_t j = 0; j < n; ++j) {
//       auto var = qc->initialLayout[j];
//       if (freeLevel.at(var) && dd->active.at(var) > maxActiveLevel) {
//         maxActiveLevel = dd->active.at(var);
//         level = j;
//       }
//     }
//     freeLevel.at(qc->initialLayout[level]) = false;
//     optimalState.optimalLevel = level;

//     if (level * 2 < n) {
//       auto startPos = level; // 记录开始的位置
//       while (level > 0) {
//         // 向下筛选
//         reducedSifting(level, dd, qc);
//         auto ddSize = root.size();

//         recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
//         if (ddSize < minSize) {
//           minSize = ddSize;
//           optimalState.optimalLevel = level - 1;
//         }
//         level -= 1;
//       }

//       while (level < n) {
//         reducedSifting(level, dd, qc, true);
//         if (level < startPos) {
//           cancelRecord(&vo);
//         } else {
//           auto ddSize = root.size();
//           recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
//           if (ddSize < minSize) {
//             minSize = ddSize;
//             optimalState.optimalLevel = level + 1;
//           }
//         }
//         level += 1;
//       }

//       while (level > optimalState.optimalLevel) {
//         reducedSifting(level, dd, qc);
//         if (level > startPos) {
//           cancelRecord(&vo);
//         } else {
//           auto ddSize = root.size();
//           recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
//         }
//         level -= 1;
//       }
//     } else {
//       auto startPos = level;
//       while (level < n) {
//         reducedSifting(level, dd, qc, true);
//         auto ddSize = root.size();

//         recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
//         if (ddSize < minSize) {
//           minSize = ddSize;
//           optimalState.optimalLevel = level + 1;
//         }
//         level += 1;
//       }

//       while (level > 0) {
//         reducedSifting(level, dd, qc);

//         if (level > startPos) {
//           cancelRecord(&vo);
//         } else {
//           auto ddSize = root.size();
//           recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
//           if (ddSize < minSize) {
//             minSize = ddSize;
//             optimalState.optimalLevel = level - 1;
//           }
//         }
//         level -= 1;
//       }

//       while (level < optimalState.optimalLevel) {
//         reducedSifting(level, dd, qc, true);

//         if (level < startPos) {
//           cancelRecord(&vo);
//         } else {
//           auto ddSize = root.size();
//           recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
//         }

//         level += 1;
//       }
//     }
//   }
// }

} // namespace dd
