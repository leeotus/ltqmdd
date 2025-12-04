/**
 * @file DDSifting.hpp
 * @author leeotus
 * @brief Sifting algorithms for QMDD's dynamic reordering
 * @date 2025-02-15
 */
#pragma once

#include "dd/DDReorder.hpp"
#include "dd/Node.hpp"
#include "dd/Edge.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"
#include "dd/DDReorder.hpp"
#include "dd/DDLinear.hpp"
#include "dd/DDDebug.hpp"
#include "dd/DDCommons.hpp"

namespace dd {

// TODO: Sifting algorithms
/**------------------------------------------------------------------------
 * !                              WARNING
 * 目前这系列函数只是根据之前写的DDLinear.hpp中的代码来修改的
 * 不过因为要实现dynamic reordering,所以这里的函数目前还未完善,可能需要将一些额外的
 * 参数输入,比如:每一个qubit register的名称,每一个targets, controls(存储在构造
 * dd的过程中的变量的index, NOTE: 由于需要改变变量的序列,所以targets, controls
 * 里面的值也应该随之进行修改)
 *------------------------------------------------------------------------**/

/**
 * @brief Sifting algorithm: exchange the adjacent variable in QMDD
 * @param qbIndex qubit index, defined in "qregs"
 * @param dd Package<>* pointer, manager of DD nodes and etc.
 * @param qc Contains information of QMDD, for example number of qubits
 * the initial and output permuation of variables.
 * @param ori decide the orientation of "Sifing" algorithm
 * @deprecated Use "reducedSifting" function instead.
 */
void sifting(Qubit qbIndex, Package<> *dd, qc::QuantumComputation *qc, bool ori=false);

// TODO: 修改weight检测, 需要引入负数权重(负数权重会使用内存对齐)
template <typename Config>
void reduced_single_sifting(mNode* node, Package<Config>* dd, int curPmtIndex,
                            const Permutation* pmt) {
  if (curPmtIndex == 0) {
    return;
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
      } else if (node->e[i].p == nullptr && node->e[i].w.approximatelyZero()) {
        for (size_t j = 0; j < NEDGE; ++j) {
          rarEdges[i][j] = Edge<mNode>::zero();
        }
      }
    } else if (node->e[i].p->v < curPmtIndex - 1) {
      for (size_t j = 0; j < NEDGE; ++j) {
        if (j == 0 || j == 3) {
          // rarEdges[i][j].w = Complex::one();
          // auto *nodeptr = dd->mMemoryManager.get();
          // assert(nodeptr->ref == 0);
          // nodeptr->e = node->e;
          // rarEdges[i][j].p = nodeptr;

          // RESEARCH: 直接复制进去是否可行?
          rarEdges[i][j] = node->e[i];
          rarEdges[i][j].w = node->e[i].w;
          // dd->incRef(rarEdges[i][j]);

          // NOTE: 暂时去掉下面这行
          // rarEdges[i][j].p = node->e[i].p;
        } else {
          rarEdges[i][j] = Edge<mNode>::zero();
        }
      }
    } else if(node->e[i].p->v == curPmtIndex - 1) {
      for (size_t j = 0; j < NEDGE; ++j) {
        rarEdges[i][j] = node->e[i].p->e[j];
        rarEdges[i][j].w = node->e[i].p->e[j].w.approximatelyZero() ? Complex::zero() : dd->cn.lookup(node->e[i].p->e[j].w * eiw);
      }
    } else if(node->e[i].p->v >= curPmtIndex) {
      // 正常不可能运行到此处
      std::cerr << "equals to current permutation index!\r\n";
    }
    // node->e[i].w =
    //     (!node->e[i].w.exactlyZero()) ? Complex::one() : Complex::zero();
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
          (es[0].w.approximatelyEquals(dd::Complex::one()) && es[1].w.approximatelyZero() &&
           es[2].w.approximatelyZero() && es[3].w.approximatelyEquals(dd::Complex::one()))) {
        auto* ptr = es[0].p;

        if (!node->e[i].isTerminal()) {
          dd->decRef(node->e[i]);
        }
        dd->mMemoryManager.returnEntry(eptr.p);
        node->e[i].p = ptr;
        if (!node->e[i].isTerminal()) {
          dd->incRef(node->e[i]);
        }
        node->e[i].w = dd->cn.lookup(eptr.w);
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
    // NOTE: 测试环境加入了ref值大小的检测
    assert(es.isTerminal() || es.p->ref != 0);
  }
}

/**
 * @brief 自己提出来的更简洁的sifting算法的一种可能实现方式，用于替换上述的sifting函数
 * @param qbIndex qubit index, defined in "qregs"
 * @param dd Package<>* pointer, manager of DD nodes and etc.
 * @param qc Contains information of QMDD, for example number of qubits
 * the initial and output permuation of variables.
 * @param ori decide the orientation of "Sifing" algorithm
 * @note 目前还在测试当中
 */
template <typename Config>
void reducedSifting(Qubit qbIndex, Package<Config> *dd, qc::QuantumComputation *qc, bool ori=false)
{
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
        reduced_single_sifting(node, dd, pmtlvl, &qc->initialLayout);
      }
      node = next;
    }
  }

  if (ori) {
    auto tmp = qc->initialLayout.at(pmtlvl);
    qc->initialLayout.at(pmtlvl) = qc->initialLayout.at(pmtlvl - 1);
    qc->initialLayout.at(pmtlvl - 1) = tmp;
  } else {
    auto tmp = qc->initialLayout.at(pmtlvl - 1);
    qc->initialLayout.at(pmtlvl - 1) = qc->initialLayout.at(pmtlvl);
    qc->initialLayout.at(pmtlvl) = tmp;
  }
}

/**
 * @brief Upper sifting algorithm
 */
void reducedUpper(Qubit index, Package<> *dd, qc::QuantumComputation *qc, bool ori=false);

/**
 * @brief 完整的sifting算法的入口函数
 * @param edge 要做sifting算法的指向根节点的边
 * @param dd 管理decision diagram的数据包
 * @param qc
 * @todo 修改存放变量序的结构
 */
// template<typename Config=dd::DDPackageConfig>
template <typename Config=dd::UnitarySimulatorDDPackageConfig>
// void DDSiftingAux(Edge<mNode> root, Package<Config>* dd, QuantumComputation *qc);
// void DDSiftingUp(Edge<mNode> root, Package<>* dd, QuantumComputation *qc);
// void DDSiftingDown(Edge<mNode> root, Package<>* dd, QuantumComputation *qc);
// template<typename Config=dd::DDPackageConfig>
void DDSiftingAux(Edge<mNode> root, Package<Config>* dd, QuantumComputation* qc, VarOrder *vo) {
  // VarOrder *vo = new VarOrder(qc);
  assert(vo != nullptr);
  size_t n = qc->getNqubits() - 1;
  std::vector<bool> freeLevel(n + 1, true);
  Qubit level{0};

  OptimalState optimalState{}; // 记录最优位置和采用的方案
  optimalState.scheme = SCHEME_NONE;
  for (size_t i = 0; i < n; ++i) {
    auto minSize = root.size();
    uint64_t maxActiveLevel{0};

    for (size_t j = 0; j < n; ++j) {
      auto var = qc->initialLayout[j];
      if (freeLevel.at(var) && dd->active.at(var) > maxActiveLevel) {
        maxActiveLevel = dd->active.at(var);
        level = j;
      }
    }
    freeLevel.at(qc->initialLayout[level]) = false;
    optimalState.optimalLevel = level;

    if (level * 2 < n) {
      auto startPos = level; // 记录开始的位置
      while (level > 0) {
        // 向下筛选
        reducedSifting(level, dd, qc);
        auto ddSize = root.size();

        recordStep(level, SCHEME_SIFTING, ddSize, false, vo);
        if (ddSize < minSize) {
          minSize = ddSize;
          optimalState.optimalLevel = level - 1;
        }
        level -= 1;
      }

      while (level < n) {
        reducedSifting(level, dd, qc, true);

        if (level < startPos) {
          cancelRecord(vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, true, vo);
          if (ddSize < minSize) {
            minSize = ddSize;
            optimalState.optimalLevel = level + 1;
          }
        }
        level += 1;
      }

      while (level > optimalState.optimalLevel) {
        reducedSifting(level, dd, qc);

        if (level > startPos) {
          cancelRecord(vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, false, vo);
        }
        level -= 1;
      }
    } else {
      auto startPos = level;
      while (level < n) {
        reducedSifting(level, dd, qc, true);

        auto ddSize = root.size();

        recordStep(level, SCHEME_SIFTING, ddSize, true, vo);
        if (ddSize < minSize) {
          minSize = ddSize;
          optimalState.optimalLevel = level + 1;
        }
        level += 1;
      }

      while (level > 0) {
        reducedSifting(level, dd, qc);

        if (level > startPos) {
          cancelRecord(vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, false, vo);
          if (ddSize < minSize) {
            minSize = ddSize;
            optimalState.optimalLevel = level - 1;
          }
        }
        level -= 1;
      }

      while (level < optimalState.optimalLevel) {
        reducedSifting(level, dd, qc, true);

        if (level < startPos) {
          cancelRecord(vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, true, vo);
        }

        level += 1;
      }
    }
  }
}

} // namespace dd
