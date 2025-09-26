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

/**
 * @brief 自己提出来的更简洁的sifting算法的一种可能实现方式，用于替换上述的sifting函数
 * @param qbIndex qubit index, defined in "qregs"
 * @param dd Package<>* pointer, manager of DD nodes and etc.
 * @param qc Contains information of QMDD, for example number of qubits
 * the initial and output permuation of variables.
 * @param ori decide the orientation of "Sifing" algorithm
 * @note 目前还在测试当中
 */
void reducedSifting(Qubit qbIndex, Package<> *dd, qc::QuantumComputation *qc, bool ori=false);

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
template<typename Config=dd::DDPackageConfig>
// void DDSiftingAux(Edge<mNode> root, Package<Config>* dd, QuantumComputation *qc);
// void DDSiftingUp(Edge<mNode> root, Package<>* dd, QuantumComputation *qc);
// void DDSiftingDown(Edge<mNode> root, Package<>* dd, QuantumComputation *qc);
// template<typename Config=dd::DDPackageConfig>
void DDSiftingAux(Edge<mNode> root, Package<Config>* dd, QuantumComputation* qc) {
  VarOrder vo(root, qc);
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

        recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
        if (ddSize < minSize) {
          minSize = ddSize;
          optimalState.optimalLevel = level - 1;
        }
        level -= 1;
      }

      while (level < n) {
        reducedSifting(level, dd, qc, true);
        if (level < startPos) {
          cancelRecord(&vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
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
          cancelRecord(&vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
        }
        level -= 1;
      }
    } else {
      auto startPos = level;
      while (level < n) {
        reducedSifting(level, dd, qc, true);
        auto ddSize = root.size();

        recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
        if (ddSize < minSize) {
          minSize = ddSize;
          optimalState.optimalLevel = level + 1;
        }
        level += 1;
      }

      while (level > 0) {
        reducedSifting(level, dd, qc);

        if (level > startPos) {
          cancelRecord(&vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, false, &vo);
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
          cancelRecord(&vo);
        } else {
          auto ddSize = root.size();
          recordStep(level, SCHEME_SIFTING, ddSize, true, &vo);
        }

        level += 1;
      }
    }
  }
}

} // namespace dd
