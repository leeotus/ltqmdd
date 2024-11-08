#pragma once

#include "dd/Edge.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Node.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>
#include <unistd.h>
#include <vector>

#define INITIAL_STEP_SIZE 100
#define STEP_GROWTH_FACTOR 2

namespace dd {

using ReorderScheme = enum {
  SCHEME_NONE,
  SCHEME_SIFTING,
  SCHEME_LTRANS_LOWER,
  SCHEME_LTRANS_UPPER,
  SCHEME_LTRANS_MIXED
};

/**
 * @brief 记录最佳位置和所采用的scheme
 */
struct OptimalState {
  Qubit optimalLevel;
  ReorderScheme scheme;
  size_t minddSize;
  bool up;
};

/**
 * @brief 记录最优筛选位置信息,方便之后进行回溯恢复
 * @param state 记录信息的结构体对象
 * @param level 记录最优位置层
 * @param scheme 记录采用的筛选方案
 * @param up 记录采用的筛选方案是向上筛选还是向下筛选
 */
void recordOptimalState(OptimalState* state, Qubit level, ReorderScheme scheme,
                        bool up);

/**
 * @brief 记录每一步变换
 */
struct ReorderStep {
  Qubit level;          // 要进行交换的层
  ReorderScheme scheme; // 采用哪种算法
  bool up;              // 跟上层交换还是跟下层交换
  size_t ddsize;        // 记录变换之后的dd大小
  ReorderStep* next;    // 仅供内存管理使用
};

class ReorderStepManager {
public:
  explicit ReorderStepManager(size_t initialSize = INITIAL_STEP_SIZE)
      : curChunkSize(initialSize) {
    auto* chunk = new std::vector<ReorderStep>(initialSize);
    stepool.push_back(chunk);
  }

  ReorderStepManager(const ReorderStepManager&) = delete;
  ReorderStepManager& operator=(const ReorderStepManager&) = delete;

  /**
   * @brief 检测当前在使用的内存块是否有空闲节点可供使用
   * @return 有空闲节点可供使用则返回true,否则返回false
   */
  [[nodiscard]] bool isCurChunkAvail() const;

  /**
   * @brief 检测空闲链表中是否有空闲节点可供使用
   * @return 有空闲节点可返回使用则返回true,否则返回false
   */
  inline bool isLinkerAvail();

  /**
   * @brief
   * 如果当前的内存块已经被全部使用完(即空闲链表为空,内存块也为空)时,向操作系统
   * 再次申请另一块内存块放到内存池中
   */
  void reallocChunk();

  /**
   * @brief 获取内存池中的一个ReorderStep对象并返回指向该对象的指针
   */
  [[nodiscard]] ReorderStep* get();

  /**
   * @brief 将申请到的对象返还内存池
   */
  void returnEntry(ReorderStep** step);

  /**
   * @todo 垃圾回收函数
   */
  void garbageCollect();

private:
  size_t cursor{0};
  size_t curChunkSize;
  ReorderStep* freeLinker{nullptr};
  std::vector<std::vector<ReorderStep>*> stepool; // 内存池
};

/**
 * @brief 记录经过不同reordering方案之后的最优变量序和对应的decision diagram大小
 * @note 目前仅考虑适用于mNode系列的函数,暂时还未考虑其他(dNode,vNode)的情况
 * @date 2024/10/25
 */
class VarOrder {
public:
  explicit VarOrder(MatrixDD mdd, qc::QuantumComputation* qc)
      : nqubits(qc->getNqubits()), mdd(mdd), qtc(qc),
        manager(new ReorderStepManager()) {}

  ~VarOrder() {
    qtc = nullptr;
    mdd.p = nullptr;
    mdd.w = Complex::zero();
  }

  // TODO:
  /**
   * @brief 格式化打印变量序
   * @param qubitName 指定量子比特的名字
   */
  void printOrder(const std::string& qubitName) const;

  /**
   * @brief 将decision diagram导出成.dot文件
   * @param fileName 输出文件的文件名
   */
  void dump2graph(std::string& fileName);

  /**
   * @brief 将采用的步骤记录下来以及得到的结果记录下来
   * @param scheme 记录采用的是哪种筛选方案
   * @param ddSize 记录decision diagram的大小
   */
  void record(Qubit level, ReorderScheme scheme, size_t ddSize, bool up) {
    auto* step = manager->get();
    assert(step != nullptr);
    step->level = level;
    step->scheme = scheme;
    step->ddsize = ddSize;
    step->up = up;
    step->next = nullptr;

    reorderSteps.push_back(step);
  }

  /**
   * @brief 取出上一次变换的步骤.
   */
  [[nodiscard]] ReorderStep* lastRecord() {
    if (!isRecordEmpty()) {
      return reorderSteps.back();
    }
    return nullptr;
  }

  /**
   * @brief 撤销上一次的操作记录
   */
  inline void popRecord() {
    assert(reorderSteps.size() != 0);
    auto* step = reorderSteps.back();
    reorderSteps.pop_back();
    manager->returnEntry(&step);
  }

  /**
   * @brief 查看记录是否为空
   */
  inline bool isRecordEmpty() { return reorderSteps.empty(); }

  /**
   * @brief 清空所有记录步骤
   */
  void clear();

  /**
   * @brief 获取目前有多少记录
   */
  inline int size() { return reorderSteps.size(); }

  /**
   * @brief 获取处于i位置的记录
   * @param i
   */
  inline ReorderStep* at(int i) {
    if (i >= 0 && i < this->size()) {
      return reorderSteps.at(i);
    }
    return nullptr;
  }

private:
  size_t nqubits; // 记录qubit数量
  MatrixDD mdd;   // 保存指向decision diagram根节点的指针
  qc::QuantumComputation* qtc;
  std::vector<ReorderStep*>
      reorderSteps; // 记录每一步进行的哪种变换,以及对应的交换层和交换方式
  ReorderStepManager* manager;
};

} // namespace dd