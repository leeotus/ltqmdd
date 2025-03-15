/**
 * @file test_sifting_algo.cpp
 * @author leeotus (leeotus@163.com)
 * @brief Test for new sifting algos.
 * @note Usage: cd ltqmdd && ./build/test/dynmic/siftingalgo-test
 * @result All passed
 */
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"

#include <math.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>

#include <ctime>

#include "dd/DDSifting.hpp"
#include "dd/DDDebug.hpp"
#include "gtest/gtest.h"

bool checkEdgesWeightCorrect(dd::mNode *node) {
  if(node == nullptr) {
    return true;
  }
  auto &es = node->e;
  for(int i=0;i<dd::NEDGE;++i) {
    if(es[i].w.approximatelyZero() || es[i].w.approximatelyEquals(dd::Complex::one())) {
      continue;
    }
    if(es[i].w > dd::Complex::one()) {
      return false;
    }
  }
  return es[0].p != es[3].p || !es[0].w.exactlyOne() || !es[3].w.exactlyOne() ||
      !es[1].w.exactlyZero() || !es[2].w.exactlyZero();
}

bool siftingalgoTest(const char *file) {
  if(strlen(file) == 0) {
    return false;
  }
  qc::QuantumComputation qc(file);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto dd = dd::buildFunctionality(&qc, *ddpackPtr);
  auto total = qc.getNqubits();

  for(int i=total-1;i>=1;--i) {
    dd::sifting(i, ddpackPtr.get(), &qc);
    debug_info_printf("current dd's size = %d", dd.size());

    auto nodes = ddpackPtr->mUniqueTable.getTableColumn(i);  // 获取每层的dd节点
    for(auto &node : nodes) {
      if(node != nullptr && node->ref != 0 && node->v == i) {
        // Firstly, check the weights of the outgoing edges:
        // if(!checkEdgesWeightCorrect(node)) {
        //   debug_error_printf("Wrong weights!");
        //   return false;
        // }

        auto parId = node->id;
        auto es = node->e;
        for(size_t i=0;i<dd::NEDGE;++i) {
          if(!es[i].isTerminal() && es[i].p->ref!=0) {
            if(es[i].p->parents.find(parId) == es[i].p->parents.end()) {
              debug_error_printf("CURRENT LEVEL:%d", i);
              debug_error_printf("Wrong parents field!");
              debug_error_printf("Suppose parent's id:%d, v%d", parId, node->v);
              auto ptr = es[i].p;
              debug_error_printf("Current node:{ id:%d; v:%d; ref:%d }", ptr->id, ptr->v, ptr->ref);
              debug_error_printf("Its parents:");
              for(auto &par : ptr->parents) {
                auto *p = par.second;
                debug_error_printf("{id:%d, v:%d, ref:%d}", p->id, p->v, p->ref);
              }
              return false;
            }
          }
        }
      }
    }
  }

  return true;
}

TEST(SiftAlgo, alu4201) {
  const char *file = "./circuits/revLib/alu4_201.real";
  EXPECT_TRUE(siftingalgoTest(file));
}

TEST(SiftAlgo, apex4202) {
  const char *testCircuit = "./circuits/revLib/apex4_202.real";
  EXPECT_TRUE(siftingalgoTest(testCircuit));
}

TEST(SiftAlgo, ex1010230) {
  const char *testCircuit = "./circuits/revLib/ex1010_230.real";
  EXPECT_TRUE(siftingalgoTest(testCircuit));
}

TEST(SiftAlgo, ham7299) {
  const char *testCircuit = "./circuits/revLib/ham7_299.real";
  EXPECT_TRUE(siftingalgoTest(testCircuit));
}

TEST(SiftAlgo, in0235) {
  const char *testCircuit = "./circuits/revLib/in0_235.real";
  EXPECT_TRUE(siftingalgoTest(testCircuit));
}

TEST(SiftAlgo, in2236) {
  const char *testCircuit = "./circuits/revLib/in2_236.real";
  EXPECT_TRUE(siftingalgoTest(testCircuit));
}
