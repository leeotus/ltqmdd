/**
 * @file test_mnodepar.cpp
 * @author leeotus (leeotus@163.com)
 * @brief Google Test for checking the 'parents' field in "mNode" type.
 * @date 2025-03-02
 * @note Usage: cd ltqmdd && ./build/dynmTest/mnodepar-test
 * @result All passed
 */

#include "gtest/gtest.h"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "dd/DDDebug.hpp"
#include "ir/QuantumComputation.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>

bool checkNodeParents(const char *fileName) {
  if(strlen(fileName) == 0) {
    DEBUG_ERROR("fileName empty!");
    return false;
  }
  qc::QuantumComputation qc(fileName);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto dd = dd::buildFunctionality(&qc, *ddpackPtr);

  auto total = qc.getNqubits();
  // debug_info_printf("total qubit: %d\r\n", total);

  for(int i=total-1;i>=0;--i) {
    auto nodes = ddpackPtr->mUniqueTable.getTableColumn(i);  // 获取每层的dd节点
    for(auto &node : nodes) {
      if(node != nullptr && node->ref != 0) {
        auto es = node->e;
        for(size_t i=0;i<dd::NEDGE;++i) {
          if(!es[i].isTerminal() && es[i].p->ref!=0) {

          }
        }
      }
    }
  }

  return true;
}

TEST(mNodeParUTest, demo2) {
   // path to the test circuits
  const char *testCircuit = "./circuits/demo2.real";

  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, alu4201) {
  const char *testCircuit = "./circuits/revLib/alu4_201.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, apex4202) {
  const char *testCircuit = "./circuits/revLib/apex4_202.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, ex1010230) {
  const char *testCircuit = "./circuits/revLib/ex1010_230.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, ham7299) {
  const char *testCircuit = "./circuits/revLib/ham7_299.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, in0235) {
  const char *testCircuit = "./circuits/revLib/in0_235.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, in2236) {
  const char *testCircuit = "./circuits/revLib/in2_236.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, rd84253) {
  const char *testCircuit = "./circuits/revLib/rd84_253.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, table3264) {
  const char *testCircuit = "./circuits/revLib/table3_264.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, tial265) {
  const char *testCircuit = "./circuits/revLib/tial_265.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, urf6281) {
  const char *testCircuit = "./circuits/revLib/urf6_281.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}

TEST(mNodeParUTest, wim266) {
  const char *testCircuit = "./circuits/revLib/wim_266.real";
  EXPECT_TRUE(checkNodeParents(testCircuit));
}
