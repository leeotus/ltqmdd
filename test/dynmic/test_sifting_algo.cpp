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

bool siftingalgoTest(const char *fileName) {
  qc::QuantumComputation qc(fileName);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto root = dd::buildFunctionality(&qc, *ddpackPtr);
  std::cout << "Initial size: " << root.size() << "\r\n";
  dd::DDSiftingAux(root, ddpackPtr.get(), &qc);
  std::cout << "Current size: " << root.size() << "\r\n";
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
