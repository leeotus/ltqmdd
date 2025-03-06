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

bool siftingalgoTest(const char *file) {
  if(strlen(file) == 0) {
    return false;
  }
  qc::QuantumComputation qc(file);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto dd = dd::buildFunctionality(&qc, *ddpackPtr);

  auto initialSize = dd.size();
  auto half = qc.getNqubits() / 2;

  dd::sifting(half, ddpackPtr.get(), &qc);
  auto oncesiftingSize = dd.size();

  debug_info_printf("Initial dd's size:%d; After one sifting: %d", initialSize, oncesiftingSize);

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
