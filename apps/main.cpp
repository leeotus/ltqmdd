/**
 * @file main.cpp
 * @author leeotus (leeotus@163.com)
 * @brief Tests for new dynamic reordering algorithms (defined in "dd/DDSifting.cpp")
 * @note 代码研究使用
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

int main(int argc, char** argv) {
  if(argc != 2)
  {
    std::cout << "Usage: " << static_cast<std::string>(argv[0]) << " <filename>\r\n";
    return 0;
  }
  std::string fileName = argv[1];
  qc::QuantumComputation qc(fileName);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto functionality = dd::buildFunctionality(&qc, *ddpackPtr);

  auto initailDDsize = functionality.size();
  std::cout << "initial dd's size:" << initailDDsize   << "\r\n";

  // 计算时间:
  clock_t start = 0;
  clock_t finish = 0;
  double totalTime{0};

  start = clock();
  // 进行sifting算法:
  dd::DDSiftingAux(functionality, ddpackPtr.get(), &qc);

  finish = clock();

  auto finalSize = functionality.size();
  totalTime = (double)(finish-start) / CLOCKS_PER_SEC;

  std::cout << "total time: " << totalTime << "s, \t";
  std::cout << "final dd's size:" << finalSize << "\r\n";

  return 0;
}
