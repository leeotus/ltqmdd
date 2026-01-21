/**
 * @file main.cpp
 * @author leeotus (leeotus@163.com)
 * @brief Tests for new dynamic reordering algorithms (defined in
 * "dd/DDSifting.cpp")
 * @note 代码研究使用
 */

#include "dd/DDCommons.hpp"
#include "dd/DDReorder.hpp"
#include "dd/DDSifting.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"

#include <cmath>
#include <ctime>
#include <iostream>
#include <math.h>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Usage: " << static_cast<std::string>(argv[0])
              << " <filename>\r\n";
    return 0;
  }
  std::string fileName = argv[1];
  qc::QuantumComputation qc(fileName);
  auto ddpackPtr = std::make_unique<dd::Package<>>();
  auto diagram = dd::buildFunctionality(&qc, *ddpackPtr);

  auto initailDDsize = diagram.size();
  std::cout << "initial dd's size:" << initailDDsize << "\r\n";

  // 计算时间:
  clock_t start = 0;
  clock_t finish = 0;
  double totalTime{0};

  start = clock();
  // 进行sifting算法:
  auto* vo = new dd::VarOrder(&qc);
  dd::DDSiftingAux<>(diagram, ddpackPtr.get(), &qc, vo);

  finish = clock();
  delete vo;

  auto finalSize = diagram.size();
  totalTime = (double)(finish - start) / CLOCKS_PER_SEC;

  std::cout << "total time: " << totalTime << "s, \t";
  std::cout << "final dd's size:" << finalSize << "\r\n";

  return 0;
}
