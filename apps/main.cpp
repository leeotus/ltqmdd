#include "algorithms/BernsteinVazirani.hpp"
#include "algorithms/Entanglement.hpp"
#include "algorithms/Grover.hpp"
#include "algorithms/QFT.hpp"
#include "algorithms/QPE.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"
#include "algorithms/WState.hpp"
#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "dd/Export.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "dd/Simulation.hpp"
#include "dd/statistics/PackageStatistics.hpp"
#include "ir/QuantumComputation.hpp"

#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>

#include "dd/DDCompletement.hpp"
#include "dd/DDLinear.hpp"
#include "dd/DDReorder.hpp"

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

  // 补全dd决策图
  dd::levelCompleteSkipped(functionality, ddpackPtr.get());

  auto afterDDsize = functionality.size();
  std::cout << "initial dd size :" << initailDDsize << "\r\n";
  std::cout << "after completing, dd size :" << afterDDsize << "\r\n";

  auto *vo = new dd::VarOrder(functionality, &qc);
  dd::reorderSelect(functionality, ddpackPtr.get(), &qc, dd::SCHEME_LTRANS_MIXED, vo);
  std::cout << "dd's size after original sifting:" << functionality.size() << "\r\n";
  // dd::checkRefValue(functionality);

  // 第一次筛选后的dd大小:
  auto curddSize = functionality.size();
  size_t cycleddSize{};
  int cnt = 10;   // 原本是12,现在尝试将该值改得更大
  for(int i=0,j=0;i<100;++i)
  {
    dd::reorderSelect(functionality, ddpackPtr.get(), &qc, dd::SCHEME_LTRANS_MIXED, vo);
    cycleddSize = functionality.size();
    if(std::abs(cycleddSize - curddSize <= 10))
    {
      j+=1;
      if(j == cnt)
      {
        break;
      }
    }
    if(cycleddSize < curddSize)
    {
      curddSize = cycleddSize;
    }
  }
  std::cout << "final dd's size:" << cycleddSize << "\r\n";
  // dd::checkRefValue(functionality);

  const std::string qubitName = "x";
  vo->printOrder(qubitName);

  return 0;
}


