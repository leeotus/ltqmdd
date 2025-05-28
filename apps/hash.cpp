#include "dd/FunctionalityConstruction.hpp"
#include "dd/Package.hpp"
#include "ir/QuantumComputation.hpp"
#include "dd/Edge.hpp"

static std::size_t hash(const dd::mNode* p) {
  static constexpr std::size_t MASK = 32787;
  std::size_t key = 0U;
  for(std::size_t i=0U; i<p->e.size(); ++i) {
    qc::hashCombine(key, std::hash<dd::Edge<dd::mNode>>{}(p->e[i]));
  }
  key &= MASK;
  return key;
}

int main(int argc, char **argv) {
  printf("hello world\r\n");
  auto *node_root_1 = new dd::mNode();
  node_root_1->v = 1;
  auto *node_root_2 = new dd::mNode();
  node_root_2->v = 1;


  auto *subnode_e1 = new dd::mNode();
  subnode_e1->v = 0;
  subnode_e1->ref = 2;
  for(int i=0;i<dd::NEDGE;++i)
  {
    if(i == 0 || i == 3) {
      subnode_e1->e[i] = dd::Edge<dd::mNode>::zero();
    } else {
      subnode_e1->e[i] = dd::Edge<dd::mNode>::one();
    }
  }
  auto *subnode1_e3 = new dd::mNode();
  for(int i=0;i<dd::NEDGE;++i) {
    if(i == 0 || i == 2) {
      subnode1_e3->e[i] = dd::Edge<dd::mNode>::zero();
    } else {
      subnode1_e3->e[i] = dd::Edge<dd::mNode>::one();
    }
  }
  auto *subnode2_e3 = new dd::mNode();
  for(int i=0;i<dd::NEDGE;++i) {
    if(i == 0 || i == 2) {
      subnode1_e3->e[i] = dd::Edge<dd::mNode>::zero();
    } else {
      subnode1_e3->e[i] = dd::Edge<dd::mNode>::one();
    }
  }

  node_root_1->e[0].p = subnode_e1;
  node_root_1->e[0].w = dd::Complex::one();

  node_root_1->e[1] = dd::Edge<dd::mNode>::zero();
  node_root_1->e[2] = dd::Edge<dd::mNode>::zero();

  node_root_1->e[3].p = subnode1_e3;
  node_root_1->e[3].w = dd::Complex::one();

  auto edge_root_1 = dd::Edge<dd::mNode>::zero();
  edge_root_1.p = node_root_1;
  edge_root_1.w = dd::Complex::one();

  node_root_2->e[0].p = subnode_e1;
  node_root_2->e[0].w = dd::Complex::one();

  node_root_2->e[1] = dd::Edge<dd::mNode>::zero();
  node_root_2->e[2] = dd::Edge<dd::mNode>::zero();

  node_root_2->e[3].p = subnode2_e3;
  node_root_2->e[3].w = dd::Complex::one();

  auto edge_root_2 = dd::Edge<dd::mNode>::zero();
  edge_root_1.p = node_root_2;
  edge_root_1.w = dd::Complex::one();

  std::cout << "hash1:" << hash(node_root_1) << "\r\n";
  std::cout << "hash2:" << hash(node_root_2) << "\r\n";

  // now change:
  node_root_2->e[3].p = subnode1_e3;
  std::cout << "now hash2:" << hash(node_root_2) << "\r\n";

  return  0;
}
