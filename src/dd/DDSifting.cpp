#include <vector>
#include "dd/DDSifting.hpp"
#include "dd/DDDebug.hpp"
#include "dd/DDLinear.hpp"

namespace dd {

std::queue<Edge<mNode>> __is_parent(mNode* par, mNode *child) {
  std::queue<Edge<mNode>> res{};
  if(par == nullptr || par->ref == 0 || par->v <= child->v) {
    return res;
  }
  const auto &es = par->e;
  for(int i=0;i<NEDGE;++i) {
    if(es[i].p == child) {
      res.push(es[i]);
    }
  }
  return res;
}

static void __reduce_from_parents(mNode* nodeptr) {
  auto& es = nodeptr->e;
  while ((es[0].p == es[3].p) && (es[0].w.exactlyOne() && es[1].w.exactlyZero() &&
                               es[2].w.exactlyZero() && es[3].w.exactlyOne())) {
    // TODO: Find out the parents of this node, then reduce this node
    for (auto parptr : nodeptr->parents) {
      auto res = __is_parent(parptr.second, nodeptr);
      if(res.empty()) {
        continue;
      }

      __reduce_from_parents(parptr.second);
    }
  }
}

static mNode* __single_skipped_sifting(Package<> *dd, const std::array<Edge<mNode>, NEDGE> es, int index) {
  auto *__node = dd->mMemoryManager.get();
  assert(__node->ref == 0);
  __node->v = index;
  __node->flags = 0;
  for(auto i=0;i<NEDGE;++i) {
    __node->e[i] = es[i];
  }
  __node = dd->mUniqueTable.lookup(__node);
  return __node;
}

/**
 * FIXME: TODO: For those reduced nodes, it's necessary to find out thier "parents" and do the sifting algo on them.
 * @param pmtlvl The subsequent permutation level of the adjacent variables.
 */
static void __checkpar_and_sifting(Package<>* dd, int adjPmtlvl, const Permutation& pmt) {
  // Take out all of the nodes from UniqueTable:
  auto nodes = dd->mUniqueTable.getTableColumn(adjPmtlvl);
  for(auto* node : nodes) {
    if((node != nullptr) && node->ref != 0) {
      for(auto &par : node->parents) {
        auto res = __is_parent(par.second, node);
        while(!res.empty()) {    // Successfully find out the parents which point to a reduced node
          // RESEARCH: Do the sifting algo.
          auto r = res.front();

          std::array<Edge<mNode>, NEDGE> es{};
          for(int i=0;i<NEDGE;++i) {
            es[i] = Edge<mNode>();
          }
          for(int j=0;j<NEDGE;++j) {
            if(node->e[j].w.exactlyOne()) {
              // Directly reassign es[j].p to the child node.
              es[j].p = node->e[j].p;
              es[j].w = Complex::one();
            } else {
              // Allocate a new node first
              auto *newNode = dd->mMemoryManager.get();
              assert(newNode->ref == 0);
              newNode->v = adjPmtlvl;
              for(int k=0;k<NEDGE;++k) {
                if(k == 0 || k == 3) {
                  newNode->e[k].p = node->e[j].p;
                  newNode->e[k].w = node->e[j].w;
                } else {
                  newNode->e[k].p = nullptr;
                  newNode->e[k].w = Complex::zero();
                }
              }
              es[j] = Edge<mNode>::normalize(newNode, newNode->e, dd->mMemoryManager, dd->cn);
              es[j].p = dd->mUniqueTable.lookup(es[j].p);
            }
          }
          auto *newpar = dd->mMemoryManager.get();
          assert(newpar->ref == 0);
          newpar->v = adjPmtlvl + 1;
          for(int i=0;i<NEDGE;++i) {
            newpar->e[i] = es[i];
            if(es[i].p != nullptr) {
              es[i].p->parents[newpar->id] = newpar;
            }
          }
          r.p = newpar;

          res.pop();
        }
      }
    }
  }
}

/**
 * @brief sifing with the lower variables
 * @param node selected nodes
 * @param dd dd manager
 * @param adj expected adjacent varibles' index, -1 means the current nodes
 * are placed in the lowest level, then it's no need to apply sifting algorithm
 * @todo Need to change the parents filed of mNode.
 */
static void __lvl_sifting(mNode *node, Package<>* dd, int curPmtIndex, const Permutation* pmt) {
  if(curPmtIndex == 0) {
    return;
  }
  // NOTE: sifting procedure, basically the same as "lvlswap" function in "dd/DDLinear.hpp"
  std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rarEdges{};
  for(size_t i=0;i<NEDGE;++i) {
    auto eiw = node->e[i].w;    // Get the weight of this edge
    if(node->e[i].isTerminal()) {
      for (size_t j = 0; j < NEDGE; ++j) {
        rarEdges[i][j] = (j == 0 || j == 3) ?
          (Edge<mNode>::one()) : (Edge<mNode>::zero());
        rarEdges[i][j].w = dd->cn.lookup(rarEdges[i][j].w * eiw);
      }
    }else if(node->e[i].p->v != curPmtIndex-1) {
      for (size_t j = 0; j < NEDGE; ++j) {
        if(j == 0 || j == 3) {
          rarEdges[i][j].w = Complex::one();
          rarEdges[i][j].p = node->e[i].p;
        } else {
          rarEdges[i][j] = Edge<mNode>::zero();
        }
        rarEdges[i][j].w = dd->cn.lookup(rarEdges[i][j].w * eiw);
      }
    } else {
      for (size_t j = 0; j < NEDGE; ++j) {
        rarEdges[i][j] = node->e[i].p->e[j];
        rarEdges[i][j].w = dd->cn.lookup(node->e[i].p->e[j].w * eiw);
      }
    }
    node->e[i].w = (!node->e[i].w.exactlyZero()) ? dd->cn.lookup(Complex::one())
                                               : dd->cn.lookup(Complex::zero());
  }

  for(size_t i = 0; i<NEDGE; ++i) {
    auto *nodeptr = dd->mMemoryManager.get();
    assert(nodeptr->ref == 0);
    nodeptr->v = static_cast<Qubit>(curPmtIndex-1);
    nodeptr->flags = 0;

    for(size_t j=0; j < NEDGE; ++j) {
      nodeptr->e[j] = rarEdges[j][i];
      if(!nodeptr->e[j].isTerminal()) {
        nodeptr->e[j].p->parents[nodeptr->id] = nodeptr;
      }
    }

    auto eptr = Edge<mNode>::normalize(nodeptr, nodeptr->e, dd->mMemoryManager, dd->cn);
    // NOTE: May need to apply reduction rules on this edge
    if(!eptr.isTerminal()) {
      const auto &es = eptr.p->e;
      if ((es[0].p == es[3].p) &&
          (es[0].w.exactlyOne() && es[1].w.exactlyZero() &&
           es[2].w.exactlyZero() && es[3].w.exactlyOne())) {
        auto* ptr = es[0].p;
        dd->mMemoryManager.returnEntry(eptr.p);
        node->e[i].p = ptr;
        node->e[i].w = eptr.w;
        if(ptr != nullptr) {
          ptr->parents[node->id] = node;
        }
        continue;
      }
    }
    eptr.p = dd->mUniqueTable.lookup(eptr.p);

    if(node->e[i].isTerminal()) {
      node->e[i] = eptr;
    } else {
      auto tmp = node->e[i];
      node->e[i] = eptr;
      dd->decRef(tmp);
    }

    // // Add *node to the parents of eptr.p
    // if (eptr.p != nullptr) {
    //   eptr.p->parents[node->id] = node;
    // }

    if(!node->e[i].isTerminal()) {
      node->e[i].p->parents[node->id] = node;
      dd->incRef(node->e[i]);
    }
  }

  node = dd->mUniqueTable.lookup(node);
}

void sifting(Qubit qubitIndex, Package<> *dd, qc::QuantumComputation *qc, bool ori)
{
  auto pmtlvl = qc->initialLayout.findPmtIndex(qubitIndex);
  assert(pmtlvl >= 0 && pmtlvl < qc->getNqubits());
  if(ori) {
    pmtlvl = pmtlvl + 1;
  } else {
  }
  if((pmtlvl == qc->getNqubits() && ori) || (pmtlvl == 0 && !ori)) {
    return;
  }

  // get all the variables in the current index
  auto table = dd->mUniqueTable.getTableColumnAndClear(pmtlvl);

  // TODO: After each step of dynamic reordering, the permutation "initialLayout" should stores the new permutation.
  // NOTE: TODO: the final permutation after dynamic reordering should be stored in the "qc->outputPermuation"

  for(auto bucket = 0; bucket < table.size(); ++bucket) {
    auto *node = table[bucket];
    while(node != nullptr) {
      auto *next = node->next;
      // RESEARCH: Throughout the period of building the DD, nodes will be "returned" to the corresponding memoryManager when its "ref == 0".
      // However, this node may have already been pushed into the "UniqueTable", then, its field "v" may be changed the next time when we get a
      // new node from the "UniqueTable".
      if(node->ref != 0 && node->v == pmtlvl) {
        __lvl_sifting(node, dd, pmtlvl, &qc->initialLayout);
        if (node != nullptr) {
          for (auto& e : node->e) {
            if (!e.isTerminal()) {
              e.p->parents[node->id] = node;
            }
          }
        }
      }
      // TODO: Need to check whether *node should be reduced?
      node = next;
    }
  }

  if(pmtlvl != qc->getNqubits()-1) {
    int adjlvl = pmtlvl - 1;
    __checkpar_and_sifting(dd, adjlvl, qc->initialLayout); // FIXME: Bad function
  }

  if(ori) {
    auto tmp = qc->initialLayout.at(pmtlvl);
    qc->initialLayout.at(pmtlvl) = qubitIndex;
    qc->initialLayout.at(pmtlvl-1) = tmp;
  } else {
    auto tmp = qc->initialLayout.at(pmtlvl-1);
    qc->initialLayout.at(pmtlvl-1) = qubitIndex;
    qc->initialLayout.at(pmtlvl) = tmp;
  }
}


// TODO: upper linear sifting algorithm
void upper(Qubit index, Package<> *dd, qc::QuantumComputation *qc, bool ori)
{

}

} // namespace dd
