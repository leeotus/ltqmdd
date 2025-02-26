#include "dd/DDSifting.hpp"
#include "dd/DDDebug.hpp"
#include "dd/DDLinear.hpp"

namespace dd {

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

// TODO: We need to figure out those nodes that point to the level "adj" and
// then apply "sifting" algorithms on them.
static void __check_and_sifting_bf(Package<>* dd, int adj, const Permutation& pmt) {
  auto pmtlvl = pmt.findPmtLevel(adj);
  for(auto itm = pmt.rend(); itm!=pmt.rbegin(); --itm) {
    if(itm->first > pmtlvl+1) {
      auto nodes = dd->mUniqueTable.getTableColumn(itm->second);
      for(auto &ptr : nodes) {
        if(ptr->ref == 0) {
          continue;
        }
        auto es = ptr->e;
        for(auto i=0;i<NEDGE;++i) {
          if(es[i].p->v == adj) {
            dd->decRefOnly(es[i]);
            // TODO: Do single sifting algorithm
            es[i].p = __single_skipped_sifting(dd, es[i].p->e, pmtlvl-1);
            dd->incRef(es[i]);
          }
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
 */
static void __lvl_sifting(mNode *node, Package<>* dd, int curPmtIndex, const Permutation* pmt) {
  if(curPmtIndex == 0) {
    DEBUG_ERROR("sifting with the lowest level varibles!");
    return;
  }
  // FIXME: "adj" is qubitIndex type, sentences like "node->e[i].p->v ==(!=) adj" is wrong
  // NOTE: sifting procedure, basically the same as "lvlswap" function in "dd/DDLinear.hpp"
  std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rearrangeEdges{};
  for(size_t i=0;i<NEDGE;++i) {
    auto eiw = node->e[i].w;    // Get the weight of this edge
    if(node->e[i].isTerminal()) {
      for (size_t j = 0; j < NEDGE; ++j) {
        rearrangeEdges[i][j] = (j == 0 || j == 3) ?
          (Edge<mNode>::one()) : (Edge<mNode>::zero());
        rearrangeEdges[i][j].w = dd->cn.lookup(rearrangeEdges[i][j].w * eiw);
      }
    }else if(node->e[i].p->v != curPmtIndex-1) {
      for (size_t j = 0; j < NEDGE; ++j) {
        if(j == 0 || j == 3) {
          rearrangeEdges[i][j] = Edge<mNode>::one();
          rearrangeEdges[i][j].p = node->e[i].p;
        } else {
          rearrangeEdges[i][j] = Edge<mNode>::zero();
        }
        rearrangeEdges[i][j].w = dd->cn.lookup(rearrangeEdges[i][j].w * eiw);
      }
    } else {
      for (size_t j = 0; j < NEDGE; ++j) {
        rearrangeEdges[i][j] = node->e[i].p->e[j];
        rearrangeEdges[i][j].w = dd->cn.lookup(node->e[i].p->e[j].w * eiw);
      }
    }
    node->e[i].w = (!node->e[i].w.exactlyZero()) ? dd->cn.lookup(Complex::one())
                                               : dd->cn.lookup(Complex::zero());
  }

  for(size_t i = 0; i<NEDGE; ++i) {
    auto *__node = dd->mMemoryManager.get();
    assert(__node->ref == 0);
    __node->v = static_cast<Qubit>(curPmtIndex-1);
    __node->flags = 0;

    for(size_t j=0; j < NEDGE; ++j) {
      __node->e[j] = rearrangeEdges[j][i];
    }

    auto __edge = Edge<mNode>::normalize(__node, __node->e, dd->mMemoryManager, dd->cn);
    // NOTE: May need to apply reduction rules on this edge
    if(!__edge.isTerminal()) {
      const auto &es = __edge.p->e;
      if ((es[0].p == es[3].p) &&
          (es[0].w.exactlyOne() && es[1].w.exactlyZero() &&
           es[2].w.exactlyZero() && es[3].w.exactlyOne())) {
        auto* ptr = es[0].p;
        dd->mMemoryManager.returnEntry(__edge.p);
        node->e[i].p = ptr;
        node->e[i].w = __edge.w;
        continue;
      }
    }
    __edge.p = dd->mUniqueTable.lookup(__edge.p);

    // TODO: Need to verify if the following code works successfully?
    if(node->e[i].isTerminal()) {
      node->e[i] = __edge;
    } else {
      dd->decRef(node->e[i]);
      node->e[i] = __edge;
    }

    if(!node->e[i].isTerminal()) {
      dd->incRef(node->e[i]);
    }
  }

  // FIXME: It may successfully find out the node which has already been in "mUniqueTable"
  // Then, it needs to return this node back to the corresponding memory manager but with
  // "node->ref" not equals to zero!
  node = dd->mUniqueTable.lookup(node);
}

void sifting(Qubit qubitIndex, Package<> *dd, qc::QuantumComputation *qc, bool ori)
{
  auto pmtlvl = qc->initialLayout.findPmtLevel(qubitIndex);
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

  // TODO: After each step of dynamic reordering, the permutation "initialLayout"
  // should stores the new permutation.
  // NOTE: TODO: the final permutation after dynamic reordering
  // should be stored in the "qc->outputPermuation"

  for(auto bucket = 0; bucket < table.size(); ++bucket) {
    auto *node = table[bucket];
    while(node != nullptr) {
      auto *next = node->next;
      if(node->ref != 0) {
        __lvl_sifting(node, dd, pmtlvl, &qc->initialLayout);
      }
      node = next;
    }
  }

  // __check_and_sifting_bf(dd, adj, qc->initialLayout); // FIXME: Bad function

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
