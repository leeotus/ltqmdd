#include "dd/DDReorder.hpp"
#include "dd/Export.hpp"

#include <cassert>
#include <cstring>
#include <unordered_map>

namespace dd {

void recordOptimalState(OptimalState *state, Qubit level, ReorderScheme scheme, bool up)
{
    assert(state != nullptr);
    state->optimalLevel = level;
    state->scheme = scheme;
    state->up = up;
}

bool ReorderStepManager::isLinkerAvail()
{
    return (freeLinker != nullptr);
}

bool ReorderStepManager::isCurChunkAvail() const
{
    return (cursor < curChunkSize);
}

void ReorderStepManager::reallocChunk()
{
    auto chunkSize = static_cast<size_t>(INITIAL_STEP_SIZE * STEP_GROWTH_FACTOR);
    auto *chunk = new std::vector<ReorderStep>(chunkSize);
    cursor = 0;
    curChunkSize = chunkSize;
    stepool.push_back(chunk);
}

ReorderStep* ReorderStepManager::get()
{
    if(isLinkerAvail())
    {
        ReorderStep *p = freeLinker;
        freeLinker = freeLinker->next;
        return p;
    }
    if(!isCurChunkAvail())
    {
        reallocChunk();
    }
    auto *curChunk = stepool.back();
    ReorderStep *p = &(curChunk->at(cursor)); 
    cursor += 1;
    return p;
}

void ReorderStepManager::returnEntry(ReorderStep **step)
{
    (*step)->next = freeLinker;
    freeLinker = *step;
    *step = nullptr;
}

void VarOrder::printOrder(const std::string &qubitName) const
{
    std::vector<std::string> ordering;
    for(Qubit i=0;i<nqubits;++i)
    {  
        auto index = qtc->outputPermutation[i];
        ordering.push_back(qubitName + std::to_string(index));
    }
    // TODO: 由于目前只有original sifting算法,所以可以直接打印出来,之后需要修改:
    for(auto &q: ordering)
    {
        std::cout << q << " ";
    }
    std::cout<< "\r\n";
    // just for test:
    std::cout << "--------------测试:---------------\r\n";
    auto permutation = qtc->initialLayout;
    for(auto &it : reorderSteps)
    {
        if(it->scheme == SCHEME_SIFTING)
        {
            auto pos = it->level;
            if(it->up)
            {
                auto tmp = permutation[pos+1];
                permutation[pos+1] = permutation[pos];
                permutation[pos] = tmp;
            } else {
                auto tmp = permutation[pos];
                permutation[pos] = permutation[pos-1];
                permutation[pos-1] = tmp;
            }
        }
    }
    for(auto i=0;i<nqubits;++i)
    {
        std::cout << qubitName << permutation[i] << " ";
    }
    std::cout << "\r\n";
}

void VarOrder::dump2graph(std::string &filename)
{
    dd::export2Dot(mdd, filename);
}

void VarOrder::clear()
{
    while(!isRecordEmpty())
    {
        popRecord();
    }
}


} // namespace dd