#include "dd/MemoryManager.hpp"

#include "dd/Node.hpp"
#include "dd/RealNumber.hpp"

#include <cassert>
#include <cstddef>

namespace dd {

template <typename T> T* MemoryManager<T>::get() {
  if (entryAvailableForReuse()) {
    return getEntryFromAvailableList();
  }

  if (!entryAvailableInChunk()) {
    allocateNewChunk();
  }

  return getEntryFromChunk();
}

template <typename T> void MemoryManager<T>::returnEntry(T* entry) noexcept {
  assert(entry != nullptr);
  assert(entry->ref == 0);		// 必须保证Ref值为0才可以回收
  entry->next = available;
  assert(available == nullptr || available->ref == 0);
  available = entry;
  // if constexpr (std::is_same_v<T, mNode>) {
  //   size_t parId = entry->id;
  //   for(size_t i=0;i<NEDGE;++i) {
  //     auto edge = entry->e[i];
  //     if(!edge.isTerminal())
  //     {
  //       assert(edge.p->parents.find(parId) != edge.p->parents.end());
  //       edge.p->parents[parId] = nullptr;
  //     }
  //   }
  // }
  stats.trackReturnedEntry();
}

template <typename T>
void MemoryManager<T>::reset(const bool resizeToTotal) noexcept {
  available = nullptr;

  auto numAllocations = stats.numAllocations;
  chunks.resize(1U);
  if (resizeToTotal) {
    chunks[0].resize(stats.numAllocated);
    ++numAllocations;
  }

  chunkIt = chunks[0].begin();
  chunkEndIt = chunks[0].end();

  stats.reset();
  stats.numAllocations = numAllocations;
  stats.numAllocated = chunks[0].size();
}

template <typename T>
T* MemoryManager<T>::getEntryFromAvailableList() noexcept {
  assert(entryAvailableForReuse());
  auto* entry = available;
  available = available->next;
  stats.trackReusedEntries();
  return entry;
}

template <typename T> void MemoryManager<T>::allocateNewChunk() {
  assert(!entryAvailableInChunk());
  const auto newChunkSize = static_cast<std::size_t>(
      GROWTH_FACTOR * static_cast<double>(chunks.back().size()));
  chunks.emplace_back(newChunkSize);
  chunkIt = chunks.back().begin();
  chunkEndIt = chunks.back().end();
  ++stats.numAllocations;
  stats.numAllocated += newChunkSize;
}

template <typename T> T* MemoryManager<T>::getEntryFromChunk() noexcept {
  assert(!entryAvailableForReuse());
  assert(entryAvailableInChunk());
  auto* entry = &(*chunkIt);
  ++chunkIt;
  stats.trackUsedEntries();
  return entry;
}

template class MemoryManager<RealNumber>;
template class MemoryManager<vNode>;
template class MemoryManager<mNode>;
template class MemoryManager<dNode>;

} // namespace dd
