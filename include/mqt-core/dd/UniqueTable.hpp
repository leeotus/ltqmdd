#pragma once

#include "Definitions.hpp"
#include "dd/Edge.hpp"
#include "dd/MemoryManager.hpp"
#include "dd/Node.hpp"
#include "dd/statistics/UniqueTableStatistics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <type_traits>
#include <vector>

namespace dd {

/**
 * @brief Data structure for uniquely storing DD nodes
 * @tparam Node class of nodes to provide/store
 * @tparam NBUCKET number of hash buckets to use (has to be a power of two)
 */
template <class Node, std::size_t NBUCKET = 32768> class UniqueTable {

  static_assert(
      std::disjunction_v<std::is_same<Node, vNode>, std::is_same<Node, mNode>,
                         std::is_same<Node, dNode>>,
      "Node type must be one of vNode, mNode, dNode");

public:
  /**
   * @brief The initial garbage collection limit.
   * @details The initial garbage collection limit is the number of entries that
   * must be present in the table before garbage collection is triggered.
   * Increasing this number reduces the number of garbage collections, but
   * increases the memory usage.
   */
  static constexpr std::size_t INITIAL_GC_LIMIT = 131072U;

  /**
   * @brief The default constructor
   * @param nv The number of variables
   * @param manager The memory manager to use for allocating new nodes.
   * @param initialGCLim The initial garbage collection limit.
   */
  explicit UniqueTable(const std::size_t nv, MemoryManager<Node>& manager,
                       std::size_t initialGCLim = INITIAL_GC_LIMIT)
      : nvars(nv), memoryManager(&manager), initialGCLimit(initialGCLim) {
    for (auto& stat : stats) {
      stat.entrySize = sizeof(Bucket);
      stat.numBuckets = NBUCKET;
    }
  }

  void resize(std::size_t nq) {
    nvars = nq;
    tables.resize(nq);
    // TODO: if the new size is smaller than the old one we might have to
    // release the unique table entries for the superfluous variables
    stats.resize(nq);
    for (auto& stat : stats) {
      stat.entrySize = sizeof(Bucket);
      stat.numBuckets = NBUCKET;
    }
  }

  /**
   * @brief The hash function for the hash table.
   * @details The hash function just combines the hashes of the edges of the
   * node. The hash value is masked to ensure that it is in the range
   * [0, NBUCKET - 1].
   * @param p The node to hash.
   * @returns The hash value of the node.
   */
  static std::size_t hash(const Node* p) {
    static constexpr std::size_t MASK = NBUCKET - 1;
    std::size_t key = 0U;
    for (std::size_t i = 0U; i < p->e.size(); ++i) {
      qc::hashCombine(key, std::hash<Edge<Node>>{}(p->e[i]));
    }
    key &= MASK;
    return key;
  }

  /// Get a reference to the table
  [[nodiscard]] const auto& getTables() const { return tables; }

  /// Get a reference to the statistics
  [[nodiscard]] const auto& getStats() const noexcept { return stats; }

  /// Get a reference to individual statistics
  [[nodiscard]] const auto& getStats(const std::size_t idx) const noexcept {
    return stats.at(idx);
  }

  /// Get a JSON object with the statistics
  [[nodiscard]] nlohmann::basic_json<>
  getStatsJson(const bool includeIndividualTables = false) const {
    if (std::all_of(stats.begin(), stats.end(),
                    [](const UniqueTableStatistics& stat) {
                      return stat.peakNumEntries == 0U;
                    })) {
      return "unused";
    }

    UniqueTableStatistics totalStats;
    for (const auto& stat : stats) {
      totalStats.entrySize = std::max(totalStats.entrySize, stat.entrySize);
      totalStats.numBuckets += stat.numBuckets;
      totalStats.numEntries += stat.numEntries;
      totalStats.peakNumEntries += stat.peakNumEntries;
      totalStats.collisions += stat.collisions;
      totalStats.hits += stat.hits;
      totalStats.lookups += stat.lookups;
      totalStats.inserts += stat.inserts;
      totalStats.numActiveEntries += stat.numActiveEntries;
      totalStats.peakNumActiveEntries += stat.peakNumActiveEntries;
      totalStats.gcRuns = std::max(totalStats.gcRuns, stat.gcRuns);
    }

    nlohmann::basic_json<> j;
    j["total"] = totalStats.json();
    if (includeIndividualTables) {
      std::size_t v = 0U;
      for (const auto& stat : stats) {
        j[std::to_string(v)] = stat.json();
        ++v;
      }
    }
    return j;
  }

  /// Get the total number of entries
  [[nodiscard]] std::size_t getNumEntries() const noexcept {
    return std::accumulate(
        stats.begin(), stats.end(), std::size_t{0},
        [](const std::size_t& sum, const UniqueTableStatistics& stat) {
          return sum + stat.numEntries;
        });
  }

  /// Get the total number of active entries
  [[nodiscard]] std::size_t getNumActiveEntries() const noexcept {
    return std::accumulate(
        stats.begin(), stats.end(), std::size_t{0},
        [](const std::size_t& sum, const UniqueTableStatistics& stat) {
          return sum + stat.numActiveEntries;
        });
  }

  /// Get the peak total number of active entries
  [[nodiscard]] std::size_t getPeakNumActiveEntries() const noexcept {
    return std::accumulate(
        stats.begin(), stats.end(), std::size_t{0},
        [](const std::size_t& sum, const UniqueTableStatistics& stat) {
          return sum + stat.peakNumActiveEntries;
        });
  }

  static bool nodesAreEqual(const Node* p, const Node* q) {
    if constexpr (std::is_same_v<Node, dNode>) {
      return (p->e == q->e && (p->flags == q->flags));
    } else {
      return p->e == q->e;
    }
  }

  // lookup a node in the unique table for the appropriate variable; insert it,
  // if it has not been found NOTE: reference counting is to be adjusted by
  // function invoking the table lookup and only normalized nodes shall be
  // stored.
  Node* lookup(Node* p) {
    // there are unique terminal nodes
    if (Node::isTerminal(p)) {
      return p;
    }

    const auto key = hash(p);
    const auto v = p->v;
    ++stats[v].lookups;

    // search bucket in table corresponding to hashed value for the given node
    // and return it if found.
    if (auto* hashedNode = searchTable(p, key); !Node::isTerminal(hashedNode)) {
      return hashedNode;
    }

    // if node not found -> add it to front of unique table bucket
    p->next = tables[v][key];
    tables[v][key] = p;
    stats[v].trackInsert();

    return p;
  }

  /**
   * @brief Increment the reference count of a node.
   * @details This is a pass-through function that calls the increment function
   * of the node. It additionally keeps track of the number of active entries
   * in the table (entries with a reference count greater than zero). Reference
   * counts saturate at the maximum value of RefCount.
   * @param p A pointer to the node to increase the reference count of.
   * @returns Whether the reference count was increased.
   * @see Node::incRef(Node*)
   */
  [[nodiscard]] bool incRef(Node* p) noexcept {
    const auto inc = ::dd::incRef(p);
    if (inc && p->ref == 1U) {
      stats[p->v].trackActiveEntry();
    }
    return inc;
  }

  /**
   * @brief Decrement the reference count of a node.
   * @details This is a pass-through function that calls the decrement function
   * of the node. It additionally keeps track of the number of active entries
   * in the table (entries with a reference count greater than zero). Reference
   * counts saturate at the maximum value of RefCount.
   * @param p A pointer to the node to decrease the reference count of.
   * @returns Whether the reference count was decreased.
   * @see Node::decRef(Node*)
   */
  [[nodiscard]] bool decRef(Node* p) noexcept {
    const auto dec = ::dd::decRef(p);
    if (dec && p->ref == 0U) {
      --stats[p->v].numActiveEntries;
    }
    return dec;
  }

  [[nodiscard]] bool possiblyNeedsCollection() const {
    return getNumEntries() >= gcLimit;
  }

  std::size_t garbageCollect(bool force = false) {
    const std::size_t numEntriesBefore = getNumEntries();
    if ((!force && numEntriesBefore < gcLimit) || numEntriesBefore == 0U) {
      return 0U;
    }

    std::size_t v = 0U;
    for (auto& table : tables) {
      auto& stat = stats[v];
      ++stat.gcRuns;
      for (auto& bucket : table) {
        Node* p = bucket;
        Node* lastp = nullptr;
        while (p != nullptr) {
          if (p->ref == 0) {
            Node* next = p->next;
            if (lastp == nullptr) {
              bucket = next;
            } else {
              lastp->next = next;
            }
            memoryManager->returnEntry(p);
            p = next;
            --stat.numEntries;
          } else {
            lastp = p;
            p = p->next;
          }
        }
      }
      stat.numActiveEntries = stat.numEntries;
      ++v;
    }

    // The garbage collection limit changes dynamically depending on the number
    // of remaining (active) nodes. If it were not changed, garbage collection
    // would run through the complete table on each successive call once the
    // number of remaining entries reaches the garbage collection limit. It is
    // increased whenever the number of remaining entries is rather close to the
    // garbage collection threshold and decreased if the number of remaining
    // entries is much lower than the current limit.
    const auto numEntries = getNumEntries();
    if (numEntries > gcLimit / 10 * 9) {
      gcLimit = numEntries + initialGCLimit;
    }
    return numEntriesBefore - numEntries;
  }

  void clear() {
    // clear unique table buckets
    for (auto& table : tables) {
      for (auto& bucket : table) {
        bucket = nullptr;
      }
    }
    gcLimit = initialGCLimit;
    for (auto& stat : stats) {
      stat.reset();
    }
  };

  void print() {
    auto q = nvars - 1U;
    for (auto it = tables.rbegin(); it != tables.rend(); ++it) {
      auto& table = *it;
      std::cout << "\tq" << q << ":"
                << "\n";
      for (std::size_t key = 0; key < table.size(); ++key) {
        auto* p = table[key];
        if (p != nullptr) {
          std::cout << "\tkey=" << key << ": ";
        }

        while (p != nullptr) {
          std::cout << "\t\t" << std::hex << reinterpret_cast<std::uintptr_t>(p)
                    << std::dec << " " << p->ref << std::hex;
          for (const auto& e : p->e) {
            std::cout << " p" << reinterpret_cast<std::uintptr_t>(e.p) << "(r"
                      << reinterpret_cast<std::uintptr_t>(e.w.r) << " i"
                      << reinterpret_cast<std::uintptr_t>(e.w.i) << ")";
          }
          std::cout << std::dec << "\n";
          p = p->next;
        }
      }
      --q;
    }
  }

private:
  /// Typedef for a bucket in the table
  using Bucket = Node*;
  /// Typedef for the table
  using Table = std::array<Bucket, NBUCKET>;

  /// The number of variables
  std::size_t nvars = 0U;
  /**
   * @brief The actual tables (one for each variable)
   * @details Each hash table is an array of buckets. Each bucket is a linked
   * list of entries. The linked list is implemented by using the next pointer
   * of the entries.
   */
  std::vector<Table> tables{nvars};

  /// A pointer to the memory manager for the nodes stored in the table.
  MemoryManager<Node>* memoryManager;

  /// A collection of statistics
  std::vector<UniqueTableStatistics> stats{nvars};

  /// The initial garbage collection limit
  std::size_t initialGCLimit;
  /// The current garbage collection limit
  std::size_t gcLimit = initialGCLimit;

  /**
  Searches for a node in the hash table with the given key.
  @param p The node to search for.
  @param key The hashed value used to search the table.
  @return The Edge<Node> found in the hash table or Edge<Node>::zero if not
  found.
  **/
  Node* searchTable(Node* p, const std::size_t& key) {
    const auto v = p->v;
    Node* bucket = tables[v][key];
    while (bucket != nullptr) {
      if (nodesAreEqual(p, bucket)) {
        // Match found
        // ! 如果可以找到相同的节点(复用节点),那么这个被分配出来的p节点就应该返还给memoryManager.
        if (p != bucket) {
          // 将节点的出边都清空
          // memset(&(p->e), 0, sizeof(p->e));
          // put node pointed to by p on available chain
          memoryManager->returnEntry(p);
        }
        ++stats[v].hits;
        return bucket;
      }
      ++stats[v].collisions;
      bucket = bucket->next;
    }

    // Node not found in bucket
    return Node::getTerminal();
  }

public:

  /**
   * @brief 更新哈希表
   * @param p 指向节点的指针
   * @param keyBefore 在更新之前p所指向节点所在的哈希桶位置
   * @note 目前这个函数是用于在补全skipped nodes之后对其父节点在哈希表上的更新
   * 也就是说目前仅支持mNode系列的哈希表,目前没有考虑vNode,dNode系列的哈希表修改
   * 不过即使是dNode或vNode对应的哈希表应该也是一样的更新方法
   * @warning 1.该函数还需要做出改进,考虑这样一种情况: 在一个DD中有两个同构的节点(其key值和v值
   * 都相同,且子节点中都包含有skipped node),那么如果我们对其中一个节点进行alterUniqueTable()处理之后
   * 另一个节点在进入alterUniqueTable时发现该其实已经被处理了(从uniqueTable)中被去除了. (已解决)
   * 2.该函数目前无法正常运行!!! (已解决)
   */
  void alterUniqueTable(Node *p, int keyBefore)
  {
    assert(keyBefore <= NBUCKET);
    // 获取节点的index值：
    const auto v = p->v;
    // TODO: 这里有误!!!

    if(tables[v][keyBefore] == nullptr) 
    {
      return;
    }
    assert(tables[v][keyBefore] != nullptr);

    Node **head = &tables[v][keyBefore];
    // 找到节点原本的桶位置
    Node** node = &tables[v][keyBefore];
    Node** pre = nullptr;

    // 需要在单链表中找到p节点指针的位置:
    while(node != nullptr && *node != nullptr)
    {
      assert(node != nullptr && *node != nullptr);
      if(*node == p)
      {
        if(*node == *head)
        {
          *head = (*head)->next;
          break;
        }
        (*pre)->next = (*node)->next;
        break;
      }
      pre = node;
      node = &((*node)->next);
    }
    // !之后只需要用lookup函数将该节点重新放回哈希表即可

    //? 因为这里仅为了移动节点在哈希表的位置,因此我个人认为统计这次移动的lookup次数没有意义???
    --stats[v].lookups;
  }


  /**
   * @brief 取出tables[index][*]中所有的哈希冲突链并返回,以便之后修改节点的哈希值
   * @note 本质上也就是取出第index层的所有节点所在的内存位置,以供之后做各种sifting变换
   */
  std::vector<Node*> getTableColumn(Qubit index)
  {
    assert(index >= 0);
    std::vector<Node*> res;
    for(auto i=0;i<NBUCKET;++i)
    {
      res.push_back(tables[index][i]);
      // 将对应列的哈希冲突链取出,并将对应列的哈希冲突链清空
      tables[index][i] = nullptr;
    }    
    // clearNextIndexTable(index-1);
    return res;
  }

  /**
   * @brief 将第index层的所有节点在哈希表的指针都清空
   */
  void clearNextIndexTable(Qubit index)
  {
    Node **p{nullptr};
    Node **pnext{nullptr};
    for(auto i=0;i<NBUCKET;++i)
    {
      p = &tables[index][i];
      while(p!=nullptr && (*p)!=nullptr)
      {
        pnext = &((*p)->next);
        *p = nullptr;
        p = pnext;
      }
    }
  }

};

} // namespace dd
