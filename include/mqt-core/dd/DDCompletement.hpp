#pragma once

#include "dd/Node.hpp"
#include "dd/Edge.hpp"
#include "dd/Package.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "ir/QuantumComputation.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <queue>

/**
 * @note 目前已经可以正常运行该补全函数,不过代码还有可改进空间,比如层序遍历的过程.
 * @date 2024/10/23
 */

// TODO: 目前为了方便起见把普通函数的声明和定义都放在一个位置,之后如果没有发现有其他问题的话再考虑将函数分离(模板函数除外)

namespace dd {

    void checkRefValue(MatrixDD root)
    {
        std::queue<MatrixDD> que;
        if(root.isTerminal())
        {
            return;
        }
        que.push(root);
        MatrixDD mdd;
        while(!que.empty())
        {
            mdd = que.front();
            que.pop();
            if(mdd.p->ref > 40000)
            {
                while(!que.empty())
                {
                    que.pop();
                }
                std::cout << "ref oveflow!\r\n";
                return;
            }
            auto edges = mdd.p->e;
            for(int i=0;i<NEDGE;++i)
            {
                if(!edges[i].isTerminal())
                {
                    que.push(edges[i]);
                }
            }
        }
        std::cout << "no ref overflow occurs! Everything ok!\r\n";
    }

    /**
     * @brief 仅作为验证测试
     */
    void checkForCorrect(MatrixDD root)
    {
        std::queue<MatrixDD> que;
        if(root.isTerminal())
        {
            return;
        }
        que.push(root);
        MatrixDD mdd;
        while(!que.empty())
        {
            mdd = que.front();
            que.pop();

            if(mdd.p->ref == 0)
            {
                // 查看ref值是否出现错误
                while(!que.empty()){
                    que.pop();
                }
                std::cout << "ref == 0 occurs in " << &mdd << "\r\n";
                return;
            }

            auto edges = mdd.p->e;
            for(int i=0;i<NEDGE;++i)
            {
                if(!edges[i].isTerminal())
                {
                    que.push(edges[i]);
                }
            }
        }
        std::cout << "ok!\r\n";
    }

    /**
     * @brief 判断是否有skipped node
     * @param mdd 输入进来的根节点
     * @return std::array<bool, NEDGE> 一个包含了NEDGE(4)的数组,其中的类型是bool,只要有一条出边指向的是
     * skipped node,其对应的数组索引位置的值即为true,否则为false
     */
    std::array<bool, NEDGE> hasSkippedSubNodes(MatrixDD &mdd)
    {
        int curIndex = mdd.p->v;
        std::array<bool, NEDGE> res{false, false, false, false};
        auto edges = mdd.p->e;
        for(size_t i=0;i<NEDGE;++i)
        {
            if(edges[i].isTerminal())
            {
                // ! 这里需要判断边的权重是1还是0,因为终端节点0和1都用(p=)nullptr来表示,以权重是0还是1(or其他复数)来表示终端节点0和1
                if(curIndex == 0 || edges[i].w.approximatelyZero())
                {
                    res[i] = false;
                } else if(curIndex > 0 && !edges[i].w.approximatelyZero()){
                    res[i] = true;
                }
                continue;
            }
            int nextIndex = edges[i].p->v;
            if(curIndex - nextIndex != 1)
            {
                // 说明出现了skipped node
                res[i] = true;
            }
        }
        return res;
    }

    /**
     * @brief 查看是否有skipped node存在
     * @param root 本质上是输入的DD树的根节点
     * @return 如果存在有skipped node则返回true,否则返回false
     * @note 目前只用该文件内的系列函数侦察mNode,即只考虑quantum matrix而不考虑quantum vector
     */
    bool isSkippedNodeExist(MatrixDD root)
    {
        std::queue<MatrixDD> que;
        if(root.isTerminal())
        {
            return false;
        }
        que.push(root);
        MatrixDD mdd;
        while(!que.empty())
        {
            mdd = que.front();
            que.pop();

            auto curIndex = mdd.p->v;
            auto edges = mdd.p->e;
            for(auto i=0;i<NEDGE;++i)
            {
                // 判断是否有skipped node:
                if(edges[i].isOneTerminal() && curIndex != 0)
                {
                    while(!que.empty())
                    {
                        que.pop();
                    }
                    return true;
                }
                if(!edges[i].isTerminal())
                {
                    auto nextIndex = edges[i].p->v;
                    if(curIndex != nextIndex + 1)
                    {
                        while (!que.empty()) {
                            que.pop();
                        }
                        return true;
                    }
                    que.push(edges[i]);
                }
            }
        }
        return false;
    }

    /**
     * @brief 打印四条出边的信息
     * @param index 父层节点的index
     * @param edges 保存四条出边的数组
     */
    void printEdgesInfo(size_t index, std::array<Edge<mNode>, NEDGE> &edges)
    {
        std::vector<std::pair<size_t, Edge<mNode>>> memEdges;
        for(size_t i=0;i<NEDGE;++i)
        {
            if(edges[i].p == nullptr)
            {
                std::cout << "  v(" << index << ")->e" << "[" << i << "] = null\r\n";
            } else {
                if(memEdges.empty())
                {
                    std::cout << "  v(" << index << ")->e" << "[" << i << "] = " 
                    << "v" << "(" << edges[i].p->v << ") with w=" << edges[i].w.toString() << "\r\n";

                    memEdges.emplace_back(i, edges[i]);
                } else {
                    // 判断出边是否都是指向同一个节点
                    for(auto & memEdge : memEdges)
                    {
                        // 判断两条边指向的节点是否是同一个(同时查看边的权重是否相同)
                        if(memEdge.second == edges[i] && memEdge.second.w.approximatelyEquals(edges[i].w))
                        {
                            std::cout << "  v(" << index << ")->e" << "[" << i << "] = "
                            << "e[" << memEdge.first << "]\r\n";
                        } else {
                            std::cout << "  v(" << index << ")->e" << "[" << i << "] = " 
                            << "v" << "(" << edges[i].p->v << ") with w=" << edges[i].w.toString() << "\r\n";
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief 打印出现了skipped node的节点及其父节点等相关信息
     */
    void printSkippedNodeInfo(MatrixDD root, int &totalSkipped)
    {
        if(root.p == nullptr)
        {
            // 如果为空树
            return;
        }
        Qubit curLevelIndex = root.p->v;
        std::array<Edge<mNode>, NEDGE> edges = root.p->e;
        for(size_t i=0;i<NEDGE;++i)
        {
            // 判断是否有出现skipped node
            if(edges[i].p == nullptr)
            {
                continue;
            }
            Qubit nextLevelIndex = edges[i].p->v;
            if(curLevelIndex-nextLevelIndex != 1)
            {
                totalSkipped += (curLevelIndex-nextLevelIndex)-1;
                std::cout << curLevelIndex << " ===> " << nextLevelIndex << ":\r\n";
                // !观察其四条出边:
                printEdgesInfo(curLevelIndex, edges);
            }
            printSkippedNodeInfo(edges[i], totalSkipped);
        }
    }

    /**
     * @deprecated 使用levelCompleteSkipped替代
     * @brief 尝试将skipped node补全,方便之后进行linear sifting变换
     * @note BUG!!! 之后可能会采用层序遍历的方式来完成这个任务
     * @date 2024/10/21
     */
    template<class Config>
    void completeSkippedNodes(MatrixDD root, Package<Config>& dd)
    {
        if(root.p == nullptr)
        {
            // 基准情况
            return;
        }
        Qubit curLevelIndex = root.p->v;
        std::array<Edge<mNode>, NEDGE> edges = root.p->e;
        for(size_t i=0;i<NEDGE;++i)
        {
            // 判断是否有出现skipped node
            if(edges[i].p == nullptr)
            {
                continue;
            }
            Qubit nextLevelIndex = edges[i].p->v;
            if(curLevelIndex-nextLevelIndex != 1)
            {
                // !说明有skipped node,需要补全:

                // 获取节点的内存管理器和哈希表
                auto nodeManager = dd.mMemoryManager;
                auto uniqueTable = dd.mUniqueTable;
                // 向内存管理器申请获得一个节点的内存
                auto p = dd->mMemoryManager.get();
                assert(p->ref == 0);        // 确保获取的节点的ref值为0
                p->v = curLevelIndex-1;
                p->flags = 0;

                // ~ 使新建的节点指向原本的skipped node:
                auto *nextNode = edges[i].p;
                edges[i].p->ref += 1;

                p->e[0].p = nextNode;
                p->e[1].p = mNode::getTerminal();
                p->e[2].p = mNode::getTerminal();
                p->e[3].p = nextNode;

                // ~ 设置p->e[0]和e->p[3]的权重为1并对该权重lookup
                p->e[0].w = dd.cn.lookup(Complex::one());
                p->e[3].w = dd.cn.lookup(Complex::one());

                // 对该新建的节点进行lookup:
                edges[i].p = uniqueTable.lookup(p);
                // 增加新建的节点的ref值
                edges[i].p->ref += 1;
            }

            // 递归处理:
            completeSkippedNodes(edges[i], dd);
        }
    }

    /************************************************************************************************** */
    /************************************************************************************************** */

    /**
     * @brief 层序遍历补全DD树
     * @param root DD的根节点
     * @param dd 管理对应DD的对象
     * @note 这里的代码还需要进行完善,尤其是对terminal节点的判断 -- 2024/10/22  √
     * @todo 修改哈希表的过程出现问题!!! -- 2024/10/23 √
     */
    template<typename Config>
    void levelCompleteSkipped(MatrixDD root, Package<Config>* dd)
    {
        std::queue<MatrixDD> nodeQue;
        lvlCmpl(root, dd, nodeQue);
    }

    /**
     * @brief 层序遍历补全DD的过程函数
     * @param root 指向根节点的边
     * @param dd 管理dd的对象
     * @param que 层序遍历要用到的队列
     * @note 补充进来了新的节点之后,其原本的父节点在哈希表中的位置应该也需要做出修改,
     * 猜测这也是之前为啥会出现莫名其妙bug的主要原因，因为没有对父节点做处理
     * 
     */
    template<typename Config>
    void lvlCmpl(MatrixDD root, Package<Config> *dd, std::queue<MatrixDD> &que)
    {
        if(root.isTerminal())
        {
            return;
        }
        que.push(root);
        MatrixDD mdd;

        while(!que.empty())
        {
            mdd = que.front();
            que.pop();

            if(mdd.isTerminal())
            {
                continue;
            }

            // 记录父节点的信息,之后要去哈希表中查找该节点并将该节点进行切换(因为其子边发生了变换,所以切换后的key值也会变换)
            std::size_t keyBefore = dd->mUniqueTable.hash(mdd.p);    // TODO: 确定要进行补全节点之后再计算key

            // 判断mdd指向的节点的子节点是否包含有skipped node:
            // ! 如果一个节点的index(v)不为0,而其出边有terminal 1,这意味着必定出现了skipped node!
            std::array<bool, NEDGE> skips = hasSkippedSubNodes(mdd);

            bool completeSkipped = false;

            for(size_t i=0;i<NEDGE;++i)
            {
                if(!skips[i])
                {
                    if(!mdd.p->e[i].isTerminal())
                    {
                        que.push(mdd.p->e[i]);
                    }
                    continue;
                }
                completeSkipped = true;
                // TODO: 需要将mdd指向的节点的四条出边的节点ref值都减1
                if(!mdd.p->e[i].isTerminal())
                {
                    if(!mdd.p->e[i].isTerminal())
                    {
                        // mdd.p->e[i].p->ref -= 1;
                        // TODO: 测试ref overflow:
                        dd->decRef(mdd.p->e[i]);
                    }
                }
                // 补全函数:
                completeSkippedNodev2(mdd, i, dd);
                que.push(mdd.p->e[i]);
            }

            if(completeSkipped)
            {
                // 子节点出现变换(即发生了补全skipped node之后才需要对父节点进行修改)
                dd->mUniqueTable.alterUniqueTable(mdd.p, keyBefore);
                mdd.p = dd->mUniqueTable.lookup(mdd.p);
            }
        }
    }

    /**
     * @brief 补全skipped nodes
     * @param root 根节点
     * @param outEdges root指向节点的哪条出边有skipped node
     * @param dd 管理节点内存,权重值,内存池的dd对象
     */
    template<typename Config>
    void completeSkippedNodev2(MatrixDD mdd, size_t outEdges, Package<Config> *dd)
    {
        if(mdd.isTerminal())
        {
            // 以防万一
            return;
        }
        Qubit supposeIndex = mdd.p->v - 1;

        // 因为要修改边指向新的节点,所以采用指针的指针
        MatrixDD *skippedEdge = &(mdd.p->e[outEdges]);

        mNode *p = dd->mMemoryManager.get();
        assert(p->ref == 0);
        p->v = supposeIndex;
        p->flags = 0;
        
        // 处理p的四条出边:
        // ? 不知道是否有误 -- 2024/10/22
        for(size_t i=0;i<NEDGE;++i)
        {
            if(i == 0 || i == 3)
            {
                // 指向原本root指向的对应节点
                p->e[i].p = skippedEdge->p;
                p->e[i].w = dd->cn.lookup(Complex::one());
                // 增加节点的ref值
                // dd->incRef(p->e[i]);
                if(!p->e[i].isTerminal())
                {
                    // p->e[i].p->ref += 1;
                    // TODO:测试ref overflow 
                    dd->incRef(p->e[i]);
                }
            } else {
                // terminal node
                p->e[i].p = nullptr;
                p->e[i].w = dd->cn.lookup(Complex::zero());
            }
        }

        skippedEdge->p = dd->mUniqueTable.lookup(p);
        // skippedEdge->p->ref += 1;
        // TODO: 测试ref overflow:
        dd->incRef(*skippedEdge);
    }


} // namespace dd