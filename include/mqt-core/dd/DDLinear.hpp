#pragma once

#include "dd/Node.hpp"
#include "dd/Edge.hpp"
#include "dd/Package.hpp"
#include "dd/FunctionalityConstruction.hpp"
#include "ir/QuantumComputation.hpp"
#include "dd/DDReorder.hpp"

#include <vector>
#include <array>

#include <cstring>
#include <cstdlib>
#include <unistd.h>

#define DEBUG_MODE 1

namespace dd {

    /**
     * @brief 选择使用哪种筛选算法的入口函数
     */
    template<typename Config>
    void reorderSelect(MatrixDD mdd, Package<Config> *dd, qc::QuantumComputation *qtc, ReorderScheme scheme, VarOrder *vo=nullptr)
    {
        switch (scheme) {
            case SCHEME_SIFTING : {
                DDOriginalSifting(mdd, dd, qtc, vo);
                break;
            }
            case SCHEME_LTRANS_LOWER : {
                DDLinearTransLower(mdd, dd, qtc, vo);
                break;
            }
            case SCHEME_LTRANS_UPPER : {
                DDLinearTransUpper(mdd, dd, qtc, vo);
                break;
            }
            case SCHEME_LTRANS_MIXED : {
                // TODO: mixed算法
                break;
            }
            case SCHEME_NONE : {
                // 无需处理
                break;
            }
        }
    }

    /**
     * @brief 根据输入的VarOrder对象对变换后的DD进行恢复操作
     * @param dd 
     * @param qtc 
     * @param vo 变量序对象指针 
     * @param pop 是否需要清理掉vo内的变换步骤
     * @note pop掉数据之后可以更加方便的对数据进行迁移
     */
    template <typename Config>
    void resetVorder(Package<Config> *dd, qc::QuantumComputation *qtc, VarOrder *vo, bool pop)
    {
        if(pop)
        {
            while(!vo->isRecordEmpty())
            {
                auto *last = vo->lastRecord();
                if(last->scheme == SCHEME_SIFTING)
                {
                    // 原始线性筛选算法
                    levelExchange(last->level, dd, qtc, last->up);
                } else {
                    // lt算法
                    linearTrans(last->level, dd, qtc, last->scheme, last->up);
                }
                vo->popRecord();
            }
        } else {
            // 利用变换步骤来恢复成原本的dd而无需pop掉
            auto i = vo->size() - 1;
            while(i >= 0)
            {
                // 获取最后一项对DD的操作,然后再次执行即可.
                auto *last = vo->at(i);
                assert(last != nullptr);
                if(last->scheme == SCHEME_SIFTING)
                {
                    levelExchange(last->level, dd, qtc, last->up);
                } else {
                    linearTrans(last->level, dd, qtc, last->scheme, last->up);
                }
                i-=1;
            }
        }
    }

    /**
     * @brief 根据输入的OptimalState对象来跟踪变量序的变换,当变量序变换到OptimalState记录的
     * 位置的时候就会退出,表明此时的dd已经被变换到最佳位置
     * @tparam Config 
     * @param dd 
     * @param qtc 
     * @param vo 
     * @param state 
     */
    template <typename Config>
    void trackVorder(Package<Config> *dd, qc::QuantumComputation *qtc, VarOrder *vo, OptimalState *state)
    {

    }

    /**
     * @brief 将筛选步骤记录到VarOrder对象中
     * @param level 指示是哪一层做变换
     * @param scheme 指示采用的是哪种变换方案
     * @param ddSize 指示采用变换方案后decision diagram的大小
     * @param up 指示是当前level层和上层还是下层节点做变换
     * @param vo 存放变换步骤的VarOrder对象指针
     */
    inline void recordStep(Qubit level, ReorderScheme scheme, size_t ddSize, bool up, VarOrder *vo)
    {
        if(vo != nullptr)
        {
            vo->record(level, scheme, ddSize, up);
        }
    }

    /**
     * @brief 撤销变换记录
     */
    inline void cancelRecord(VarOrder *vo)
    {
        if(vo!=nullptr && !vo->isRecordEmpty())
        {
            vo->popRecord();
        }
    }

    /**
     * @brief original sifting 算法, 实现第i层和第i-1层节点之间的交换
     * @param index 需要处理的哪一层节点
     * @param dd 管理decision diagram中节点和对应哈希表的dd管理器
     * @param qtc 负责解析传入的文件,之后将获取的信息传递给dd管理器以完成decision diagram的构建,
     * @param up 当该值为true时表示将index层的节点和index+1层的节点进行交换,否则与index-1层的节点交换(默认为false)
     * @return 成功返回0,失败返回-1
     * @copyright Leejxian
     */
    template<typename Config>
    void levelExchange(Qubit index, Package<Config> *dd, qc::QuantumComputation *qtc, bool up=false)
    {
        if(up)
        {
            index = index + 1;
        }
        assert(index > 0 && index < qtc->getNqubits());

        // 获取对应index的哈希冲突链
        auto table = dd->mUniqueTable.getTableColumn(index); 

        // 与下层做交换
        auto tmp = qtc->outputPermutation[index];
        qtc->outputPermutation[index] = qtc->outputPermutation[index-1];
        qtc->outputPermutation[index-1] = tmp;

        // 开始遍历哈希桶和哈希冲突链
        for(auto bucket=0;bucket<table.size();++bucket)
        {
            auto *node = table[bucket];      // 节点指针
            while(node != nullptr)
            {
                auto *next = node->next;
                if(node->ref != 0)
                {
                    // levelExchange2(node, dd);
                    lvlswap(node, dd);
                }
                node = next;
            }
        }
    }

    /**
     * @brief 实现lt算法的函数
     * @param index 被处理的层
     * @param dd 管理节点
     * @param qtc
     * @param scheme 指定使用哪种方式来     
     * @param up 指定当前index层的节点和下层节点做lt还是上层节点和当前index层的节点做lt变换
     */
    template<typename Config>
    void linearTrans(Qubit index, Package<Config> *dd, qc::QuantumComputation *qtc, ReorderScheme scheme, bool up=false)
    {
        if(scheme == SCHEME_SIFTING)
        {
            // 如果是original sifting算法,直接调用先前写好的算法函数即可
            return levelExchange(index, dd, qtc, up);
        }

        if(up)
        {
            index = index + 1;
        }
        assert(index > 0 && index < qtc->getNqubits());

        // 获取对应的index的哈希冲突链
        auto table = dd->mUniqueTable.getTableColumn(index);

        // upper和lower筛选算法不需要修改permutation

        for(auto bucket=0;bucket<table.size();++bucket)
        {
            auto *node = table[bucket];
            while(node != nullptr)
            {
                auto *next = node->next;
                if(node->ref != 0)
                {
                    if(scheme == SCHEME_LTRANS_UPPER)
                    {
                        upperlvlswap(node, dd);
                    } else if(scheme == SCHEME_LTRANS_LOWER)
                    {
                        lowerlvlswp(node, dd);
                    }
                }
                node = next;
            }
        }
    }

    /**
     * @brief 见levelExchange2函数注释
     */
    template<typename Config>
    void lvlswap(mNode *node, Package<Config> *dd)
    {
        /* 
         *   获取node指针指向的所有子节点的所有四条出边,存放到数组之中,(tips:这里可以用草稿纸演算下
         *   ,看看变量序从[x0,x1]==>[x1,x0]之后矩阵是怎么变换的,以及decision diagram的那些出边
         *   的变换规律是如何,之后再看下面这个for循环就可以明白了)
        */
        std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rearrangeEdges{};
        for(size_t i=0;i<NEDGE;++i)
        {
            auto eiw = node->e[i].w;
            for(size_t j=0;j<NEDGE;++j)
            {
                if(node->e[i].isTerminal())
                {
                    rearrangeEdges[j][i] = node->e[i];
                } else {
                    auto eijw = node->e[i].p->e[j].w;

                    rearrangeEdges[j][i] = node->e[i].p->e[j];
                    if(!eiw.exactlyOne())
                    {
                        rearrangeEdges[j][i].w = dd->cn.lookup(eiw * eijw);
                    }
                }
            }
            node->e[i].w = dd->cn.lookup(Complex::one());
        }

        for(size_t i=0;i<NEDGE;++i)
        {
            auto *newNode = dd->mMemoryManager.get();
            assert(newNode->ref == 0);
            newNode->v = node->v-1;
            newNode->flags = 0;

            for(size_t j=0;j<NEDGE;++j)
            {
                newNode->e[j] = rearrangeEdges[i][j];
            }

            auto newEdge = Edge<mNode>::normalize(newNode, newNode->e, dd->mMemoryManager, dd->cn);
            newEdge.p = dd->mUniqueTable.lookup(newEdge.p);

            if(node->e[i].isTerminal())
            {
                node->e[i] = newEdge;
            } else {
                dd->decRef(node->e[i]);
                node->e[i] = newEdge;
            }

            if(!node->e[i].isTerminal())
            {
                dd->incRef(node->e[i]);
            }
        }

        // TODO: 应该选择自己来把这个节点写入哈希表中,而不是用lookup -- 2024/10/29
        node = dd->mUniqueTable.lookup(node);
    }

    /**
     * @brief 实现upper变换的基本步骤
     */
    template <typename Config>
    void upperlvlswap(mNode *node, Package<Config> *dd)
    {
        std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rearrangeEdges{};
        size_t row = 0;
        for(size_t i=0; i<NEDGE; ++i)
        {
            auto eiw = node->e[i].w;
            for(size_t j=0;j<NEDGE;++j)
            {
                // 先判断这条应该要放在矩阵的哪个位置:
                row = j ^ i;

                if(node->e[i].isTerminal())
                {
                    rearrangeEdges[row][j] = node->e[i];
                } else {
                    auto eijw = node->e[i].p->e[j].w;

                    rearrangeEdges[row][j] = node->e[i].p->e[j];
                    if(!eiw.exactlyOne())
                    {
                        rearrangeEdges[row][j].w = dd->cn.lookup(eiw * eijw);
                    }
                }
            }
            node->e[i].w = dd->cn.lookup(Complex::one());
        }

        for(size_t i=0;i<NEDGE;++i)
        {
            auto *newNode = dd->mMemoryManager.get();
            assert(newNode->ref == 0);
            newNode->v = node->v-1;
            newNode->flags = 0;

            for(size_t j=0;j<NEDGE;++j)
            {
                newNode->e[j] = rearrangeEdges[i][j];
            }

            auto newEdge = Edge<mNode>::normalize(newNode, newNode->e, dd->mMemoryManager, dd->cn);
            newEdge.p = dd->mUniqueTable.lookup(newEdge.p);

            if(node->e[i].isTerminal())
            {
                node->e[i] = newEdge;
            } else {
                dd->decRef(node->e[i]);
                node->e[i] = newEdge;
            }

            if(!node->e[i].isTerminal())
            {
                dd->incRef(node->e[i]);
            }
        }
        node = dd->mUniqueTable.lookup(node);
    }

    /**
     * @brief 实现lower变换的基本步骤
     */
    template <typename Config>
    void lowerlvlswp(mNode *node, Package<Config> *dd)
    {
        std::array<std::array<Edge<mNode>, NEDGE>, NEDGE> rearrangeEdges{};
        size_t col = 0;
        for(size_t i=0; i<NEDGE; ++i)
        {
            auto eiw = node->e[i].w;
            for(size_t j=0;j<NEDGE;++j)
            {
                // 先判断这条应该要放在矩阵的哪个位置:
                col = j ^ i;

                if(node->e[i].isTerminal())
                {
                    rearrangeEdges[i][col] = node->e[i];
                } else {
                    auto eijw = node->e[i].p->e[j].w;

                    rearrangeEdges[i][col] = node->e[i].p->e[j];
                    if(!eiw.exactlyOne())
                    {
                        rearrangeEdges[i][col].w = dd->cn.lookup(eiw * eijw);
                    }
                }
            }
            node->e[i].w = dd->cn.lookup(Complex::one());
        }

        for(size_t i=0;i<NEDGE;++i)
        {
            auto *newNode = dd->mMemoryManager.get();
            assert(newNode->ref == 0);
            newNode->v = node->v-1;
            newNode->flags = 0;

            for(size_t j=0;j<NEDGE;++j)
            {
                newNode->e[j] = rearrangeEdges[i][j];
            }

            auto newEdge = Edge<mNode>::normalize(newNode, newNode->e, dd->mMemoryManager, dd->cn);
            newEdge.p = dd->mUniqueTable.lookup(newEdge.p);

            if(node->e[i].isTerminal())
            {
                node->e[i] = newEdge;
            } else {
                dd->decRef(node->e[i]);
                node->e[i] = newEdge;
            }

            if(!node->e[i].isTerminal())
            {
                dd->incRef(node->e[i]);
            }

        }
        node = dd->mUniqueTable.lookup(node);
    }

    /**
     * @brief 专为upper算法设计的linear transformation算法向上筛选的过程函数
     * @param mdd decision diagram的根节点边 
     * @param curLevel 当前选定的要向上筛选的起始层数level
     * @param scheme 仅能是SCHEME_LTRANS_LOWER,...UPPER和...MIXED
     * @param dd 管理节点的dd对象
     * @param qtc 
     * @param state 存储最佳位置的状态
     * @param vo 存储变换步骤的对象指针
     */
    template <typename Config>
    void linearTransUpper2Top(
        MatrixDD mdd,
        Qubit curLevel,
        Package<Config> *dd,
        qc::QuantumComputation *qtc,
        OptimalState *state,
        VarOrder *vo
    )
    {
        auto n = qtc->getNqubits() - 1;
        Qubit level = curLevel;

        while(level < n)
        {
            // step1. 向上交换permutation,之后看变换之后的dd大小:
            levelExchange(level, dd, qtc, true);
            auto osddSize = mdd.size();
            recordStep(level, SCHEME_SIFTING, osddSize, true, vo);

            // step2. 做upper变换之后记录dd大小:
            linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER, true);
            auto upddSize = mdd.size();
            recordStep(level, SCHEME_LTRANS_UPPER, upddSize, true, vo);

            // step3. 判断是哪一种方案比较好
            if(state->minddSize <= std::min(osddSize, upddSize))    
            {
                // 第一种情况,原本的dd要更小,需要取消upper变换,只记录层交换变换:
                linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER, true);
                cancelRecord(vo);
            } else if(osddSize <= std::min(state->minddSize, upddSize))
            {
                // 第二种情况,说明交换层节点之后的dd更小,撤销upper变换
                linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER, true);
                cancelRecord(vo);

                state->minddSize = osddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_SIFTING;
                state->up = true;
            } else {
                state->minddSize = upddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_LTRANS_UPPER;
                state->up = true;
            }
            level += 1;
        }
    }

    /**
     * @brief linear transformation算法向下筛选的过程函数 
     * @param mdd decision diagram的根节点边
     * @param curLevel 当前选定的要向下筛选的起始层数level
     * @param scheme 仅能是SCHEME_LTRANS_LOWER,...UPPER和...MIXED方案
     * @param dd 管理节点的dd对象
     * @param qtc 
     * @param state 记录最佳位置的状态
     * @param vo 存储变换步骤的对象指针
     * @note 目前只能作用于upper算法,之后需要设计成可以应用其他lt方案
     */
    template <typename Config>
    void linearTransUpper2Bottom(
        MatrixDD mdd,
        Qubit curLevel,
        Package<Config> *dd,
        qc::QuantumComputation *qtc,
        OptimalState *state,
        VarOrder *vo
    )
    {
        auto level = curLevel;
        while(level > 0)
        {
            // step1. 先交换层,看变换之后的dd大小
            levelExchange(level, dd, qtc);
            auto osddSize = mdd.size();
            recordStep(level, SCHEME_SIFTING, osddSize, false, vo);

            // step2. 记录使用upper算法之后的dd大小:
            linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER);
            auto upddSize = mdd.size();
            recordStep(level, SCHEME_LTRANS_UPPER, upddSize, false, vo);

            // step3. 判断是哪种方案比较好:
            if(state->minddSize <= std::min(osddSize, upddSize))
            {
                // 第一种情况,upper算法不行,需要撤销它,继续使用Original Sfiting算法移动当前层
                linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER);
                cancelRecord(vo);
            } else if(osddSize <= std::min(state->minddSize, upddSize))
            {
                // 第二情况:说明Original Sifting变换更好
                // 撤销upper
                linearTrans(level, dd, qtc, SCHEME_LTRANS_UPPER);
                cancelRecord(vo);

                state->minddSize = osddSize;
                state->scheme = SCHEME_SIFTING;
                state->optimalLevel = level;
                state->up = false;
            } else {
                state->minddSize = upddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_LTRANS_UPPER;
                state->up = false;
            }

            level -= 1;
        }
    }

    /**
     * @brief 使用lower算法逐步向上筛选
     * @param mdd 指向decision diagram的根节点的边 
     * @param curLevel 当前选定的要向上筛选的起始层数level
     * @param dd 管理节点的dd对象
     * @param qtc 
     * @param state 存储筛选过程中的最佳dd情况 
     * @param vo 存储筛选过程中的步骤
     * @date 2024/11/6
     */
    template <typename Config>
    void linearTransLower2Top(
        MatrixDD mdd,
        Qubit curLevel,
        Package<Config> *dd,
        qc::QuantumComputation *qtc,
        OptimalState *state,
        VarOrder *vo
    )
    {
        auto n = qtc->getNqubits() - 1;
        Qubit level = curLevel;

        while(level < n)
        {
            // step1. 向上交换层,之后看变换之后的dd大小
            levelExchange(level, dd, qtc, true);
            auto osddSize = mdd.size();
            recordStep(level, SCHEME_SIFTING, osddSize, true, vo);

            // step2. 向上做lower变换之后记录dd大小
            linearTrans(level, dd, qtc, SCHEME_LTRANS_LOWER, true);
            auto lwddSize = mdd.size();
            recordStep(level, SCHEME_LTRANS_LOWER, lwddSize, true, vo);

            // step3. 判断是哪一种方案比较好
            if(state->minddSize <= std::min(osddSize, lwddSize))
            {
                linearTrans(level, dd, qtc, SCHEME_LTRANS_LOWER, true);
                cancelRecord(vo);

                state->minddSize = osddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_SIFTING;
                state->up = true;
            } else {
                state->minddSize = lwddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_LTRANS_LOWER;
                state->up = true;
            }

            level += 1;
        }
    }

    //TODO: lower算法向下筛选的过程
    template <typename Config>
    void linearTransLower2Bottom(
        MatrixDD mdd,
        Qubit curLevel,
        Package<Config> *dd,
        qc::QuantumComputation *qtc,
        OptimalState *state,
        VarOrder *vo
    )
    {
        // 记录当前的起始层数
        auto level = curLevel;
        while(level > 0)
        {
            // step1. 先交换层,看变换之后的dd大小:
            levelExchange(level, dd, qtc);
            auto osddSize = mdd.size();
            recordStep(level, SCHEME_SIFTING, osddSize, false, vo);

            // step2. 记录使用lower算法之后的dd大小
            linearTrans(level, dd, qtc, SCHEME_LTRANS_LOWER);
            auto lwddSize = mdd.size();
            recordStep(level, SCHEME_LTRANS_LOWER, lwddSize, false, vo);

            // step3. 判断是哪种方案比较好
            if(state->minddSize <= std::min(osddSize, lwddSize))
            {
                // 第一种情况,lower算法效果比较差,还是原始的dd比较好,不过之后还是需要不断使用levelExchange
                // 函数不断交换层以找到最优位置,所以level交换的步骤需要保存,撤销lower算法的步骤:
                linearTrans(level, dd, qtc, SCHEME_LTRANS_LOWER);
                cancelRecord(vo);
            } else if(osddSize <= std::min(state->minddSize, lwddSize))
            {
                // 第二种情况,说明进行levelExchange之后的效果更好,lower算法在此处的效果较差,
                // 只需要记录Sifting,之后继续在适当位置调用lower算法查看能否获取更好的dd
                // 撤销lower算法:
                linearTrans(level, dd, qtc, SCHEME_LTRANS_LOWER);
                cancelRecord(vo);

                state->minddSize = osddSize;
                state->scheme = SCHEME_SIFTING;
                state->optimalLevel = level;
                state->up = false;
            } else {
                // 最后一种情况:说明我们的lower算法获得了不错的效果,保存起来
                state->minddSize = lwddSize;
                state->optimalLevel = level;
                state->scheme = SCHEME_LTRANS_LOWER;
                state->up = false;
            }

            level -= 1;
        }
    }

    /**
     * @deprecated 有bug,已经弃用 -- 2024/10/27
     * @brief original sifting 算法的实现函数
     * @param mdd 指向decision diagram的root edge
     * @param dd 
     * @param qtc 
     * @param vo 存储变换期间的步骤和dd大小
     * @note 该函数的某些循环逻辑还可以做出改进 -- 2024/10/26
     */
    template<typename Config>
    void DDOriginalSifting(MatrixDD mdd, Package<Config> *dd, qc::QuantumComputation *qtc, VarOrder *vo=nullptr)
    {
        size_t n = qtc->getNqubits()-1;
        std::vector<bool> freeLevel(n+1, true);
        Qubit level{0};

        OptimalState optimalState{};        // 记录最优位置和采用的方案
        optimalState.scheme = SCHEME_SIFTING;  // 该函数中采用的最优方案永远都是OriginalSifting

        for(size_t i=0;i<n;++i)
        {
            auto minSize = mdd.size();
            uint64_t maxActiveLevel = 0;

            for(size_t j=0;j<n;++j)
            {
                auto var = qtc->outputPermutation[j];
                if(freeLevel.at(var) && dd->active.at(var) > maxActiveLevel)
                {
                    maxActiveLevel = dd->active.at(var);
                    level = j;
                }
            }
            freeLevel.at(qtc->outputPermutation[level]) = false;

            optimalState.optimalLevel = level;

            if(level * 2 < n)
            {
                auto startPos = level;        // 记录开始的位置
                while(level > 0)
                {
                    levelExchange(level, dd, qtc);
                    auto ddSize = mdd.size();

                    recordStep(
                        level,
                        SCHEME_SIFTING,
                        ddSize,
                        false,
                        vo
                    );

                    if(ddSize < minSize)
                    {
                        minSize = ddSize;
                        optimalState.optimalLevel = level-1;
                    }
                    level -= 1;
                }

                while(level < n)
                {
                    levelExchange(level, dd, qtc, true);

                    if(level < startPos)
                    {
                        // 撤销记录
                        cancelRecord(vo);
                    } else {
                        // 记录步骤
                        auto ddSize = mdd.size();
                        recordStep(
                            level,
                            SCHEME_SIFTING,
                            ddSize,
                            true,
                            vo
                        );
                        if(ddSize < minSize)
                        {
                            minSize = ddSize;
                            optimalState.optimalLevel = level+1;
                        }
                    }

                    level += 1;
                }

                while(level > optimalState.optimalLevel)
                {
                    levelExchange(level, dd, qtc);

                    if(level > startPos)
                    {
                        // 撤销记录:
                        cancelRecord(vo);
                    } else {
                        // 记录数据
                        auto ddSize = mdd.size();
                        recordStep(
                            level,
                            SCHEME_SIFTING,
                            ddSize,
                            false,
                            vo
                        );
                    }

                    level -= 1;
                }
            } else {
                auto startPos = level;

                while(level < n)
                {
                    levelExchange(level, dd, qtc, true);
                    auto ddSize = mdd.size();

                    recordStep(
                        level,
                        SCHEME_SIFTING,
                        ddSize,
                        true,
                        vo
                    );

                    if(ddSize < minSize)
                    {
                        minSize = ddSize;
                        optimalState.optimalLevel = level+1;
                    }
                    level += 1;
                }

                while(level > 0)
                {
                    levelExchange(level, dd, qtc);

                    if(level > startPos)
                    {
                        // 撤销记录
                        cancelRecord(vo);
                    } else {
                        // 记录步骤:
                        auto ddSize = mdd.size();
                        recordStep(
                            level,
                            SCHEME_SIFTING,
                            ddSize,
                            false,
                            vo
                        );
                        if(ddSize < minSize)
                        {
                            minSize = ddSize;
                            optimalState.optimalLevel = level-1;
                        }
                    }

                    level -= 1;
                }

                while(level < optimalState.optimalLevel)
                {
                    levelExchange(level, dd, qtc, true);

                    if(level < startPos)
                    {
                        // 撤销记录
                        cancelRecord(vo);
                    } else {
                        // 记录步骤
                        auto ddSize = mdd.size();
                        recordStep(
                            level,
                            SCHEME_SIFTING,
                            ddSize,
                            true,
                            vo
                        );
                    }
                    
                    ++level;
                }
            }
        }

    }

    /**
     * @brief linear transformation -- upper算法的实现函数
     * @param mdd
     * @param dd
     * @param qtc
     * @param vo 存储变换期间的步骤和dd大小等信息
     * @date 2024/10/26 - 2024/11/3
     * @copyright leejxian
     */
    template<typename Config>
    void DDLinearTransUpper(MatrixDD mdd, Package<Config> *dd, qc::QuantumComputation *qtc, VarOrder *vo)
    {
        size_t n = qtc->getNqubits() - 1;
        std::vector<bool> freeLevel(n+1, true);

        Qubit level{0};

        for(size_t i=0;i<n;++i)
        {
            // 记录当前的decision diagram大小
            
            OptimalState optimalState{};        // 记录当前选中的层(level)最优筛选位置在哪以及所采用的是何种变换算法
            uint64_t maxActiveLevel = 0;

            VarOrder voDown(mdd, qtc);      // 记录向下筛选过程中的变换步骤
            VarOrder voUp(mdd, qtc);        // 记录向上筛选过程中的变换步骤

            // TODO: 将initialLayout改为outputpermutation之后效果反而不好,需要分析!!
            for(int j=0;j<n;++j)
            {
                auto var = qtc->outputPermutation[j];
                if(freeLevel.at(var) && dd->active.at(var) > maxActiveLevel)
                {
                    maxActiveLevel = dd->active.at(var);
                    level = j;
                }
            }
            freeLevel.at(qtc->outputPermutation[level]) = false;

            // 初始化optimalState对象 
            optimalState.minddSize = mdd.size();
            optimalState.optimalLevel = level;
            optimalState.scheme = SCHEME_NONE;

            if(level == 0)
            {
                // 刚刚好选中最底层的节点来做变换
                // 向上筛选:
                linearTransUpper2Top(mdd, level, dd, qtc, &optimalState, &voUp);
                // 找到最优的位置之后
                while(!voUp.isRecordEmpty())
                {
                    auto *last = voUp.lastRecord();
                    if(
                        last->level == optimalState.optimalLevel && 
                        last->up == optimalState.up && 
                        last->scheme == optimalState.scheme
                    )
                    {
#if DEBUG_MODE
                      if (mdd.size() != optimalState.minddSize) {
                        std::cout << "in line " << __LINE__
                                  << ", mdd.size() != "
                                     "optimalState.minddSize\r\n";
                      }
#endif
                      break;
                    }
                    if(last->scheme == SCHEME_SIFTING) {
                        levelExchange(last->level, dd, qtc, last->up);
                    } else {
                        linearTrans(last->level, dd, qtc, last->scheme, last->up);
                    }
                    voUp.popRecord();
                } 
            } else if(level == n)
            {
                // 如果刚好选中的是顶层节点:
                linearTransUpper2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                while(!voDown.isRecordEmpty())
                {
                    auto *last = voDown.lastRecord();
                    if(
                        last->level == optimalState.optimalLevel &&
                        last->up == optimalState.up && 
                        last->scheme == optimalState.scheme
                    )
                    {
#if DEBUG_MODE
                      if (mdd.size() != optimalState.minddSize) {
                        std::cout << "in line " << __LINE__
                                  << ", mdd.size() != "
                                     "optimalState.minddSize\r\n";
                      }
#endif
                      break;
                    }
                    // 筛选方案只可能是SCHEME_SIFTING或者SCHEME_LTRANS_...
                    if(last->scheme == SCHEME_SIFTING)
                    {
                        levelExchange(last->level, dd, qtc, last->up);
                    } else {
                        linearTrans(last->level, dd, qtc, last->scheme, last->up);
                    }
                    voDown.popRecord();
                }
            } else if(level * 2 < n)
            {
                auto startLevel = level;        // 记录初始位置
                // 向下筛选
                linearTransUpper2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                // 利用voDown来恢复成原本的变量序:(步骤不需要pop掉)
                resetVorder(dd, qtc, &voDown, false);
                // 向上筛选
                linearTransUpper2Top(mdd, level, dd, qtc, &optimalState, &voUp);

                // 最后需要根据optimalState的状态来恢复dd以取得最小dd
                if(optimalState.optimalLevel > startLevel)
                {
                    // 在向上筛选的过程中找到了最优方法,所以无需将整个voUp全部恢复:
                    while(!voUp.isRecordEmpty())
                    {
                        auto *last = voUp.lastRecord();
                        // 判断是否达到最好的位置
                        if(
                            last->level == optimalState.optimalLevel &&
                            last->scheme == optimalState.scheme &&
                            last->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                          if (mdd.size() != optimalState.minddSize) {
                            std::cout << "in line " << __LINE__
                                      << ", mdd.size() != "
                                         "optimalState.minddSize\r\n";
                          }
#endif

                          break;
                        }
                        if(last->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(last->level, dd, qtc, last->up);
                        } else {
                            linearTrans(last->level, dd, qtc, last->scheme, last->up);
                        }
                        voUp.popRecord();
                    }
                    // TODO 将最后的结果保存到vo指针中
                } else if(optimalState.optimalLevel == startLevel && optimalState.scheme == SCHEME_NONE)
                {
                    // 如果一开始不筛选的dd更好
                    resetVorder(dd, qtc, &voUp, true);
                } else {
                    resetVorder(dd, qtc, &voUp, true);

                    int k = 0;
                    while(k < voDown.size())
                    {
                        auto *record = voDown.at(k);
                        if(record->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(record->level, dd, qtc, record->up);
                        } else {
                            linearTrans(record->level, dd, qtc, record->scheme, record->up);
                        }
                        // 如果做完变换后发现和optimalState一致,则说明已经找到最佳位置:
                        if(
                            record->level == optimalState.optimalLevel &&
                            record->scheme == optimalState.scheme && 
                            record->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                          if (mdd.size() != optimalState.minddSize) {
                            std::cout << "in line " << __LINE__
                                      << ", mdd.size() != "
                                         "optimalState.minddSize\r\n";
                          }
#endif

                          break;
                        }
                        k += 1;
                        // TODO: 标记此处k位置,将其他的步骤数据全部pop掉,之后将结果保存在vo指针中
                    }
                }

            } else {
                // 更贴近上限,先向上筛选
                auto startLevel = level;
                // 向上筛选:
                linearTransUpper2Top(mdd, level, dd, qtc, &optimalState, &voUp); 
                // 利用voUp来恢复成原本的dd:
                resetVorder(dd, qtc, &voUp, false);
                // 向下筛选的过程:      -- 2024/11/1
                linearTransUpper2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                // 最后需要根据optimalState的记录来获取最小dd:
                if(optimalState.optimalLevel < startLevel)
                {
                    // 说明是在第二步骤向下筛选的时候找到的最优位置:
                    while(!voDown.isRecordEmpty())
                    {
                        auto *last = voDown.lastRecord();
                        if(
                            last->level == optimalState.optimalLevel &&
                            last->scheme == optimalState.scheme &&
                            last->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                            if (mdd.size() != optimalState.minddSize) {
                                std::cout << "in line " << __LINE__
                                        << ", mdd.size() != "
                                            "optimalState.minddSize\r\n";
                            }
#endif

                            break;
                        }
                        if(last->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(last->level, dd, qtc, last->up);
                        } else {
                            linearTrans(last->level, dd, qtc, last->scheme);
                        }
                        voDown.popRecord();
                    }
                } else if(optimalState.optimalLevel == startLevel && optimalState.scheme == SCHEME_NONE)
                {
                    // 只需要将所有的voDown全部恢复即可
                    resetVorder(dd, qtc, &voDown, true);
                } else {
                    resetVorder(dd, qtc, &voDown, true);
                    int k = 0;
                    while(k < voUp.size())
                    {
                        // 最优位置在向上筛选的过程中被发现
                        auto *record = voUp.at(k);
                        if(record->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(record->level, dd, qtc, record->up);
                        } else {
                            linearTrans(record->level, dd, qtc, record->scheme, record->up);
                        }
                        // 如果做完变换之后就发现和optimalState一致,那么说明已经到达最佳位置:
                        if(
                            record->level == optimalState.optimalLevel &&
                            record->scheme == optimalState.scheme &&
                            record->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                          if (mdd.size() != optimalState.minddSize) {
                            std::cout << "in line " << __LINE__
                                      << ", mdd.size() != "
                                         "optimalState.minddSize\r\n";
                          }
#endif

                          break;
                        }
                        k += 1;
                    }
                }
            }
        }
    }
    
    /**
     * @brief lower算法的实现 
     * @param mdd 
     * @param dd 
     * @param qtc 
     * @param vo 
     */
    template <typename Config>
    void DDLinearTransLower(MatrixDD mdd, Package<Config> *dd, qc::QuantumComputation *qtc, VarOrder *vo)
    {
        auto n = qtc->getNqubits() - 1;
        std::vector<bool> freeLevel(n+1, true);

        Qubit level{0};

        for(size_t i=0;i<n;++i)
        {
            OptimalState optimalState{};

            uint64_t maxActiveLevel = 0;

            VarOrder voDown(mdd, qtc);
            VarOrder voUp(mdd, qtc);

            for(int j=0;j<n;++j)
            {
                auto var = qtc->outputPermutation[j];
                if(freeLevel.at(var) && dd->active.at(var) > maxActiveLevel)
                {
                    maxActiveLevel = dd->active.at(var);
                    level = j;
                }
            }
            freeLevel.at(qtc->outputPermutation[level]) = false;

            // 初始化optimalState对象
            optimalState.optimalLevel = level;
            optimalState.scheme = SCHEME_NONE;
            optimalState.minddSize = mdd.size();

            if(level == 0)
            {
                // 刚好选中最底层来做变换,需要向上筛选.
                // 向上筛选:
                linearTransLower2Top(mdd, level, dd, qtc, &optimalState, &voUp);
                // 找到最佳位置之后:
                while(!voUp.isRecordEmpty())
                {
                    auto *last = voUp.lastRecord();
                    if(
                        last->level == optimalState.optimalLevel &&
                        last->up == optimalState.up &&
                        last->scheme == optimalState.scheme
                    )
                    {
#if DEBUG_MODE
                      if (mdd.size() != optimalState.minddSize) {
                        std::cout << "in line " << __LINE__
                                  << ", mdd.size() != "
                                     "optimalState.minddSize\r\n";
                      }
#endif
                      break;
                    }

                    // 在未达到最佳位置之前,需要一步步的回溯
                    if(last->scheme == SCHEME_SIFTING)
                    {
                        levelExchange(last->level, dd, qtc, last->up);
                    } else {
                        linearTrans(last->level, dd, qtc, last->scheme, last->up);
                    }
                    voUp.popRecord();
                }
            } else if(level == n)
            {
                // 如果选中的刚好是顶层节点:
                linearTransLower2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                while(!voDown.isRecordEmpty())
                {
                    auto *last = voDown.lastRecord();
                    if(
                        last->level == optimalState.optimalLevel &&
                        last->up == optimalState.up &&
                        last->scheme == optimalState.scheme
                    )
                    {
#if DEBUG_MODE
                      if (mdd.size() != optimalState.minddSize) {
                        std::cout << "in line " << __LINE__
                                  << ", mdd.size() != "
                                     "optimalState.minddSize\r\n";
                      }
#endif

                      break;
                    }

                    // 一步步回溯并pop掉步骤:
                    if(last->scheme == SCHEME_SIFTING)
                    {
                        levelExchange(last->level, dd, qtc, last->up);
                    } else {
                        linearTrans(last->level, dd, qtc, last->scheme, last->up);
                    }
                    voDown.popRecord();
                }
            } else if(level * 2 < n)
            {
                // 选中的层偏下,需要先向下筛选:
                auto startLevel = level;
                // 向下筛选:
                linearTransLower2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                // 接下来利用voDown来恢复成原本的dd(需要注意:此过程不能将voDown里保存的步骤pop掉)
                resetVorder(dd, qtc, &voDown, false);
                // 开始向上筛选:
                linearTransLower2Top(mdd, level, dd, qtc, &optimalState, &voUp);

                // 最后需要根据optimalState状态来恢复至最佳dd:
                if(optimalState.optimalLevel > startLevel)
                {
                    // 在向上筛选的过程中找到了最优方法,所以不用将整个voUp全部pop掉:
                    while(!voUp.isRecordEmpty())
                    {
                        auto *last = voUp.lastRecord();
                        // 判断是否达到最好的位置
                        if(
                            last->level == optimalState.optimalLevel &&
                            last->scheme == optimalState.scheme &&
                            last->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                          if (mdd.size() != optimalState.minddSize) {
                            std::cout << "in line " << __LINE__
                                      << ", mdd.size() != "
                                         "optimalState.minddSize\r\n";
                          }
#endif
                          break;
                        }
                        if(last->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(last->level, dd, qtc, last->up);
                        } else {
                            linearTrans(last->level, dd, qtc, last->scheme, last->up);
                        }
                        voUp.popRecord();
                    }
                } else if(optimalState.optimalLevel == startLevel && optimalState.scheme == SCHEME_NONE)
                {
                    // 如果一开始没经过筛选的dd更好,那么直接按照voUp恢复即可
                    resetVorder(dd, qtc, &voUp, true);

                    int k = 0;
                    while(k < voDown.size())
                    {
                        auto *record = voDown.at(k);
                        if(record->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(record->level, dd, qtc, record->up);
                        } else {
                            linearTrans(record->level, dd, qtc, record->scheme, record->up);
                        }
                        // 如果做完变换后发现和optimalState一致,则说明已经找到最佳位置:
                        if(
                            record->level == optimalState.optimalLevel &&
                            record->scheme == optimalState.scheme && 
                            record->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                            if (mdd.size() != optimalState.minddSize) {
                                std::cout << "in line " << __LINE__
                                        << ", mdd.size() != "
                                            "optimalState.minddSize\r\n";
                            }
#endif
                          break;
                        }
                        k += 1;

                        // TODO: 之后需要将保存好的步骤移动到vo对象中
                    }
                }
            } else {
                // 更贴近上限,先向上筛选
                auto startLevel = level;
                // 向上筛选:
                linearTransLower2Top(mdd, level, dd, qtc, &optimalState, &voUp);
                // 利用voUp来恢复成原本的dd:
                resetVorder(dd, qtc, &voUp, false);
                // 向下筛选:
                linearTransLower2Bottom(mdd, level, dd, qtc, &optimalState, &voDown);
                // 最后需要根据optimalState的记录来获取最小dd:
                if(optimalState.optimalLevel < startLevel)
                {
                    // 说明最优方案是在向下筛选的过程中找到的:
                    while(!voDown.isRecordEmpty())
                    {
                        auto *last = voDown.lastRecord();
                        if(
                            last->level == optimalState.optimalLevel &&
                            last->scheme == optimalState.scheme &&
                            last->up == optimalState.up
                        )
                        {
#if DEBUG_MODE
                            if (mdd.size() != optimalState.minddSize) {
                                std::cout << "in line " << __LINE__
                                        << ", mdd.size() != "
                                            "optimalState.minddSize\r\n";
                            }
#endif
                            break;
                        }
                        if(last->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(last->level, dd, qtc, last->up);
                        } else {
                            linearTrans(last->level, dd, qtc, last->scheme);
                        }
                        voDown.popRecord();
                    }
                } else if(optimalState.optimalLevel == startLevel && optimalState.scheme == SCHEME_NONE)
                {
                    // 只需要将所有的voDown全部恢复即可
                    resetVorder(dd, qtc, &voDown, true);
                } else {
                    resetVorder(dd, qtc, &voDown, true);

                    int k = 0;
                    while(k < voUp.size())
                    {
                        // 最优位置在向上筛选的过程中被发现
                        auto *record = voUp.at(k);
                        if(record->scheme == SCHEME_SIFTING)
                        {
                            levelExchange(record->level, dd, qtc, record->up);
                        } else {
                            linearTrans(record->level, dd, qtc, record->scheme, record->up);
                        }
                        // 如果做完变换之后就发现和optimalState一致,那么说明已经到达最佳位置:
                        if(
                            record->level == optimalState.optimalLevel &&
                            record->scheme == optimalState.scheme &&
                            record->up == optimalState.up
                        )
                        {                   
#if DEBUG_MODE
                            if (mdd.size() != optimalState.minddSize) {
                                std::cout << "in line " << __LINE__
                                        << ", mdd.size() != "
                                            "optimalState.minddSize\r\n";
                            }
#endif

                            break;
                        }
                        k += 1;
                    }
                    
                }
            }
        }
    }


} // namespace dd