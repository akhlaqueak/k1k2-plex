#include <cassert>
#include <functional>
#include "utils.h"
#include "timer.h"
#include "MRandSet.h"
#include "MStack.h"
#include "degeneracy.h"
#include "graph.h"
#include "graphIO.h"
#include "sequence.h"
#include <chrono>

using namespace std::chrono;

int k;
int lb; // size下界
int bd; // lb-k plex中度的下界
// thread_local long cntT=0;
thread_local size_t cntT = 0;
thread_local size_t taskCnt = 0;
int overCnt = 0;

namespace ListPlex
{
    using VtxSet = RandSet<int>;
    using Graph = graph<int>;
    enum : uint8_t
    { // 奇数表示连接
        UNLINK2LESS = 0,
        LINK2LESS = 1,
        UNLINK2EQUAL = 2,
        LINK2EQUAL = 3,
        UNLINK2MORE = 4,
        LINK2MORE = 5
    };
    struct PlexEmitor;
    struct DecomposeMaxChecker;
    template <typename MaximalityChecker, typename Emitor>
    struct KplexListor
    {
        const Graph &subg;
        MaximalityChecker *plexMaxChecker;
        Emitor *emitor;                   // 输出Plex
        int *neiInP;                      // 各个顶点在P中的度
        int *const neiInG;                // 各个顶点在P+Cand中的度
        Stack<int> *const plex;           // Plex
        VtxSet *const cand1;              // 种子点的一阶邻居
        VtxSet *const cand2;              // 种子点的二阶邻居
        VtxSet *excl;                     // 排除的顶点
        VtxSet *const exclBK;             // 临时保存excl，用于恢复
        Stack<int> *const exclStack;      // 临时保存excl，用于恢复
        const int hopSz, hop2Sz;          // 种子点+他的一阶邻居size，种子点的二阶邻居size
        const vector<uint8_t> &commonMtx; // 邻接矩阵

        int *nonadjInP;

        KplexListor(const Graph &_subg, MaximalityChecker *_plexMaxChecker, Emitor *_emitor,
                    Stack<int> *_plex, VtxSet *_cand1, VtxSet *_cand2, VtxSet *_excl, VtxSet *_exclBK, Stack<int> *_exclStack,
                    int *_neiInP, int *_neiInG, int *_nonadjInP, int _hopSz, int _hop2Sz, const vector<uint8_t> &_commonMtx)
            : subg(_subg), plexMaxChecker(_plexMaxChecker), emitor(_emitor), plex(_plex), cand1(_cand1), cand2(_cand2), excl(_excl), exclBK(_exclBK), exclStack(_exclStack), neiInP(_neiInP), neiInG(_neiInG), nonadjInP(_nonadjInP), hopSz(_hopSz), hop2Sz(_hop2Sz), commonMtx(_commonMtx)
        {
        }

        KplexListor(const KplexListor &o)
            : subg(o.subg), plexMaxChecker(o.plexMaxChecker), emitor(o.emitor), plex(new Stack<int>(*(o.plex))), cand1(new VtxSet(*(o.cand1))), cand2(nullptr), excl(new VtxSet(*(o.excl))), exclBK(nullptr), exclStack(new Stack<int>(o.exclStack->cap)), neiInP(new int[subg.n]), neiInG(new int[subg.n]), nonadjInP(new int[subg.n]), hopSz(o.hopSz), hop2Sz(o.hop2Sz), commonMtx(o.commonMtx)
        {
            memcpy(neiInP, o.neiInP, sizeof(int) * subg.n);
            memcpy(neiInG, o.neiInG, sizeof(int) * subg.n);
        }
        void del()
        {
            delete plex;
            delete cand1;
            delete excl;
            delete exclStack;
            delete[] neiInP;
            delete[] neiInG;
            delete[] nonadjInP;
        }
        // 获取邻接矩阵中的位置
        int getIdx(const int v1, const int v2)
        {
            return v1 * subg.n + v2;
        }
        // 是否相连
        bool isAdjMtx(const int v1, const int v2)
        {
            return commonMtx[getIdx(v1, v2)] & 1;
        }
        // 当顶点u从Excl加入P+Cand1中时更新其邻居在P+Cand中的度
        void addG(const int u)
        {
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInG[nei]++;
            }
        }
        // 当顶点u从P+Cand1移入Excl中时更新其邻居在P+Cand中的度
        void subG(const int u)
        {
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInG[nei]--;
            }
        }
        // 当顶点u加入P时，更新他的邻居在P中的度
        void addP(const int u)
        {
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInP[nei]++;
            }
        }
        // 当顶点u移出P时，更新他的邻居在P中的度
        void subP(const int u)
        {
            for (int i = 0; i < subg.V[u].degree; i++)
            {
                const int nei = subg.V[u].Neighbors[i];
                neiInP[nei]--;
            }
        }

        // 目前效果最好的upper bound
        bool upperbound(const int u)
        {
            for (int t = 0; t < plex->sz; t++)
            {
                const int v = plex->members[t];
                nonadjInP[t] = plex->sz - neiInP[v];
            }
            int cnt = k + neiInP[u];
            for (int i = 0; i < cand1->sz; i++)
            {
                int max_noadj = -1;
                int max_index = 0;
                const int ele = cand1->members[i];
                if (isAdjMtx(u, ele))
                {
                    for (int j = 0; j < plex->sz; j++)
                    {
                        const int v = plex->members[j];
                        if (!isAdjMtx(ele, v) && nonadjInP[j] > max_noadj)
                        {
                            max_noadj = nonadjInP[j];
                            max_index = j;
                        }
                    }
                    if (max_noadj < k)
                    {
                        cnt++;
                        nonadjInP[max_index]++;
                        if (cnt >= lb)
                            return true;
                    }
                }
            }
            if (cnt >= lb)
                return true;
            return false;
        }

        // 用于构建seed set的upper_bound
        bool upperboundK()
        {
            int cnt = plex->sz;
            for (int k = 0; k < plex->sz; k++)
            {
                const int v = plex->members[k];
                nonadjInP[k] = plex->sz - neiInP[v];
            }
            for (int i = 0; i < cand1->sz; i++)
            {
                int max_noadj = -1;
                int max_index = 0;
                const int ele = cand1->members[i];
                for (int j = 1; j < plex->sz; j++)
                {
                    const int v = plex->members[j];
                    if (!isAdjMtx(ele, v) && nonadjInP[j] > max_noadj)
                    {
                        max_noadj = nonadjInP[j];
                        max_index = j;
                    }
                }
                if (max_noadj < k)
                {
                    cnt++;
                    nonadjInP[max_index]++;
                }
                if (cnt >= lb)
                    return true;
            }
            if (cnt >= lb)
                return true;
            return false;
        }

        // 从Cand1加入到Plex中
        KplexListor *cand1ToPlex(const int v)
        {
            plex->push(v);
            cand1->remove(v);
            addP(v);
            return this;
        }
        // 从plex恢复至Cand1中
        KplexListor *plexToCand1()
        {
            assert(plex->sz > 0);
            const int u = plex->top();
            cand1->add(u);
            plex->pop();
            subP(u);
            return this;
        }
        // 从Cand2中加入Plex
        int cand2BackToPlex()
        {
            const int v = cand2->pop_back();
            plex->push(v);
            addP(v);
            addG(v);
            return v;
        }
        // 从plex恢复至Cand2
        KplexListor *plexToCand2()
        {
            assert(plex->sz > 0);
            const int u = plex->top();
            cand2->add(u);
            plex->pop();
            subP(u);
            subG(u);
            return this;
        }
        // 从Cand1加入Excl
        KplexListor *cand1ToExcl(const int v)
        {
            cand1->remove(v);
            excl->add(v);
            subG(v);
            return this;
        }
        // 从Excl恢复至Cadn1
        KplexListor *exclToCand1(const int v)
        {
            excl->remove(v);
            cand1->add(v);
            addG(v);
            return this;
        }
        // 从Cand2加入到Excl
        int cand2BackToExcl()
        {
            const int u = cand2->pop_back();
            excl->add(u);
            return u;
        }
        // 从Excl恢复至Cand2
        KplexListor *exclToCand2(const int v)
        {
            excl->remove(v);
            cand2->add(v);
            return this;
        }
        // 当点v2add加入到plex中时更新Excl
        KplexListor *updateExcl(int &recExcl, const int v2add)
        {
            recExcl = excl->sz;
            for (int i = 0; i < excl->sz;)
            {
                const int ele = excl->members[i];
                if ((subg.proper && UNLINK2MORE > commonMtx[getIdx(ele, v2add)]) || !canFormPlex(ele, 1))
                {
                    excl->remove(ele);
                    exclStack->push(ele);
                }
                else
                    ++i;
            }
            recExcl -= excl->sz;
            return this;
        }
        // 当点v2add加入到plex中时更新Excl 且不用判断是否和plex构成新的plex，用于构建seed set时
        KplexListor *updateExclK(int &recExcl, const int v2add)
        {
            recExcl = excl->sz;
            for (int i = 0; i < excl->sz;)
            {
                int ele = excl->members[i];
                if (subg.proper && UNLINK2MORE > commonMtx[getIdx(ele, v2add)])
                {
                    excl->remove(ele);
                    exclStack->push(ele);
                }
                else
                    ++i;
            }
            recExcl -= excl->sz;
            return this;
        }
        // 当点v2add加入到plex中时更新Cand1
        KplexListor *updateCand1(int &recCand1, const int v2add)
        {
            recCand1 = cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                const int ele = cand1->members[i];
                if ((subg.proper && UNLINK2EQUAL > commonMtx[getIdx(ele, v2add)]) || !canFormPlex(ele, 0))
                {
                    cand1->remove(ele);
                    subG(ele);
                }
                else
                    ++i;
            }
            recCand1 -= cand1->sz;
            return this;
        }
        // 当点v2add加入到plex中时更新Cand1，但仍可在Cand1中查询到该点
        KplexListor *updateCand1Fake(int &recCand1, const int v2add)
        {
            recCand1 = cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                int ele = cand1->members[i];
                if ((subg.proper && UNLINK2EQUAL > commonMtx[getIdx(v2add, ele)]) || !canFormPlex(ele, 0))
                {
                    cand1->fakeRemove(ele);
                    subG(ele);
                }
                else
                    ++i;
            }
            recCand1 -= cand1->sz;
            return this;
        }

        // 当v2add加入到plex中时更新Cand1，仍可在Cand1中查询到该点，且不用判断是否能和plex形成更大plex，用于构建seed set时更新Cand1
        KplexListor *updateCand1KFake(int &recCand1, const int v2add)
        {
            recCand1 = cand1->sz;
            for (int i = 0; i < cand1->sz;)
            {
                int ele = cand1->members[i];
                if (subg.proper && UNLINK2EQUAL > commonMtx[getIdx(v2add, ele)])
                {
                    cand1->fakeRemove(ele);
                    subG(ele);
                }
                else
                    ++i;
            }
            recCand1 -= cand1->sz;
            return this;
        }
        // 当点v2add加入到plex中时更新Cand2，但仍可在Cand2中查询到该点
        KplexListor *updateCand2Fake(int &recCand2, const int v2add)
        {
            recCand2 = cand2->sz;
            for (int i = 0; i < cand2->sz;)
            {
                int ele = cand2->members[i];
                if (subg.proper && UNLINK2EQUAL > commonMtx[getIdx(v2add, ele)])
                {
                    cand2->fakeRemove(ele);
                }
                else
                    ++i;
            }
            recCand2 -= cand2->sz;
            return this;
        }
        // 判断u是否能和plex构成新的plex
        bool canFormPlex(const int u, const int extra)
        {
            if (neiInP[u] + k < plex->sz + 1 || neiInG[u] + k < max(lb, plex->sz) + extra)
                return false;

            for (int i = 0; i < plex->sz; i++)
            {
                const int v = plex->members[i];
                assert(v != u);
                if (neiInP[v] + k == plex->sz && !isAdjMtx(v, u))
                { // v is saturated by v,u are not neighbors
                    return false;
                }
            }
            return true;
        }
        // 递归找到所有的seed set
        void kSearch(const int res)
        {
            int recExcl = 0, recExclTmp;
            int recCand1[K_LIMIT], recCand2[K_LIMIT];
            if (cand2->sz == 0)
            {
                if (plex->sz + cand1->sz < lb)
                {
                    return;
                }
                if (plex->sz > 1 && !upperboundK())
                {
                    return;
                }
                listByCase(steady_clock::now());
                return;
            }
            // br0
            int v2delete = cand2BackToExcl();
            kSearch(res);
            exclToCand2(v2delete);
            int br = 1;
            for (; br < res; br++)
            {
                const int v2add = cand2BackToPlex();
                updateCand1KFake(recCand1[br], v2add);
                updateCand2Fake(recCand2[br], v2add);
                updateExclK(recExclTmp, v2add);
                recExcl += recExclTmp;
                if (cand2->sz)
                {
                    v2delete = cand2BackToExcl();
                    kSearch(res - br);
                    exclToCand2(v2delete);
                }
                else
                {
                    kSearch(res - br);
                    break;
                }
            }
            if (br == res)
            { // 最多能再加入1个2阶邻居
                recCand2[br] = 0;
                const int v2add = cand2BackToPlex();
                updateCand1Fake(recCand1[br], v2add);
                if (plex->sz + cand1->sz < lb)
                {
                    goto restore;
                }
                if (plex->sz > 1 && !upperboundK())
                {
                    goto restore;
                }
                VtxSet *tmp = excl;
                excl = exclBK;
                listByCase(steady_clock::now());
                excl = tmp;
            }
        restore:
            for (int i = br; i >= 1; i--)
            {
                cand1->fakeRecoverAdd(recCand1[i], this);
                cand2->fakeRecover(recCand2[i]);
                plexToCand2();
            }
            for (int i = 0; i < recExcl; ++i)
            {
                excl->add(exclStack->top());
                exclStack->pop();
            }
        }
        // 判断P+Cand1是否是一个plex
        void checkAsMaxPlex()
        {
            memcpy(plex->members + plex->sz, cand1->members, cand1->sz * sizeof(int));
            plex->sz += cand1->sz;
            int *tmp = neiInP;
            neiInP = neiInG;
            // assume plex+cand as a plex

            bool flag = true;
            for (int i = 0; i < excl->sz; ++i)
            {
                if (canFormPlex(excl->members[i], 1))
                {
                    flag = false;
                    break;
                }
            }
            if (flag && plexMaxChecker->isMaximal(plex->members, plex->sz, neiInP))
            {
                emitor->emitPlex();
            }
            plex->sz -= cand1->sz;
            neiInP = tmp;
        }
        bool isTimeout(auto start_t)
        {
            return duration_cast<milliseconds>(steady_clock::now() - start_t).count() > TIMEOUT_THRESH;
        }
        // 分支
        void branchInCandBase(const int pivot, auto start_t)
        {
            // In the first branch, select pivot
            int recCand1 = 0, recExcl = 0;
            cand1ToPlex(pivot)->updateCand1Fake(recCand1, pivot); //->updateExcl(recExcl,pivot)->listByCase();
            if (upperbound(pivot))
            {
                updateExcl(recExcl, pivot)->listByCase(start_t);
            }
            // Recover
            cand1->fakeRecoverAdd(recCand1, this);
            plexToCand1();
            while (recExcl)
            {
                excl->add(exclStack->top());
                exclStack->pop();
                recExcl--;
            }
            // In the second branch, remove pivot
            cand1ToExcl(pivot)->listByCase(start_t);
            exclToCand1(pivot);
        }

        // 分支
        void branchInCand(const int pivot, auto start_t)
        {
            // In the first branch, select pivot

            taskCnt++;
            KplexListor *listor = new KplexListor(*this);
#pragma omp task firstprivate(listor, pivot, start_t)
            {
                int recCand1 = 0, recExcl = 0;
                listor->cand1ToPlex(pivot)->updateCand1Fake(recCand1, pivot); //->updateExcl(recExcl,pivot)->listByCase();
                if (listor->upperbound(pivot))
                {
                    listor->updateExcl(recExcl, pivot)->listByCase(start_t);
                }
                listor->del();
                delete listor;
                // overCnt++;
            }

            // Recover
            //  cand1->fakeRecoverAdd(recCand1,this);
            //  plexToCand1();
            //  while(recExcl){
            //      excl->add(exclStack->top());
            //      exclStack->pop();
            //      recExcl--;
            //  }
            // In the second branch, remove pivot

            KplexListor *listor2 = new KplexListor(*this);
            taskCnt++;
#pragma omp task firstprivate(listor2, pivot, start_t)
            {
                listor2->cand1ToExcl(pivot)->listByCase(start_t);
                listor2->del();
                delete listor2;
                // overCnt++;
            }

            // cand1ToExcl(pivot)->listByCase(start_t);
            // exclToCand1(pivot);
        }

        // 找pivot并判断分支(CIKM)
        void listByCase(auto start_t)
        {
            if (plex->sz + cand1->sz < lb)
            {
                return;
            }
            if (cand1->sz == 0)
            {
                if (excl->sz == 0 && plexMaxChecker->isMaximal(plex->members, plex->sz, neiInP))
                {
                    emitor->emitPlex();
                }
                return;
            }
            int minnei_Plex = INT_MAX;
            int pivot;
            int minnei_Cand = INT_MAX;

            // line 8 9 10
            // 在P中找P∪C中有最小度的点p
            auto maxswap_Plex = [&](int u)
            {
                if (neiInG[u] < minnei_Plex)
                {
                    minnei_Plex = neiInG[u];
                    pivot = u;
                }
                else if (neiInG[u] == minnei_Cand && neiInP[u] < neiInP[pivot])
                { // TODO: impossible condition??
                    pivot = u;
                }
            };
            plex->for_each(maxswap_Plex);
            int pivot_Plex = pivot;
            if (minnei_Plex + k < max(lb, plex->sz))
                return;
            int mintt = INT_MAX;

            constexpr int grainSize = 10;

            // line 11
            // TODO: line 17-19??
            if (minnei_Plex + k < plex->sz + cand1->sz)
            {

                // 在C中找P∪C中有最小度的且不是p的邻居的点
                minnei_Cand = INT_MAX;
                auto restart = [&](int u)
                {
                    if (!isAdjMtx(u, pivot_Plex))
                    {
                        if (neiInG[u] < minnei_Cand)
                        {
                            minnei_Cand = neiInG[u];
                            pivot = u;
                        }
                        else if (neiInG[u] == minnei_Cand && neiInP[u] < neiInP[pivot])
                        {
                            pivot = u;
                        }
                    }
                };
                cand1->for_each(restart);

                // line 20-23
                if (isTimeout(start_t) && cand1->sz > grainSize)
                    branchInCand(pivot, start_t);
                else
                    branchInCandBase(pivot, start_t);
                return;
            }
            // line 13

            // 继续在P∪C中找P∪C有最小度的点p
            int minnei = minnei_Plex;
            auto maxswap_Cand = [&](int u)
            {
                if (neiInG[u] < minnei)
                {
                    minnei = neiInG[u];
                    pivot = u;
                }
                else if (neiInG[u] == minnei && neiInP[u] < neiInP[pivot])
                {
                    pivot = u;
                }
            };
            cand1->for_each(maxswap_Cand);

            // line 14-16
            if (minnei >= plex->sz + cand1->sz - k)
            {
                checkAsMaxPlex(); // P+C is a k-plex
                return;
            }
            // line 20-23
            if (isTimeout(start_t) && cand1->sz > grainSize)
                branchInCand(pivot, start_t);
            else
                branchInCandBase(pivot, start_t);
        }
    };
}