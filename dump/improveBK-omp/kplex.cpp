#include "../common/utils.h"
#include "../common/command_line.h"
#define PuCSize (P.size() + C.size())
#define MAXSUBG 1 << 12
#include <omp.h>
#define ITERATIVE_PRUNE
#define BRANCHING
#define LOOKAHEAD
#define CTCP
typedef vector<ui> VecUI;

ui k1, k2, q;
Graph *raw_g;

thread_local VecUI looka, lookb, lookc, lookd;
thread_local ui kplexes;
thread_local RandList twoHopG;

enum CommonNeighbors
{
    PM,
    MM,
    MP,
    PP
};

enum Direction
{
    Out,
    In
};

class KPlexGraph
{
    // ui kplexes;
    vector<char> pruned;    // todo change it to bitset for memory efficiency
    vector<char> prePruned; // todo change it to bitset for memory efficiency
    VecUI inDegree;
    VecUI outDegree;
    vector<VecUI> cnPP;
    vector<VecUI> cnMP;
    vector<VecUI> cnPM;
    vector<VecUI> cnMM;
    VecUI look1, look2;
    vector<MBitSet> deletedOutEdge;
    vector<pair<ui, ui>> Qe;
    VecUI Qv;

public:
    vector<MBitSet> edgeIn, edgeOut;
    VecUI degenOrder;
    VecUI recode;
    VecUI peelSeq;
  
    KPlexGraph() : pruned(raw_g->V), peelSeq(raw_g->V),
                            inDegree(raw_g->V), outDegree(raw_g->V),
                            deletedOutEdge(raw_g->V), cnPP(raw_g->V), cnPM(raw_g->V),
                            cnMP(raw_g->V), cnMM(raw_g->V), look1(raw_g->V), look2(raw_g->V),
                            recode(raw_g->V)
    {

        for (ui i = 0; i < raw_g->V; i++)
        {
            deletedOutEdge[i] = MBitSet(raw_g->nsOut[i].size());
        }
        Qe.reserve(raw_g->E / 10);

        // reset();
    }

    void decrementCN(CommonNeighbors CN, ui u, ui vind)
    {
        if (deletedOutEdge[u].test(vind))
            return;
        if (CN == PM)
        {
            // Theorem 3

            // cnPM[u][vind]--;
            if (cnPM[u][vind]-- + k1 + k2 == q)
            {
                Qe.push_back({u, vind});
            }
        }
        else if (CN == PP)
        {
            // Theorem 4

            if (cnPP[u][vind]-- + k1 + k1 == q)
            {
                Qe.push_back({u, vind});
            }
        }
        else if (CN == MP)
        {
            // Theorem 5

            if (cnMP[u][vind]-- + k1 + k2 == q)
            {
                Qe.push_back({u, vind});
            }
        }
        else if (CN == MM)
        {
            // Theorem 6
            // cout << u << " " << vind << " " << cnMM[u].size() << endl;

            if (cnMM[u][vind]-- + k2 + k2 == q)
            {
                Qe.push_back({u, vind});
            }
        }
    }
    void applyCoreTrussPruning()
    {
        calculateCNs();
        // edges are intialized already in Qe by calculateCNs()
        prePruned = pruned;
        while (true)
        {
            if (!Qe.empty())
                trussPrune();
            else if (!Qv.empty())
            {
                removeVertex();
            }
            else
                break;
        }

        compactAdjListsWithRemovedEdges();
    }

    void deleteEdge(ui u, ui vIndu)
    {
        // the edge might already have been deleted...
        if (deletedOutEdge[u].test(vIndu))
            return;

        deletedOutEdge[u].set(vIndu);

        ui v = raw_g->nsOut[u][vIndu];
        // u -> v is an edge that we want to remove...

        if (not prePruned[u] and outDegree[u]-- + k1 == q)
        {
            Qv.push_back(u);
            prePruned[u] = 1;
        }
        if (not prePruned[v] and inDegree[v]-- + k2 == q)
        {
            Qv.push_back(v);
            prePruned[v] = 1;
        }
        if (prePruned[u] or prePruned[v])
            return;
        // updating triangles now...
        // following loop check the u->w edges, that form a triangle with v->w and w->v
        // auto outLookup = getLookup(raw_g->nsOut[u]);
        // auto inLookup = getLookup(raw_g->nsIn[u]);
        Lookup outLookup(look1, raw_g->nsOut[u]);
        Lookup inLookup(look2, raw_g->nsIn[u]);

        for (ui i = 0; i < raw_g->nsOut[v].size(); i++)
        {
            ui w = raw_g->nsOut[v][i];
            if (prePruned[w] or deletedOutEdge[v].test(i))
                continue;
            ui ind = outLookup[w];
            if (ind and not deletedOutEdge[u].test(ind - 1)) // w is found in out-neigh of u
            {
                decrementCN(MM, v, i);
                decrementCN(PM, u, ind - 1);
            }
            if (inLookup[w]) // w is found in in-neigh of u
            {
                ui uIndw = getLowerBound(raw_g->nsOut[w], u);
                if (not deletedOutEdge[w].test(uIndw))
                {
                    decrementCN(MP, v, i);
                    decrementCN(MP, w, uIndw);
                }
            }
        }
        for (ui i = 0; i < raw_g->nsIn[v].size(); i++)
        {
            ui w = raw_g->nsIn[v][i];
            if (prePruned[w])
                continue;
            ui vIndw = getLowerBound(raw_g->nsOut[w], v);
            if (deletedOutEdge[w].test(vIndw))
                continue;
            if (inLookup[w]) // w is found in in-neigh of u
            {
                ui uIndw = getLowerBound(raw_g->nsOut[w], u);
                if (not deletedOutEdge[w].test(uIndw))
                {
                    decrementCN(PM, w, vIndw);
                    decrementCN(PP, w, uIndw);
                }
            }
            ui ind = outLookup[w];
            if (ind and not deletedOutEdge[u].test(ind - 1)) // w is found in out-neigh of u
            {
                decrementCN(MM, w, vIndw);
                decrementCN(PP, u, ind - 1);
            }
        }
    }

    void trussPrune()
    {

        for (ui i = 0; i < Qe.size(); i++)
        {
            ui u = Qe[i].first;
            ui vInd = Qe[i].second;
            // if (vInd >= raw_g->nsOut[u].size())
            // {
            //     cout << u << " " << vInd << " " << raw_g->nsOut[u].size() << endl;
            // }
            deleteEdge(u, vInd);
        }

        // all edges have been processed, so clear the contents
        // cout << Qe.size() << " Edges to be removed... " << endl;
        Qe.clear();
    }

    void calculateCNs()
    {
        // create space, sine each out-edge has to store number of common neighbors it have...

        for (ui u = 0; u < raw_g->V; u++)
        {
            cnPP[u] = VecUI(raw_g->nsOut[u].size(), 0);
            cnPM[u] = VecUI(raw_g->nsOut[u].size(), 0);
            cnMP[u] = VecUI(raw_g->nsOut[u].size(), 0);
            cnMM[u] = VecUI(raw_g->nsOut[u].size(), 0);
        }
        for (ui u = 0; u < raw_g->V; u++)
        {
            if (pruned[u])
                continue;
            // VecUI outLookup = getLookup(raw_g->nsOut[u]);
            // VecUI inLookup = getLookup(raw_g->nsIn[u]);

            Lookup outLookup(look1, raw_g->nsOut[u]);
            Lookup inLookup(look2, raw_g->nsIn[u]);

            for (ui j = 0; j < raw_g->nsOut[u].size(); j++)
            {
                ui v = raw_g->nsOut[u][j];
                if (pruned[v])
                    continue;
                // u -> v there is an edge
                // intersectSize returns the number of elements found in intersection
                for (ui k = 0; k < raw_g->nsOut[v].size(); k++)
                {
                    ui w = raw_g->nsOut[v][k];
                    if (pruned[w])
                        continue;
                    if (outLookup[w])
                    {
                        cnPP[u][j]++;
                        cnMM[v][k]++;
                        cnPM[u][outLookup[w] - 1]++;
                    }
                    if (inLookup[w])
                        cnMP[u][j]++;
                }
            }
        }

        for (ui u = 0; u < raw_g->V; u++)
        {
            for (ui j = 0; j < raw_g->nsOut[u].size(); j++)
            {

                if (pruned[u] or pruned[raw_g->nsOut[u][j]])
                    continue;
                // Theorem 3
                if (cnPM[u][j] + k1 + k2 < q)
                    Qe.push_back({u, j});
                // Theorem 4
                else if (cnPP[u][j] + 2 * k1 < q)
                    Qe.push_back({u, j});
                // Theorem 5
                else if (cnMP[u][j] + k1 + k2 < q)
                    Qe.push_back({u, j});
                // Theorem 6
                else if (cnPM[u][j] + 2 * k2 < q)
                    Qe.push_back({u, j});
            }
        }
    }
    void degenerate()
    {
        ui kc1 = q - k1, kc2 = q - k2, pn = 0, tail = 0;
        vector<char> peeled = pruned;
        auto outDegreeTemp = outDegree;
        auto inDegreeTemp = inDegree;

        // check how many vertices have been pruned by k1k2pruning
        for (ui i = 0; i < raw_g->V; i++)
        {
            if (peeled[i])
                pn++;
        }
        cout << "k1, k2 pruned vertices: " << pn << endl;
        while (pn + degenOrder.size() < raw_g->V)
        {

            for (ui i = 0; i < raw_g->V; i++)
            {
                if (peeled[i])
                    continue;
                if (outDegreeTemp[i] == kc1 or inDegreeTemp[i] == kc2)
                {
                    degenOrder.push_back(i);
                    peeled[i] = 1;
                }
            }
            for (ui i = tail; i < degenOrder.size(); i++)
            {
                ui v = degenOrder[i];
                for (ui u : raw_g->nsIn[v])
                {
                    outDegreeTemp[u]--;
                    if (outDegreeTemp[u] == kc1 and not peeled[u])
                    {
                        degenOrder.push_back(u);
                        peeled[u] = 1;
                    }
                }
                for (ui u : raw_g->nsOut[v])
                {

                    inDegreeTemp[u]--;
                    if (inDegreeTemp[u] == kc2 and not peeled[u])
                    {
                        degenOrder.push_back(u);
                        peeled[u] = 1;
                    }
                }
            }
            kc1++;
            kc2++;
            tail = degenOrder.size();
        }
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            peelSeq[degenOrder[i]] = i;
        }
    }

    void k1k2CorePrune()
    {
        for (ui i = 0; i < raw_g->V; i++)
        {
            inDegree[i] = raw_g->nsIn[i].size();
            outDegree[i] = raw_g->nsOut[i].size();
            if (outDegree[i] + k1 < q or inDegree[i] + k2 < q)
            {
                pruned[i] = 1;
                Qv.push_back(i);
            }
        }

        for (ui i = 0; i < Qv.size(); i++)
        {
            ui v = Qv[i];
            for (ui u : raw_g->nsIn[v])
            {
                outDegree[u]--;
                if (outDegree[u] + k1 < q and not pruned[u])
                {
                    Qv.push_back(u);
                    pruned[u] = 1;
                }
            }
            for (ui u : raw_g->nsOut[v])
            {
                inDegree[u]--;
                if (inDegree[u] + k2 < q and not pruned[u])
                {
                    Qv.push_back(u);
                    pruned[u] = 1;
                }
            }
            // verex is is pruned, so no need to keep in/out degree vertices for it
            raw_g->nsIn[v].clear();
            raw_g->nsOut[v].clear();
        }
        Qv.clear();
        compactAdjLists();
    }

    void removeVertex()
    {

        ui u = Qv.back();
        Qv.pop_back();
        if (pruned[u])
            return;
        pruned[u] = 1;

        // VecUI outLookup = getLookup(raw_g->nsOut[u]);
        // VecUI inLookup = getLookup(raw_g->nsIn[u]);

        Lookup outLookup(look1, raw_g->nsOut[u]);
        Lookup inLookup(look2, raw_g->nsIn[u]);
        // All out-edges are deleted

        for (ui i = 0; i < raw_g->nsOut[u].size(); i++)
        {
            ui v = raw_g->nsOut[u][i];
            // u->v is an edge
            if (pruned[v] or deletedOutEdge[u].test(i))
                continue;
            // adding now
            deletedOutEdge[u].set(i);
            if (not prePruned[v] and inDegree[v]-- + k2 == q)
            {
                Qv.push_back(v);
                prePruned[v] = 1;
                continue;
            }

            for (ui j = 0; j < raw_g->nsOut[v].size(); j++)
            {
                ui w = raw_g->nsOut[v][j];
                if (pruned[w] or deletedOutEdge[v].test(j))
                    continue;
                ui ind = outLookup[w];
                if (ind and not deletedOutEdge[u].test(ind - 1))
                    decrementCN(MM, v, j);
                if (inLookup[w])
                {
                    ind = getLowerBound(raw_g->nsOut[w], u);
                    if (not deletedOutEdge[w].test(ind))
                        decrementCN(MP, v, j);
                }
            }
        }
        for (ui i = 0; i < raw_g->nsIn[u].size(); i++)
        {
            ui v = raw_g->nsIn[u][i];
            // u->v is an edge
            if (pruned[v])
                continue;
            ui uIndv = getLowerBound(raw_g->nsOut[v], u);
            if (deletedOutEdge[v].test(uIndv))
                continue;
            deletedOutEdge[v].set(uIndv);
            if (not prePruned[v] and outDegree[v]-- + k1 == q)
            {
                Qv.push_back(v);
                prePruned[v] = 1;
                continue;
            }

            for (ui j = 0; j < raw_g->nsOut[v].size(); j++)
            {
                ui w = raw_g->nsOut[v][j];
                if (pruned[w] or deletedOutEdge[v].test(j))
                    continue;
                ui ind = outLookup[w];
                if (ind and not deletedOutEdge[u].test(ind - 1))
                    decrementCN(PM, v, j);
                if (inLookup[w])
                {
                    ind = getLowerBound(raw_g->nsOut[w], u);
                    if (not deletedOutEdge[w].test(ind))
                        decrementCN(PP, v, j);
                }
            }
        }

        // Here calling erase explicitly, because in/out edge list of u is being cleared afterwards...
        inLookup.erase();
        outLookup.erase();

        raw_g->nsIn[u].clear();
        raw_g->nsOut[u].clear();
        inDegree[u] = outDegree[u] = 0;
    }

    void compact(auto &adjList, auto &deletedEdges)
    {
        ui ind = 0;
        for (ui j = 0; j < adjList.size(); j++)
        {
            ui v = adjList[j];
            if (pruned[v] or deletedEdges.test(j))
                continue;
            adjList[ind++] = v;
        }
        adjList.resize(ind);
    }

    void compact(auto &adjList)
    {
        ui ind = 0;
        for (ui j = 0; j < adjList.size(); j++)
        {
            ui v = adjList[j];
            if (not pruned[v])
                adjList[ind++] = v;
        }
        adjList.resize(ind);
    }

    void compactAdjLists()
    {
        for (ui u = 0; u < raw_g->V; u++)
        {
            if (!pruned[u])
            {
                // out nieghbors are pruned by deleted edges as well
                compact(raw_g->nsOut[u]);
                compact(raw_g->nsIn[u]);
            }
        }
    }

    void compactAdjListsWithRemovedEdges()
    {
        for (ui u = 0; u < raw_g->V; u++)
            raw_g->nsIn[u].clear();

        for (ui u = 0; u < raw_g->V; u++)
        {
            if (pruned[u])
            {
                continue;
            }
            // out nieghbors are pruned by deleted edges as well
            // compact(raw_g->nsOut[u], deletedOutEdge[u]);
            compact(raw_g->nsOut[u], deletedOutEdge[u]);

            for (ui v : raw_g->nsOut[u])
            {
                raw_g->nsIn[v].push_back(u);
            }
        }
        for (ui u = 0; u < raw_g->V; u++)
            sort(raw_g->nsIn[u].begin(), raw_g->nsIn[u].end());
        ui k = 0;
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            ui u = degenOrder[i];
            if (pruned[u])
                continue;
            degenOrder[k++] = u;
        }
        degenOrder.resize(k);
        initNeighborsMapping();
    }

    void initNeighborsMapping()
    {
        ui VV = degenOrder.size();
        cout << "Remaining vertices to process..." << VV << endl;
        edgeOut.resize(VV);
        edgeIn.resize(VV);
        for (ui i = 0; i < VV; i++)
        {
            edgeOut[i] = MBitSet(VV);
            edgeIn[i] = MBitSet(VV);
            recode[degenOrder[i]] = i;
        }
        for (ui i = 0; i < VV; i++)
        {
            ui u = degenOrder[i];
            ui ru = recode[u];
            for (ui j = 0; j < raw_g->nsOut[u].size(); j++)
                edgeOut[ru].set(recode[raw_g->nsOut[u][j]]);

            for (ui j = 0; j < raw_g->nsIn[u].size(); j++)
                edgeIn[ru].set(recode[raw_g->nsIn[u][j]]);
        }
    }

    // calculates two-hop iterative pruned graph according to Algo 2
    void getTwoHopIterativePrunedG(ui s)
    {
        auto &nsIn = raw_g->nsIn[s];
        auto &nsOut = raw_g->nsOut[s];
        addTo2HopG(s);
        // auto inLookup = getLookup(nsIn);
        // auto outLookup = getLookup(nsOut);
        Lookup inLookup(looka, nsIn);
        Lookup outLookup(lookb, nsOut);
        // both in-out neighbors
        VecUI nsInOut;
        // exclusive in/out neighbors
        VecUI I;
        VecUI O;
        nsInOut.reserve(min(nsIn.size(), nsOut.size()));
        I.reserve(nsIn.size());
        O.reserve(nsOut.size());

        // VecUI existsIn(raw_g->V, 0);
        // VecUI existsOut(raw_g->V, 0);
        Lookup existsIn(lookc, nsIn);
        Lookup existsOut(lookd, nsOut);
        ui round = 1;

        for (ui u : nsIn)
        {
            if (outLookup[u])
                nsInOut.push_back(u);
            else
            {
                I.push_back(u);
                existsIn[u] = 1;
            }
        }
        for (ui u : nsOut)
        {
            if (not inLookup[u])
            {
                O.push_back(u);
                existsOut[u] = 1;
            }
        }
        ui inSize, outSize;
        do
        {
            // cout << s << " I: " << I.size() << " O: "<< O.size() << endl;
            inSize = I.size();
            outSize = O.size();

            VecUI tempIn, tempOut;
            tempIn.reserve(inSize);
            tempOut.reserve(outSize);

            for (auto v : nsInOut)
            {
                I.push_back(v);
                O.push_back(v);
                // Now these vectors are like SO and SI
            }
            for (auto u : O)
            {
                for (auto v : raw_g->nsOut[u])
                {
                    if (existsIn[v] == round)
                    {
                        existsIn[v]++;
                        tempIn.push_back(v);
                    }
                }
            }

            for (auto u : I)
            {
                for (auto v : raw_g->nsIn[u])
                {
                    // does it exist in last iteration of O
                    if (existsOut[v] == round)
                    {
                        existsOut[v]++;
                        tempOut.push_back(v);
                    }
                }
            }

            round++;
            I.swap(tempIn);
            O.swap(tempOut);
        } while (I.size() < inSize or O.size() < outSize);

        // add Exclusive-out neighbors...
        for (auto u : O)
        {
            addTo2HopG(u);
        }
        // add Exclusive-In neighbors...
        for (auto u : I)
        {
            addTo2HopG(u);
        }

        // add bi-directional nieghbors...
        // Adding in/out nieghbors to exclusive neighs.
        for (auto u : nsInOut)
        {
            I.push_back(u);
            O.push_back(u);
            addTo2HopG(u);
        }
        // erasing previous lookups so that we can reutilize the memory
        existsIn.erase();
        existsOut.erase();
        // since in/out neighbors are added, now onwards I is SI and O is SO
        // Calculate and B to two hop graph
        VecUI temp;
        Lookup intersect(lookd, temp);
        for (auto u : O)
        {
            for (auto v : raw_g->nsOut[u])
            {
                // v is not already added in 2hop graph
                if (!in2HopG(v))
                {
                    temp.push_back(v);
                    intersect[v] = 2;
                }
            }
        }
        for (auto u : O)
        {
            for (auto v : raw_g->nsIn[u])
            {
                if (intersect[v] == 2)
                {
                    intersect[v] = 3;
                }
            }
        }
        for (auto u : I)
        {
            for (auto v : raw_g->nsOut[u])
            {
                if (intersect[v] == 3)
                    intersect[v] = 4;
            }
        }
        for (auto u : I)
        {
            for (auto v : raw_g->nsIn[u])
            {
                if (intersect[v] == 4)
                {
                    addTo2HopG(v);
                }
            }
        }

        intersect.erase();
    }
    bool in2HopG(ui u)
    {
        return twoHopG.contains(u);
    }
    void addTo2HopG(ui u)
    {

        if (pruned[u] or in2HopG(u))
            return;

        twoHopG.add(u); // thread private data member
    }
    void print(string msg, auto &vec)
    {
        cout << msg;
        for (auto u : vec)
            if (u)
                cout << u << " ";
        cout << endl;
    }
};

KPlexGraph *g;

class EnumKPlex
{
    ui vi; // current vertex in degeneracy order for which we are searching kplex
    VecUI dPin;
    VecUI dPout;
    // G is a graph induced by P \cup C
    VecUI dGin;
    VecUI dGout;
    RandList C;
    RandList X;
    RandList P;
    VecUI rC, rX;
    VecUI block;

public:
    Lookup inBlock;
    EnumKPlex(VecUI blk, ui vi) : block(blk),
                                  dPin(blk.size()), dPout(blk.size()),
                                  dGin(blk.size()), dGout(blk.size()),
                                  C(blk.size()), X(blk.size()), P(blk.size()),
                                  inBlock(looka, block)

    {
        rC.reserve(blk.size()*10);
        rX.reserve(blk.size()*10);
        for (ui i=0;i<blk.size(); i++)
        {
            ui u = blk[i];
            if (g->peelSeq[u] < g->peelSeq[vi])
                X.add(i);
            else
                addToC(i);
        }
    }

    EnumKPlex(const EnumKPlex &o) : dPin(o.dPin), dPout(o.dPout),
                                    dGin(o.dGin), dGout(o.dGout),
                                    C(o.C), X(o.X), P(o.P),
                                    rC(o.rC), rX(o.rX), 
                                    block(o.block), inBlock(looka, o.block)
    {
    }
    void recurSearch(ui u)
    {
        if (PuCSize < q)
            return;
        CToP(u);
        // any vertex removed by X or C will be appended in rX, rC so that it can be recovered later

        ui rc = updateC();
        ui rx = updateX();
        // cout<<P.size()<< ", "<<C.size()<<" "<<dGout[u]<<" "<<dGin[u]<< " : ";

        branch();

        // recover
        PToC(u);
        recoverC(rc);
        recoverX(rx);
    }

    void branch()
    {
        // cout<<"("<<P.size()<< ","<<C.size()<<")";

        if (PuCSize < q)
            return;
        if (C.empty())
        {
            if (X.empty())
                reportSolution();
            return;
        }
#ifdef BRANCHING
        VecUI MOut, MIn;
        MOut.reserve(P.size());
        MIn.reserve(P.size());
        calculateM(MOut, MIn);

        if (MOut.size() or MIn.size())
        {
            Direction dir;
            ui vp = pickvp(MOut, MIn, dir);
            {
                multiRecurSearch(vp, dir);
            }
            return;
        }
        ui vpIn, vpOut;
        ui vc = minDegreeC(vpOut, vpIn);
#else
        ui vpIn, vpOut;
        ui vc = minDegreePuC(vpOut, vpIn);
#endif

// if solution is found, it is also reported in the same function
#ifdef LOOKAHEAD
        if (lookAheadSolutionExists(vpOut, vpIn))
            return;
#endif

        // Branching in Cand...
        EnumKPlex kp = EnumKPlex(*this);
#pragma omp task firstprivate(kp, vc)
        {
            kp.recurSearch(vc);
        }

        // EnumKPlex *kp = new EnumKPlex(*this);

#pragma omp task firstprivate(kp, vc)
        {

            kp.CToX(vc);
            kp.branch();
            // recover
            kp.XToC(vc);
        }
        // other branch where P contains u
    }

    void multiRecurSearch(ui vp, Direction dir)
    {
        if (PuCSize < q)
            return;
        ui p;
        VecUI vpNN; // It stores {u1, u2, ..., ud} vertices
        vpNN.reserve(C.size());

        auto getNonNeigh = [&](auto &adj)
        {
            for (ui i = 0; i < C.size(); i++)
            {
                ui u = C[i];
                ui ru = g->recode[block[u]];
                if (!adj.test(ru))
                    vpNN.push_back(u);
            }
        };
        ui rvp = g->recode[block[vp]];
        if (dir == Out)
        {
            getNonNeigh(g->edgeOut[rvp]);
            p = k1 - (P.size() - dPout[vp]);
        }
        else
        {
            getNonNeigh(g->edgeIn[rvp]);
            p = k2 - (P.size() - dPin[vp]);
        }
        // cout<<vpNN.size()<<" "<<p<<endl;
        if (vpNN.size() <= p)
        {
            // ! this condition should never be satisfied, check this bug...
            cout << "d<=p: " << vpNN.size() << " " << p << " ";
            return;
        }

        ui rc = 0;
        ui rx = 0;
        // Branch 0
        ui u = vpNN[0];
        CToX(u);
        branch();
        XToC(u);
        // Branches 1...p
        for (ui i = 1; i < p; i++)
        {

            ui u = vpNN[i];     // u_i of algo
            ui v = vpNN[i - 1]; // u_(i-1) of algo
            if (!C.contains(v))
                continue;
            CToP(v);
            rc += updateC();
            rx += updateX();
            // cout << rC.size() << " " << rX.size() << " . ";
            if (C.contains(u))
            {
                CToX(u);
                branch();
                XToC(u);
            }
            else
            {
                X.add(u);
                branch();
                X.remove(u);
            }
        }

        // p+1th last branch.
        // so far 1...(p-2) vertices are moved from C to P
        // now move p...d vertices from C to X

        for (ui i = p; i < vpNN.size(); i++)
        {
            ui u = vpNN[i];
            if (C.contains(u))
            {
                // cout<<"*";
                removeFromC(u);
                rC.emplace_back(u);
                rc++;
            }
        }
        u = vpNN[p - 1];
        if (C.contains(u))
        {
            CToP(u);
            rc += updateC();
            rx += updateX();
            branch();
        }
        // cout<<vpNN.size()<<" "<<p<<" "<<P.size()<<" "<<C.size()<<endl;

        // recover C from P
        for (ui i = 0; i < p; i++)
        {
            ui u = vpNN[i];
            if (P.contains(u))
                PToC(u);
        }
        // recover C and X
        recoverC(rc);
        recoverX(rx);
    }
    void calculateM(VecUI &MOut, VecUI &MIn)
    {
        for (ui i = 0; i < P.size(); i++)
        {
            ui u = P[i];
            if (dGout[u] + k1 < PuCSize)
                MOut.push_back(u);
            if (dGin[u] + k2 < PuCSize)
                MIn.push_back(u);
        }
    }
    ui pickvp(VecUI &MOut, VecUI &MIn, Direction &dir)
    {
        ui vpOut = -1, vpIn = -1;
        for (ui u : MOut)
        {
            if (vpOut == -1 or dGout[u] < dGout[vpOut])
                vpOut = u;
        }
        for (ui u : MIn)
        {
            if (vpIn == -1 or dGin[u] < dGin[vpIn])
                vpIn = u;
        }
        if (MOut.size() and MIn.size())
        {
            ui outSupport = k1 + dPout[vpOut] - P.size();
            ui inSupport = k2 + dPin[vpIn] - P.size();
            if (outSupport < inSupport)
                dir = Out;
            else
                dir = In;
        }
        else if (MOut.size())
            dir = Out;
        else if (MIn.size())
            dir = In;
        else
            cout << "#"; // if this happens there must be a problem...

        if (dir == Out)
            return vpOut;
        else
            return vpIn;
    }

    ui pickvp(ui vpOut, ui vpIn, Direction &dir)
    {

        // if both of vertices doesn't support G as a kplex
        if (dGout[vpOut] + k1 < PuCSize and dGin[vpIn] + k2 < PuCSize)
        {
            ui outSupport = k1 - (P.size() - dPout[vpOut]);
            ui inSupport = k2 - (P.size() - dPin[vpIn]);
            if (outSupport < inSupport)
                dir = Out;
            else
                dir = In;
        }
        else if (dGout[vpOut] + k1 < PuCSize)
            dir = Out;
        else if (dGin[vpIn] + k2 < PuCSize)
            dir = In;
        else
            cout << "#"; // if this happens there must be a problem...

        if (dir == Out)
            return vpOut;
        else
            return vpIn;
    }

    ui pickvpFromP(ui vp, Direction &dir)
    {

        ui vpOut = -1, vpIn = -1;
        ui vpOutDegree = 0, vpInDegree = 0;

        for (ui i = 0; i < P.size(); i++)
        {
            ui u = P[i];
            if (dGout[u] + k1 < PuCSize)
            {
                if (vpOut == -1 or dGout[u] < vpOutDegree)
                {
                    vpOut = u;
                    vpOutDegree = dGout[u];
                }
            }

            if (dGin[u] + k2 < PuCSize)
            {
                if (vpIn == -1 or dGin[u] < vpInDegree)
                {
                    vpIn = u;
                    vpInDegree = dGin[u];
                }
            }
        }
        if (vpOut != -1 and vpIn != -1)
        {
            ui outSupport = k1 - (P.size() - dPout[vpOut]);
            ui inSupport = k2 - (P.size() - dPin[vpIn]);
            if (outSupport < inSupport)
            {
                dir = Out;
                vp = vpOut;
            }
            else
            {
                dir = In;
                vp = vpIn;
            }
        }
        else if (vpOut != -1)
        {
            vp = vpOut;
            dir = Out;
        }
        else if (vpIn != -1)
        {
            vp = vpIn;
            dir = In;
        }
        // else vp is unchanged...

        return vp;
    }

    ui minDegreeC(ui &vpOut, ui &vpIn)
    {
        vpOut = vpIn = C[0];
        for (ui i = 1; i < C.size(); i++)
        {
            ui u = C[i];
            if (dGin[u] < dGin[vpIn])
                vpIn = u;
            if (dGout[u] < dGout[vpOut])
                vpOut = u;
        }
        return dGin[vpIn] < dGout[vpOut] ? vpIn : vpOut;
    }

    ui minDegreePuC(ui &vpOut, ui &vpIn)
    {
        // Find min degree vertex...
        vpOut = vpIn = P[0];

        for (ui i = 1; i < P.size(); i++)
        {
            ui u = P[i];
            if (dGin[u] < dGin[vpIn])
                vpIn = u;
            if (dGout[u] < dGout[vpOut])
                vpOut = u;
        }
        ui vc = C[0];
        for (ui i = 0; i < C.size(); i++)
        {
            ui u = C[i];
            if (dGin[u] < dGin[vpIn])
                vpIn = u;
            if (dGout[u] < dGout[vpOut])
                vpOut = u;
            if (dGout[u] < dGout[vc] or dGin[u] < dGin[vc])
                vc = u;
        }
        return vc;
    }
    void minDegreeP(ui &vpOut, ui &vpIn)
    {
        // Find min degree vertex...
        vpOut = vpIn = P[0];

        for (ui i = 1; i < P.size(); i++)
        {
            ui u = P[i];
            if (dGin[u] < dGin[vpIn])
                vpIn = u;
            if (dGout[u] < dGout[vpOut])
                vpOut = u;
        }
    }
    bool lookAheadSolutionExists(ui vpOut, ui vpIn)
    {

        if (dGout[vpOut] + k1 < PuCSize or
            dGin[vpIn] + k2 < PuCSize)
            return false;

        // get a backup of C for reovery
        VecUI ctemp = C.getData();

        for (auto u : ctemp)
            CToP(u);

        // update X and see if it's empty...
        ui rx;
        rx = updateX();
        if (X.empty())
            reportSolution();

        recoverX(rx);
        for (ui u : ctemp)
            PToC(u);

        return true;
    }

    void reportSolution()
    {
        // print("KPlex: ", P);
        kplexes++;
        return;
    }

    ui updateC()
    {
        auto it = rC.end();
        for (ui i = 0; i < C.size(); i++)
        {
            ui u = C[i];
            if (!canMoveToP(u))
            {
                rC.emplace_back(u);
            }
        }

        ui sz = distance(it, rC.end());

        for (; it != rC.end(); it++)
            removeFromC(*it);

        return sz;
    }

    ui updateX()
    {
        auto it = rX.end();
        for (ui i = 0; i < X.size(); i++)
        {
            ui u = X[i];
            if (!canMoveToP(u))
            {
                rX.emplace_back(u);
            }
        }

        ui sz = distance(it, rX.end());
        for (; it != rX.end(); it++)
            X.remove(*it);
        return sz;
    }

    void recoverX(ui sz)
    {
        for (ui i = 0; i < sz; i++)
        {
            X.add(rX.back());
            rX.pop_back();
        }
    }
    void recoverC(ui sz)
    {
        for (ui i = 0; i < sz; i++)
        {
            addToC(rC.back());
            rC.pop_back();
        }
    }


    void print(string msg, auto &vec)
    {
        cout << msg;
        for (auto u : vec)
            if (u)
                cout << u << " ";
        cout << endl;
    }

    void addToC(ui u)
    {
        C.add(u);

        for (ui v : raw_g->nsOut[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dGin[rv-1]++;
        }
        for (ui v : raw_g->nsIn[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dGout[rv-1]++;
        }
    }

    void removeFromC(ui u)
    {
        C.remove(u);
        for (ui v : raw_g->nsOut[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dGin[rv-1]--;
        }
        for (ui v : raw_g->nsIn[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dGout[rv-1]--;
        }
    }

    void PToC(ui u)
    {
        P.remove(u);
        C.add(u);
        for (ui v : raw_g->nsOut[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dPin[rv-1]--;
        }
        for (ui v : raw_g->nsIn[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dPout[rv-1]--;
        }
    }

    void CToP(const ui &u)
    {
        // assert(Cand.contains(u));
        C.remove(u);
        P.add(u);
        for (ui v : raw_g->nsOut[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dPin[rv-1]++;
        }
        for (ui v : raw_g->nsIn[block[u]]){
            ui rv = inBlock[v];
            if(rv)
            dPout[rv-1]++;
        }
    }

    void CToX(const ui &u)
    {
        removeFromC(u);
        X.add(u);
    }
    void XToC(const ui &u)
    {
        X.remove(u);
        addToC(u);
    }

    bool canMoveToP(ui u)
    {
        // u is not yet in P, hence checking <=
        if (dPout[u] + k1 <= P.size() or dPin[u] + k2 <= P.size())
            return false;
        for (ui i = 0; i < P.size(); i++)
        {
            ui v = P.get(i);

            ui ru = g->recode[block[u]];
            ui rv = g->recode[block[v]];
            // should be in-connected to every out-boundary vertex
            if (dPout[v] + k1 == P.size() && !g->edgeIn[ru].test(rv))
                return false;
            // should be out-connected to every in-boundary vertex
            if (dPin[v] + k2 == P.size() && !g->edgeOut[ru].test(rv))
                return false;
        }
        return true;
    }
};

// thread_local RandList EnumKPlex::dp

int main(int argc, char *argv[])
{
    CommandLine cmd(argc, argv);
    std::string data_file = cmd.GetOptionValue("-g", "");
    q = cmd.GetOptionIntValue("-q", 2);
    k1 = cmd.GetOptionIntValue("-k1", 1);
    k2 = cmd.GetOptionIntValue("-k2", 1);

    if (data_file == "")
    {
        cout << "Please provide data file" << endl;
        exit(-1);
    }

    if (k1 < 2 or k2 < 2)
    {
        cout << "k1, k2 should be at least 2" << endl;
        return 0;
    }

    if (q + 1 < 2 * k1 or q + 1 < 2 * k2)
    {
        cout << "q should be at least 2*k1-1, 2*k2-1" << endl;
        return 0;
    }

    cout << "Loading Started" << endl;
    raw_g= new Graph(data_file);

    cout << "n=" << raw_g->V << endl;
    cout << "Loading Done" << endl;
    g = new KPlexGraph();

    // remove vertices v with outdeg(v)+k1 < q OR indeg(v)+k2 < q
    g->k1k2CorePrune();
    // find degeneracy order, the result is degenOrder vector
    g->degenerate();
    auto tick = chrono::steady_clock::now();
    g->applyCoreTrussPruning();

    cout << " CTCP time: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - tick).count() << " ms" << endl;
    cout << "starting kplex" << endl;
    tick = chrono::steady_clock::now();
    omp_set_num_threads(20);
#pragma omp parallel
    {
        // cout<<"N: "<<omp_get_num_threads()<<endl;
        // cout<<"id: "<<omp_get_thread_num()<<endl;
        looka.resize(raw_g->V);
        lookb.resize(raw_g->V);
        lookc.resize(raw_g->V);
        lookd.resize(raw_g->V);
        twoHopG.init(raw_g->V);
#pragma omp for schedule(dynamic)
        for (ui i = 0; i < g->degenOrder.size(); i++)
        {
            ui vi = g->degenOrder[i];
            twoHopG.clear();
            g->getTwoHopIterativePrunedG(vi);
            EnumKPlex kp =  EnumKPlex(twoHopG.getData(), vi);
            vi = kp.inBlock[vi]-1;
#pragma omp taskgroup
            {
                kp.recurSearch(vi);
            }
        }
    }
    ui total = 0;
#pragma omp parallel reduction(+ : total)
    {
        total += kplexes;
    }

    cout << "Total (" << k1 << "," << k2 << ")-plexes of at least " << q << " size: " << total << endl;


    cout << data_file << " Enumeration time: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - tick).count() << " ms" << endl;

    cout << endl;
    return 0;
}