
#include "graph.h"
// #define RAPIDS
#include "command_line.h"
#include "../common/utils.h"
enum CommonNeighbors
{
    PM,
    MM,
    MP,
    PP
};
class EnumKPlex
{
    Graph &g;
    ui kplexes;
    ui k1, k2, q;
    vector<ui> vBoundaryIn;
    vector<ui> vBoundaryOut;
    vector<ui> dPin;
    vector<ui> dPout;

    // G is a graph induced by P \cup C
    vector<ui> dGin;
    vector<ui> dGout;

    vector<char> pruned;    // todo change it to bitset for memory efficiency
    vector<char> prePruned; // todo change it to bitset for memory efficiency
    vector<char> in2HopG;
    vector<ui> degenOrder;
    vector<ui> peelSeq;
    vector<ui> inDegree;
    vector<ui> outDegree;
    vector<vector<ui>> cnPP;
    vector<vector<ui>> cnMP;
    vector<vector<ui>> cnPM;
    vector<vector<ui>> cnMM;

    vector<MBitSet> deletedOutEdge;
    vector<pair<ui, ui>> Qe;
    vector<ui> Qv;
    vector<ui> counts;

    vector<ui> look1, look2, look3, look4;

    ui vi; // current vertex in degeneracy order for which we are searching kplex

    RandList C;
    RandList X;
    RandList P;

public:
    void enumerate()
    {
        // remove vertices v with outdeg(v)+k1 < q OR indeg(v)+k2 < q
        k1k2CorePrune();
        // find degeneracy order, the result is degenOrder vector
        degenerate();
        auto tick = chrono::steady_clock::now();
        applyCoreTrussPruning();
        cout << " CTCP time: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - tick).count() << " ms" << endl;
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            // getTwoHopG selects vertices in 2 hop neigbhors of i, and appends them in X, C
            // B = two hop neighbors of i
            // C = vertices u in B such that u<i
            // X = vertices u in B such that u>i
            vi = degenOrder[i];
            // getTwoHopG(vi);
            addToC(vi);
            getTwoHopIterativePrunedG(vi);
            // cout << u << " ... " << endl;

            recurSearch(vi);

            // print("aC: ", C);
            // print("aX: ", X);
            reset();
        }
        cout << "Total (" << k1 << ", " << k2 << ")-plexes of at least " << q << " size: " << kplexes << endl;
        for (ui i = 0; i < counts.size(); i++)
            if (counts[i])
                cout << "kplexes of size: " << i + 1 << " = " << counts[i] << endl;
    }

    void recurSearch(ui u)
    {
        if (C.size() + P.size() < q)
            return;

        CToP(u);
        // any vertex removed by X or C will be appended in rX, rC so that it can be recovered later
        vector<ui> rC = updateC();
        vector<ui> rX = updateX();

        branch();

        // recover
        PToC(u);
        for (auto u : rC)
            addToC(u);
        for (auto u : rX)
            X.add(u);
    }

    void branch()
    {

        if (C.size() + P.size() < q)
            return;
        if (C.empty())
        {
            if (X.empty())
                reportSolution();
            return;
        }

        // ui minu = getMinDegreeVertex();
        // if (X.empty() and lookAheadSolutionExists(minu))
        // {
        //     // P U C is the solution, so add all C to P
        //     vector<ui> Ctemp=C;
        //     for (ui u : C)
        //         CToP(u);
        //     reportSolution();
        //     // recover...
        //     for(ui u:Ctemp)
        //         PToC(u);
        //     return;
        // }
        // todo find min degree vertex in P U C
        // todo lookahead solution check
        // todo branching...
        ui u = C.get(0);
        // one branch that contains u
        recurSearch(u);
        // other branch that doesn't contain u
        CToX(u);
        branch();
        // recover
        XToC(u);
    }

    void reportSolution()
    {
        counts[P.size() - 1]++;
        // print("KPlex: ", P);
        kplexes++;
        return;
    }

    EnumKPlex(Graph &_g, ui _k1, ui _k2, ui _q) : pruned(_g.V), peelSeq(_g.V),
                                                  g(_g), inDegree(_g.V), outDegree(_g.V),
                                                  in2HopG(_g.V), dPin(_g.V), dPout(_g.V),
                                                  dGin(_g.V), dGout(_g.V),
                                                  k1(_k1), k2(_k2), q(_q), kplexes(0),
                                                  deletedOutEdge(_g.V), cnPP(_g.V), cnPM(_g.V),
                                                  cnMP(_g.V), cnMM(_g.V), look1(_g.V), look2(_g.V),
                                                  look3(_g.V), look4(_g.V),
                                                  P(_g.V), C(_g.V), X(_g.V),
                                                  counts(1000)
    {
        vBoundaryIn.reserve(_g.V);
        vBoundaryOut.reserve(_g.V);

        degenOrder.reserve(_g.V);

        for (ui i = 0; i < g.V; i++)
        {
            deletedOutEdge[i] = MBitSet(g.nsOut[i].size());
        }
        Qe.reserve(g.E / 10);
        reset();
    }

private:
    inline ui getLowerBound(auto &vec, ui x)
    {
        auto it = lower_bound(vec.begin(), vec.end(), x);

        if (it == vec.end() || *it != x)
            return vec.size();
        return it - vec.begin();
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

        compactAdjListsWithRemomvedEdges();
    }
    void printCN()
    {
        for (ui i = 0; i < g.V; i++)
        {
            cout << "************** " << i << " **************" << endl;
            print("adj: ", g.nsOut[i]);
            print("PP ", cnPP[i]);
            print("PM ", cnPM[i]);
            print("MP ", cnMP[i]);
            print("MM ", cnMM[i]);
        }
    }
    void deleteEdge(ui u, ui vIndu)
    {
        // the edge might already have been deleted...
        if (deletedOutEdge[u].test(vIndu))
            return;

        deletedOutEdge[u].set(vIndu);

        ui v = g.nsOut[u][vIndu];
        // u -> v is an edge that we want to remove...
        if (prePruned[u] or prePruned[v])
            return;

        if (outDegree[u]-- + k1 == q)
        {
            Qv.push_back(u);
            prePruned[u] = 1;
        }
        if (inDegree[v]-- + k2 == q)
        {
            Qv.push_back(v);
            prePruned[v] = 1;
        }
        // updating triangles now...
        // following loop check the u->w edges, that form a triangle with v->w and w->v
        // auto outLookup = getLookup(g.nsOut[u]);
        // auto inLookup = getLookup(g.nsIn[u]);
        Lookup outLookup(look1, g.nsOut[u]);
        Lookup inLookup(look2, g.nsIn[u]);

        for (ui i = 0; i < g.nsOut[v].size(); i++)
        {
            ui w = g.nsOut[v][i];
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
                ui uIndw = getLowerBound(g.nsOut[w], u);
                if (not deletedOutEdge[w].test(uIndw))
                {
                    decrementCN(MP, v, i);
                    decrementCN(MP, w, uIndw);
                }
            }
        }
        for (ui i = 0; i < g.nsIn[v].size(); i++)
        {
            ui w = g.nsIn[v][i];
            if (prePruned[w])
                continue;
            ui vIndw = getLowerBound(g.nsOut[w], v);
            if (deletedOutEdge[w].test(vIndw))
                continue;
            if (inLookup[w]) // w is found in in-neigh of u
            {
                ui uIndw = getLowerBound(g.nsOut[w], u);
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
            if (vInd >= g.nsOut[u].size())
            {
                cout << u << " " << vInd << " " << g.nsOut[u].size() << endl;
            }
            deleteEdge(u, vInd);
        }

        // all edges have been processed, so clear the contents
        // cout << Qe.size() << " Edges to be removed... " << endl;
        Qe.clear();
    }

    ui intersectSize(auto &X, auto &Y)
    {
        ui n = 0;
        if (X.size() < Y.size())
            for (auto x : X)
            {
                if (pruned[x])
                    continue;
                if (binary_search(Y.begin(), Y.end(), x))
                    n++;
            }
        else
            for (auto y : Y)
            {
                if (pruned[y])
                    continue;
                if (binary_search(X.begin(), X.end(), y))
                    n++;
            }
        return n;
    }
    void calculateCNs()
    {
        // create space, sine each out-edge has to store number of common neighbors it have...

        for (ui u = 0; u < g.V; u++)
        {
            cnPP[u] = vector<ui>(g.nsOut[u].size(), 0);
            cnPM[u] = vector<ui>(g.nsOut[u].size(), 0);
            cnMP[u] = vector<ui>(g.nsOut[u].size(), 0);
            cnMM[u] = vector<ui>(g.nsOut[u].size(), 0);
        }
        for (ui u = 0; u < g.V; u++)
        {
            if (pruned[u])
                continue;
            // vector<ui> outLookup = getLookup(g.nsOut[u]);
            // vector<ui> inLookup = getLookup(g.nsIn[u]);

            Lookup outLookup(look1, g.nsOut[u]);
            Lookup inLookup(look2, g.nsIn[u]);

            for (ui j = 0; j < g.nsOut[u].size(); j++)
            {
                ui v = g.nsOut[u][j];
                if (pruned[v])
                    continue;
                // u -> v there is an edge
                // intersectSize returns the number of elements found in intersection
                for (ui k = 0; k < g.nsOut[v].size(); k++)
                {
                    ui w = g.nsOut[v][k];
                    if (pruned[w])
                        continue;
                    if (outLookup[w])
                    {
                        cnPP[u][j]++;
                        cnMM[v][k]++;
                        cnPM[u][outLookup[w] - 1]++;
                    }
                    if (inLookup[w])
                    {
                        cnMP[u][j]++;
                        // cnMP[v][k]++;
                        // ui uIndw = getLowerBound(g.nsOut[w], u);
                        // cnMP[w][uIndw]++;
                    }
                }
            }
            // for (ui j = 0; j < g.nsIn[u].size(); j++)
            // {
            //     ui v = g.nsIn[u][j];
            //     if (pruned[v])
            //         continue;
            //     ui uIndv = getLowerBound(g.nsOut[v], u);
            //     // v->u there is an edge
            //     // intersectSize returns the number of elements found in intersection
            //     for (ui k = 0; k < g.nsOut[v].size(); k++)
            //     {
            //         ui w = g.nsOut[v][k];
            //         if (outLookup[w])
            //         {
            //             cnPP[v][uIndv]++;
            //             cnPM[v][k]++;
            //             cnMM[u][outLookup[w] - 1]++;
            //         }
            //         if (inLookup[w])
            //         {
            //             cnPM[v][uIndv]++;
            //             cnPP[v][k]++;
            //             ui uIndw = getLowerBound(g.nsOut[w], u);
            //             cnMM[w][uIndw]++;
            //         }
            //     }
            // }
        }

        for (ui u = 0; u < g.V; u++)
        {
            for (ui j = 0; j < g.nsOut[u].size(); j++)
            {

                if (pruned[u] or pruned[g.nsOut[u][j]])
                    continue;
                // Theorem 3
                if (cnPM[u][j] + k1 + k2 < q)
                {
                    Qe.push_back({u, j});
                }
                // Theorem 4
                else if (cnPP[u][j] + 2 * k1 < q)
                {
                    Qe.push_back({u, j});
                }
                // Theorem 5
                else if (cnMP[u][j] + k1 + k2 < q)
                {
                    Qe.push_back({u, j});
                }
                // Theorem 6
                else if (cnPM[u][j] + 2 * k2 < q)
                {
                    Qe.push_back({u, j});
                }
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
        for (ui i = 0; i < g.V; i++)
        {
            if (peeled[i])
                pn++;
        }
        cout << "k1, k2 pruned vertices: " << pn << endl;
        while (pn + degenOrder.size() < g.V)
        {

            for (ui i = 0; i < g.V; i++)
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
                for (ui u : g.nsIn[v])
                {
                    outDegreeTemp[u]--;
                    if (outDegreeTemp[u] == kc1 and not peeled[u])
                    {
                        degenOrder.push_back(u);
                        peeled[u] = 1;
                    }
                }
                for (ui u : g.nsOut[v])
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

    void getNeighbors(auto &neigh)
    {
        for (ui i = 0; i < neigh.size(); i++)
        {
            ui v = neigh[i];
            if (pruned[v] or in2HopG[v])
                continue;
            // vi is the vertex for which we are building a kplex in current iteration
            if (peelSeq[v] < peelSeq[vi])
            {
                X.add(v);
            }
            else
                addToC(v);
            in2HopG[v] = 1;
        }
    }

    void getFirstHop(ui u)
    {
        getNeighbors(g.nsOut[u]);
        getNeighbors(g.nsIn[u]);
    }

    void getSecondHop(ui c1h, ui x1h)
    {
        for (size_t i = 0; i < c1h; i++)
        {
            getFirstHop(C[i]);
        }
        for (size_t i = 0; i < x1h; i++)
        {
            getFirstHop(X[i]);
        }
    }

    void getTwoHopG(ui u)
    {
        // reset in2HopG to 0

        in2HopG[u] = 1;
        // get 1st hop neighbors, they are appended to C and X
        getFirstHop(u);

        // get 2nd hop neighbors, first hop neighbors are found until c1h, x1h
        // 2nd hop neighbors are appended afterwards
        getSecondHop(C.size(), X.size());
        // getSecondHop(c1h, x1h);
    }

    void reset()
    {
        ui sz = C.size();
        for (ui i = 0; i < sz; i++)
        {
            removeFromC(C[0]);
        }
    }

    void k1k2CorePrune()
    {
        for (ui i = 0; i < g.V; i++)
        {
            inDegree[i] = g.nsIn[i].size();
            outDegree[i] = g.nsOut[i].size();
            if (outDegree[i] + k1 < q or inDegree[i] + k2 < q)
            {
                pruned[i] = 1;
                Qv.push_back(i);
            }
        }

        for (ui i = 0; i < Qv.size(); i++)
        {
            ui v = Qv[i];
            for (ui u : g.nsIn[v])
            {
                outDegree[u]--;
                if (outDegree[u] + k1 + 1 == q and not pruned[u])
                {
                    Qv.push_back(u);
                    pruned[u] = 1;
                }
            }
            for (ui u : g.nsOut[v])
            {
                inDegree[u]--;
                if (inDegree[u] + k2 + 1 == q and not pruned[u])
                {
                    Qv.push_back(u);
                    pruned[u] = 1;
                }
            }
            // verex is is pruned, so no need to keep in/out degree vertices for it
            g.nsIn[v].clear();
            g.nsOut[v].clear();
        }
        Qv.clear();
        compactAdjLists();
    }

    vector<ui> getLookup(auto &adjList)
    {
        // lookups contain index+1 if a vertex found, else zero
        vector<ui> lookup(g.V, 0);
        for (ui i = 0; i < adjList.size(); i++)
        {
            lookup[adjList[i]] = i + 1;
        }
        return lookup;
    }
    void removeVertex()
    {

        ui u = Qv.back();
        Qv.pop_back();
        if (pruned[u])
            return;
        pruned[u] = 1;

        // vector<ui> outLookup = getLookup(g.nsOut[u]);
        // vector<ui> inLookup = getLookup(g.nsIn[u]);

        Lookup outLookup(look1, g.nsOut[u]);
        Lookup inLookup(look2, g.nsIn[u]);
        // All out-edges are deleted

        for (ui i = 0; i < g.nsOut[u].size(); i++)
        {
            ui v = g.nsOut[u][i];
            // u->v is an edge
            if (pruned[v] or deletedOutEdge[u].test(i))
                continue;
            if (not prePruned[v] and inDegree[v]-- + k2 == q)
            {
                Qv.push_back(v);
                prePruned[v] = 1;
                continue;
            }

            for (ui j = 0; j < g.nsOut[v].size(); j++)
            {
                ui w = g.nsOut[v][j];
                if (pruned[w] or deletedOutEdge[v].test(j))
                    continue;
                ui ind = outLookup[w];
                if (ind and not deletedOutEdge[u].test(ind - 1))
                    decrementCN(MM, v, j);
                if (inLookup[w])
                {
                    ind = getLowerBound(g.nsOut[w], u);
                    if (not deletedOutEdge[w].test(ind))
                        decrementCN(MP, v, j);
                }
            }
        }
        for (ui i = 0; i < g.nsIn[u].size(); i++)
        {
            ui v = g.nsIn[u][i];
            // u->v is an edge
            if (pruned[v])
                continue;
            ui uIndv = getLowerBound(g.nsOut[v], u);
            if (deletedOutEdge[v].test(uIndv))
                continue;
            if (not prePruned[v] and outDegree[v]-- + k1 == q)
            {
                Qv.push_back(v);
                prePruned[v] = 1;
                continue;
            }

            for (ui j = 0; j < g.nsOut[v].size(); j++)
            {
                ui w = g.nsOut[v][j];
                if (pruned[w] or deletedOutEdge[v].test(j))
                    continue;
                ui ind = outLookup[w];
                if (ind and not deletedOutEdge[u].test(ind - 1))
                    decrementCN(PM, v, j);
                if (inLookup[w])
                {
                    ind = getLowerBound(g.nsOut[w], u);
                    if (not deletedOutEdge[w].test(ind))
                        decrementCN(PP, v, j);
                }
            }
        }

        // Here calling erase explicitly, because in/out edge list of u is being cleared afterwards...
        inLookup.erase();
        outLookup.erase();

        g.nsIn[u].clear();
        g.nsOut[u].clear();
        // All out-edges are deleted now
        deletedOutEdge[u].setAll();
        // All in-edges are deleted...
        for (ui v : g.nsIn[u])
        {
            ui uIndv = getLowerBound(g.nsOut[v], u);
            deletedOutEdge[v].set(uIndv);
        }

        // cout << Qv.size() << " vertices to be removed" << endl;
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
        for (ui u = 0; u < g.V; u++)
        {

            if (!pruned[u])
            {
                // out nieghbors are pruned by deleted edges as well
                compact(g.nsOut[u]);
                compact(g.nsIn[u]);
            }
        }
    }

    void compactAdjListsWithRemomvedEdges()
    {
        for (ui u = 0; u < g.V; u++)
            g.nsIn[u].clear();

        for (ui u = 0; u < g.V; u++)
        {
            if (pruned[u])
                continue;
            // out nieghbors are pruned by deleted edges as well
            compact(g.nsOut[u], deletedOutEdge[u]);

            for (ui v : g.nsOut[u])
            {
                g.nsIn[v].push_back(u);
            }
        }
        for (ui u = 0; u < g.V; u++)
            sort(g.nsIn[u].begin(), g.nsIn[u].end());
    }

    bool intersectsAll(auto &X, auto &Y)
    {
        // checkes that every element in X is in Y
        // ? do we need to care about pruned vertices? probably not, because pruned vertices can't be in X
        for (auto x : X)
        {
            if (!binary_search(Y.begin(), Y.end(), x))
                return false;
        }
        return true;
    }

    vector<ui> updateC()
    {
        vector<ui> rC;
        rC.reserve(C.size());
        for (ui i = 0; i < C.size(); i++)
        {
            ui u=C[i];
            if (!canMoveToP(u))
            {
                rC.push_back(u);
            }
        }

        for (auto u : rC)
            removeFromC(u);
        return rC;
    }
    vector<ui> updateX()
    {
        vector<ui> rX;
        rX.reserve(X.size());

        for (ui i = 0; i < X.size(); i++)
        {
            ui u=X[i];
            if (!canMoveToP(u))
            {
                rX.push_back(u);
            }
        }
        for (auto u : rX)
                X.remove(u);
        return rX;
    }
    // calculates two-hop iterative pruned graph according to Algo 2
    void getTwoHopIterativePrunedG(ui s)
    {
        in2HopG[s] = 1;
        auto &nsIn = g.nsIn[s];
        auto &nsOut = g.nsOut[s];
        if (nsIn.size() == 0 or nsOut.size() == 0)
            return;
        // auto inLookup = getLookup(nsIn);
        // auto outLookup = getLookup(nsOut);
        Lookup inLookup(look1, nsIn);
        Lookup outLookup(look2, nsOut);
        // both in-out neighbors
        vector<ui> nsInOut;
        nsInOut.reserve(min(nsIn.size(), nsOut.size()));
        // exclusive in/out neighbors
        vector<ui> I;
        vector<ui> O;
        I.reserve(nsIn.size());
        O.reserve(nsOut.size());

        // vector<ui> existsIn(g.V, 0);
        // vector<ui> existsOut(g.V, 0);
        Lookup existsIn(look3, nsIn);
        Lookup existsOut(look4, nsOut);
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
            vector<ui> tempIn, tempOut;
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
                for (auto v : g.nsOut[u])
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
                for (auto v : g.nsIn[u])
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
        // since in/out neighbors are added, now onwards I is SI and O is SO

        // Calculate and B to two hop graph

        for (auto u : O)
        {
            for (auto v : g.nsOut[u])
            {
                // v is not already added in 2hop graph
                if (!in2HopG[v])
                    in2HopG[v] = 2;
            }
        }
        for (auto u : O)
        {
            for (auto v : g.nsIn[u])
            {
                if (in2HopG[v] == 2)
                    in2HopG[v] = 3;
            }
        }
        for (auto u : I)
        {
            for (auto v : g.nsOut[u])
            {
                if (in2HopG[v] == 3)
                    in2HopG[v] = 4;
            }
        }
        for (auto u : I)
        {
            for (auto v : g.nsIn[u])
            {
                if (in2HopG[v] == 4)
                {
                    addTo2HopG(v);
                    // cout<<"added "<<v<<endl;
                }
            }
        }
        // cout<<s<<" : ";
        // print("X ", X);
        // print("C ", C);
    }

    void addTo2HopG(ui u)
    {
        if (in2HopG[u] == 1)
            return;
        in2HopG[u] = 1;
        // cout<<P.front()<<endl;
        if (peelSeq[u] < peelSeq[vi])
        {
            X.add(u);
        }
        else
        {
            addToC(u);
        }
    }
    void print(string msg, auto &vec)
    {
        cout << msg;
        for (auto u : vec)
            cout << u << " ";
        cout << endl;
    }

    void addToC(ui u)
    {
        C.add(u);
        for (ui v : g.nsOut[u])
        {
            dGin[v]++;
        }
        for (ui v : g.nsIn[u])
        {
            dGout[v]++;
        }
    }

    void removeFromC(ui u)
    {
        C.remove(u);
        for (ui v : g.nsOut[u])
        {
            dGin[v]--;
        }
        for (ui v : g.nsIn[u])
        {
            dGout[v]--;
        }
    }

    void PToC(ui u)
    {
        P.remove(u);
        C.add(u);
        for (ui v : g.nsOut[u])
        {
            dPin[v]--;
        }
        for (ui v : g.nsIn[u])
        {
            dPout[v]--;
        }
    }

    void CToP(const ui &u)
    {
        // assert(Cand.contains(u));
        C.remove(u);
        P.add(u);
        for (ui v : g.nsOut[u])
        {
            dPin[v]++;
        }
        for (ui v : g.nsIn[u])
        {
            dPout[v]++;
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
        {
            return false;
        }
        for (ui i = 0; i < P.size(); i++)
        {
            ui v = P.get(i);
            // assert(neiInP[v] + k >= P.size());
            // should be in-connected to every out-boundary vertex
            if (dPout[v] + k1 == P.size() && !binary_search(g.nsIn[u].begin(), g.nsIn[u].end(), v))
            {
                return false;
            }
            // should be out-connected to every in-boundary vertex
            if (dPin[v] + k2 == P.size() && !binary_search(g.nsOut[u].begin(), g.nsOut[u].end(), v))
            {
                return false;
            }
        }
        return true;
    }
};

int main(int argc, char *argv[])
{
    CommandLine cmd(argc, argv);
    std::string data_file = cmd.GetOptionValue("-g", "");
    ui q = cmd.GetOptionIntValue("-q", 2);
    ui k1 = cmd.GetOptionIntValue("-k1", 1);
    ui k2 = cmd.GetOptionIntValue("-k2", 1);

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
    Graph g(data_file);
    cout << "n=" << g.V << endl;
    cout << "Loading Done" << endl;

    cout << "starting kplex" << endl;
    auto tick = chrono::steady_clock::now();

    EnumKPlex kp(g, k1, k2, q);
    kp.enumerate();
    cout << data_file << " Enumeration time: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - tick).count() << " ms" << endl;

    cout << endl;
    return 0;
}