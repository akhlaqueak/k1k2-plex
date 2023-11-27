#include "../common/utils.h"
#include "../common/command_line.h"
#define PuCSize (P.size() + C.size())

// one of these three options should be selected, deciding how Gi is calculated
// #define ITERATIVE_PRUNE
// #define TWO_HOP
#define NAIVE

// #define BRANCHING
// #define LOOKAHEAD
#define CTCP

#define TIME_NOW chrono::steady_clock::now()
// cutoff time is in minutes
#define CUTOFF_TIME (360)
auto clk = TIME_NOW;
#define CUTOFF (chrono::duration_cast<chrono::minutes>(TIME_NOW - clk).count() > CUTOFF_TIME)

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

class EnumKPlex
{
    Graph &g;
    ui kplexes;
    ui k1, k2, q;

    vector<ui> dPin;
    vector<ui> dPout;

    // G is a graph induced by P \cup C
    vector<ui> dGin;
    vector<ui> dGout;

    vector<char> pruned;    // todo change it to bitset for memory efficiency
    vector<char> prePruned; // todo change it to bitset for memory efficiency
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
    pair<ui, ui> ctcpPr;


    vector<ui> look1, look2, look3, look4;

    ui vi; // current vertex in degeneracy order for which we are searching kplex

    RandList C;
    RandList X;
    RandList P;
    RandList block;
    vector<vector<ui>> giIn, giOut; // two-hop graph
    vector<vector<ui>> GIn, GOut;   // shrinked graph

    vector<ui> rC, rX;

    ui p;
    Direction dir;
    ui vn;

public:
    void enumerate()
    {
        auto prep = TIME_NOW;
        cout<<"core pruning started"<<endl;
        // remove vertices v with outdeg(v)+k1 < q OR indeg(v)+k2 < q
        k1k2CorePrune();
        // find degeneracy order, the result is degenOrder vector
        cout<<"degeneracy order started"<<endl;
        degenerate();

#ifdef CTCP
        auto tick = TIME_NOW;
        cout<<"ctcp started"<<endl;
        applyCoreTrussPruning();
        cout << "ctcp cost (ms): " << chrono::duration_cast<chrono::milliseconds>(TIME_NOW - tick).count() << endl;
#endif
        // some nodes are pruned in k1k2CorePrune, so graph will be shrinked in case ctcp is not applied
        shrinkGraph();
        cout << "preprocessing cost (ms): " << chrono::duration_cast<chrono::milliseconds>(TIME_NOW - prep).count() << endl;
        init(); // initilize G, and other vectors

        ui pruningCost = 0;

        if (GOut.size() < q)
            return;

        for (ui i = 0; i < GOut.size() - q + 1; i++)
        {
            vi = i;
#ifdef ITERATIVE_PRUNE
            auto t1 = TIME_NOW;
            getTwoHopIterativePrunedG(vi);
            pruningCost += chrono::duration_cast<chrono::microseconds>(TIME_NOW - t1).count();
            recurSearch(0);
#else
#ifdef TWO_HOP
            auto t1 = TIME_NOW;
            getTwoHopG(vi);
            pruningCost += chrono::duration_cast<chrono::microseconds>(TIME_NOW - t1).count();
            recurSearch(0);
#else
#ifdef NAIVE
            buildBlockFromG();
            recurSearch(vi);
#endif
#endif
#endif
            auto t2 = TIME_NOW;
            // start from first vertex
            reset(); // clears C and X
        }
        cout << "iterative pruning cost (ms): " << pruningCost / 1000 << endl;
        cout << "Total (" << k1 << "," << k2 << ")-plexes of at least " << q << " size: " << kplexes << endl;
        for (ui i = 0; i < counts.size(); i++)
            if (counts[i])
                cout << "kplexes of size: " << i + 1 << " = " << counts[i] << endl;
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
        if (CUTOFF)
            return;
        if (PuCSize < q)
            return;
        if (C.empty())
        {
            if (X.empty())
                reportSolution();
            return;
        }
#ifdef BRANCHING
        vector<ui> MOut, MIn;
        MOut.reserve(P.size());
        MIn.reserve(P.size());
        calculateM(MOut, MIn);

        if (MOut.size() or MIn.size())
        {
            Direction dir;
            ui vp = pickvp(MOut, MIn, dir);
            multiRecurSearch(vp, dir);
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

        recurSearch(vc);
        CToX(vc);
        branch();
        // recover
        XToC(vc);
        // other branch where P contains u
    }

    void multiRecurSearch(ui vp, Direction dir)
    {
        if (PuCSize < q)
            return;
        ui p;
        vector<ui> vpNN; // It stores {u1, u2, ..., ud} vertices
        vpNN.reserve(C.size());

        auto getNonNeigh = [&](auto &adj)
        {
            for (ui i = 0; i < C.size(); i++)
            {
                ui u = C[i];
                if (!binary_search(adj.begin(), adj.end(), u))
                    vpNN.push_back(u);
            }
        };
        if (dir == Out)
        {
            getNonNeigh(giOut[vp]);
            p = k1 - (P.size() - dPout[vp]);
        }
        else
        {
            getNonNeigh(giIn[vp]);
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
    void calculateM(vector<ui> &MOut, vector<ui> &MIn)
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
    ui pickvp(vector<ui> &MOut, vector<ui> &MIn, Direction &dir)
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
        vector<ui> ctemp = C.getData();

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
        counts[P.size() - 1]++;
        // print("KPlex: ", P);
        kplexes++;
        return;
    }

    EnumKPlex(Graph &_g, ui _k1, ui _k2, ui _q) : pruned(_g.V), peelSeq(_g.V, -1),
                                                  g(_g), inDegree(_g.V), outDegree(_g.V),

                                                  k1(_k1), k2(_k2), q(_q), kplexes(0),
                                                  deletedOutEdge(_g.V),
                                                  cnPP(_g.V), cnPM(_g.V),
                                                  cnMP(_g.V), cnMM(_g.V),
                                                  look1(_g.V), look2(_g.V),
                                                  look3(_g.V), look4(_g.V),
                                                  counts(1000)
    {

        rC.reserve(g.V);
        rX.reserve(g.V);

        for (ui i = 0; i < g.V; i++)
        {
            deletedOutEdge[i] = MBitSet(g.nsOut[i].size());
        }
        Qe.reserve(g.V);
        // reset();
    }
    void getNeighbors(auto &neigh)
    {
        for (ui i = 0; i < neigh.size(); i++)
            addTo2HopG(neigh[i]);
    }

    void getTwoHopG(ui u)
    {
        addTo2HopG(u);

        getNeighbors(GOut[u]);
        getNeighbors(GIn[u]);

        // get 2nd hop neighbors, first hop neighbors are already placed in block
        // 2nd hop neighbors are appended afterwards
        ui sz = block.size();
        // block[0] is u, so not including it for 2hop
        for (size_t i = 1; i < sz; i++)
        {
            ui v = block[i];
            getNeighbors(GOut[v]);
            getNeighbors(GIn[v]);
        }
        buildBlock();
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
        cout<<"CTCP pruned vertices: "<<ctcpPr.first<<endl;
        cout<<"CTCP pruned edges: "<<ctcpPr.second<<endl;
    }

    void deleteEdge(ui u, ui vIndu)
    {
        // the edge might already have been deleted...
        if (deletedOutEdge[u].test(vIndu))
            return;
        ctcpPr.second++;
        deletedOutEdge[u].set(vIndu);

        ui v = g.nsOut[u][vIndu];
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
            // if (vInd >= g.nsOut[u].size())
            // {
            //     cout << u << " " << vInd << " " << g.nsOut[u].size() << endl;
            // }
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
                        cnMP[u][j]++;
                }
            }
        }

        for (ui u = 0; u < g.V; u++)
        {
            for (ui j = 0; j < g.nsOut[u].size(); j++)
            {

                if (pruned[u] or pruned[g.nsOut[u][j]])
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

    void reset()
    {
        ui sz = C.size();
        for (ui i = 0; i < sz; i++)
            removeFromC(C[0]);

        X.clear();
        block.clear();
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
                if (outDegree[u] + k1 < q and not pruned[u])
                {
                    Qv.push_back(u);
                    pruned[u] = 1;
                }
            }
            for (ui u : g.nsOut[v])
            {
                inDegree[u]--;
                if (inDegree[u] + k2 < q and not pruned[u])
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
        ctcpPr.first++;
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
            // adding now
            deletedOutEdge[u].set(i);
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
            deletedOutEdge[v].set(uIndv);
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
    void init()
    {
        ui ds = GOut.size();
        dPin.resize(ds);
        dPout.resize(ds);
        dGin.resize(ds);
        dGout.resize(ds);

        C.init(ds);
        X.init(ds);
        P.init(ds);
        block.init(ds);

#ifdef NAIVE
        giIn = GIn;
        giOut = GOut;
#endif
    }
    void shrinkGraph()
    {
        ui k = 0;
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            ui u = degenOrder[i];
            if (pruned[u])
            {
                continue;
            }
            peelSeq[u] = k;
            degenOrder[k] = u;
            k++;
        }
        cout << "vertices in shrinked graph: " << k << endl;
        degenOrder.resize(k);
        GOut.resize(k);
        GIn.resize(k);
        for (ui u = 0; u < g.V; u++)
        {
            if (pruned[u])
                continue;

            for (ui j = 0; j < g.nsOut[u].size(); j++)
            {
                ui v = g.nsOut[u][j];
                if (pruned[v] or deletedOutEdge[u].test(j))
                    continue;
                ui ru = peelSeq[u];
                ui rv = peelSeq[v];
                // G is degeneracy ordered graph...
                GOut[ru].push_back(rv);
                GIn[rv].push_back(ru);
            }
        }

        for (auto &adj : GOut)
        {
            sort(adj.begin(), adj.end());
            // print("out: ", adj);
        }
        for (auto &adj : GIn)
        {
            sort(adj.begin(), adj.end());
            // print("in: ", adj);
        }
        cout << "Shrinked graph |V|: " << GOut.size() << endl;
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
    // calculates two-hop iterative pruned graph according to Algo 2
    void getTwoHopIterativePrunedG(ui s)
    {
        auto &nsIn = GIn[s];
        auto &nsOut = GOut[s];
        addTo2HopG(s);
        // auto inLookup = getLookup(nsIn);
        // auto outLookup = getLookup(nsOut);
        Lookup inLookup(look1, nsIn);
        Lookup outLookup(look2, nsOut);
        // both in-out neighbors
        vector<ui> nsInOut;
        // exclusive in/out neighbors
        vector<ui> I;
        vector<ui> O;
        nsInOut.reserve(min(nsIn.size(), nsOut.size()));
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
                for (auto v : GOut[u])
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
                for (auto v : GIn[u])
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
        vector<ui> temp;
        Lookup intersect(look4, temp);
        for (auto u : O)
        {
            for (auto v : GOut[u])
            {
                // v is not already added in 2hop graph
                if (!inBlock(v))
                {
                    temp.push_back(v);
                    intersect[v] = 2;
                }
            }
        }
        for (auto u : O)
        {
            for (auto v : GIn[u])
            {
                if (intersect[v] == 2)
                {
                    intersect[v] = 3;
                }
            }
        }
        for (auto u : I)
        {
            for (auto v : GOut[u])
            {
                if (intersect[v] == 3)
                    intersect[v] = 4;
            }
        }
        for (auto u : I)
        {
            for (auto v : GIn[u])
            {
                if (intersect[v] == 4)
                {
                    addTo2HopG(v);
                }
            }
        }

        intersect.erase();
        buildBlock();
    }
    bool inBlock(ui u)
    {
        return block.contains(u);
    }
    void addTo2HopG(ui u)
    {
        if (inBlock(u))
            return;
        block.add(u);
    }

    void buildBlockFromG()
    {
        for (ui i = 0; i < giOut.size(); i++)
        {
            if (i < vi)
                X.add(i);
            else
                addToC(i);
        }
    }

    void buildBlock()
    {
        giIn.resize(block.size());
        giOut.resize(block.size());
        for (auto &adj : giIn)
        {
            adj.clear();
            adj.reserve(block.size());
        }
        for (auto &adj : giOut)
        {
            adj.reserve(block.size());
            adj.clear();
        }
        for (ui i = 0; i < block.size(); i++)
        {
            ui u = block[i];
            for (ui v : GOut[u])
            {
                if (inBlock(v))
                    giOut[i].push_back(block.getIndex(v));
            }

            for (ui v : GIn[u])
            {
                if (inBlock(v))
                    giIn[i].push_back(block.getIndex(v));
            }
        }
        for (auto &adj : giIn)
            sort(adj.begin(), adj.end());
        for (auto &adj : giOut)
            sort(adj.begin(), adj.end());
        for (ui i = 0; i < block.size(); i++)
        {
            ui u = block[i];
            if (u < vi)
                X.add(i);
            else
                addToC(i);
        }
        // cout<<P.front()<<endl;
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
        for (ui v : giOut[u])
            dGin[v]++;
        for (ui v : giIn[u])
            dGout[v]++;
    }

    void removeFromC(ui u)
    {
        C.remove(u);
        for (ui v : giOut[u])
            dGin[v]--;
        for (ui v : giIn[u])
            dGout[v]--;
    }

    void PToC(ui u)
    {
        P.remove(u);
        C.add(u);
        for (ui v : giOut[u])
            dPin[v]--;
        for (ui v : giIn[u])
            dPout[v]--;
    }

    void CToP(const ui &u)
    {
        // assert(Cand.contains(u));
        C.remove(u);
        P.add(u);
        for (ui v : giOut[u])
            dPin[v]++;
        for (ui v : giIn[u])
            dPout[v]++;
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
            // ui ru = recode[u];
            // ui rv = recode[v];
            // should be in-connected to every out-boundary vertex
            // if (dPout[v] + k1 == P.size() && !edgeIn[ru].test(rv))
            if (dPout[v] + k1 == P.size() && !binary_search(giIn[u].begin(), giIn[u].end(), v))
                return false;
            // should be out-connected to every in-boundary vertex
            if (dPin[v] + k2 == P.size() && !binary_search(giOut[u].begin(), giOut[u].end(), v))
                return false;
        }
        return true;
    }
};

int main(int argc, char *argv[])
{
    CommandLine cmd(argc, argv);
    std::string file = cmd.GetOptionValue("-g", "");
    ui q = cmd.GetOptionIntValue("-q", 2);
    ui k1 = cmd.GetOptionIntValue("-k1", 1);
    ui k2 = cmd.GetOptionIntValue("-k2", 1);
    ui sr = cmd.GetOptionIntValue("-sr", 1);

    if (file == "")
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
    Graph g;
    size_t ind = file.find_last_of(".");
    string ext = file.substr(ind, file.size());
    if (ext == string(".bin"))
        g.readBinFile(file);
    else
    {
        g.readTextFile(file);
        string binfile = file.substr(0, ind) + ".bin";
        g.writeBinFile(binfile);
    }
    cout << "n=" << g.V << endl;
    cout << "Loading Done" << endl;

    cout << "starting kplex" << endl;
    auto tick = TIME_NOW;

    EnumKPlex kp(g, k1, k2, q);
    kp.enumerate();
    if (CUTOFF)
        cout << sr << " " << file << " Timed Out" << endl;
    else
        cout << sr << " " << file << " Enumeration time (ms): " << chrono::duration_cast<chrono::milliseconds>(TIME_NOW - tick).count() << endl;

    cout << endl;
    return 0;
}