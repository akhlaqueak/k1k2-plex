
#include "./graph.h"
// #define RAPIDS
#include "command_line.h"
class EnumKPlex
{
    Graph g;
    ui kplexes;
    ui k1, k2, q;
    vector<ui> vBoundaryIn;
    vector<ui> vBoundaryOut;
    vector<ui> neighPIn;
    vector<ui> neighPOut;
    vector<ui> P;
    vector<char> pruned; // todo change it to bitset for memory efficiency
    vector<char> inTwoHopG;
    vector<ui> degenOrder;
    vector<ui> rid;
    vector<ui> inDegree;
    vector<ui> outDegree;
    vector<ui> Cinit;
    vector<ui> Xinit;

    void print(string msg, auto &vec)
    {
        cout << msg;
        for (auto u : vec)
            cout << u << " ";
        cout << endl;
    }

public:
    void enumerate()
    {
        // remove vertices v with outdeg(v)+k1 < q OR indeg(v)+k2 < q
        k1k2CorePrune();
        cout << std::flush;
        // find degeneracy order, the result is degenOrder vector
        degenerate();
        // recode IDs so that intersection in following loop is avoided
        recodeIDs();
        // print("degen: ", degenOrder);
        // print("rid: ", rid);
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            // getTwoHopG selects vertices in 2 hop neigbhors of i, and appends them in Xinit, Cinit
            // B = two hop neighbors of i
            // C = vertices u in B such that u<i
            // X = vertices u in B such that u>i
            ui u = degenOrder[i];
            // cout << u << " ... " << endl;
            addToP(u);
            getTwoHopG(u);
            // print("C: ", Cinit);
            // print("X: ", Xinit);
            kplex(Cinit, Xinit);
            reset();
        }
        cout << "Total (" << k1 << "-" << k2 << ")-plexes " << kplexes << endl;
    }
    EnumKPlex(Graph &_g, ui _k1, ui _k2, ui _q) : pruned(_g.V), rid(_g.V),
                                                  g(_g), inDegree(_g.V), outDegree(_g.V),
                                                  inTwoHopG(_g.V), neighPIn(_g.V), neighPOut(_g.V),
                                                  k1(_k1), k2(_k2), q(_q), kplexes(0)
    {
        vBoundaryIn.reserve(_g.V);
        vBoundaryOut.reserve(_g.V);

        P.reserve(_g.V);
        Cinit.reserve(_g.V);
        Xinit.reserve(_g.V);

        degenOrder.reserve(_g.V);
        reset();
    }

private:
    void degenerate()
    {
        ui kc1 = q - k1, kc2 = q - k2, pn = 0, tail = 0;
        vector<char> peeled = pruned;
        cout<<"peeled size "<<peeled.size();
        cout<<"pruned size "<<pruned.size();
        // peeled.resize(g.V);
        cout<<"a peeled size "<<peeled.size();
        // check how many vertices have been pruned by k1k2pruning
        for (ui i = 0; i < g.V; i++)
        {
            if (peeled[i])
                pn++;
        }
        cout << "pruned vertices: " << pn << endl;
        while (pn + degenOrder.size() < g.V)
        {

            for (ui i = 0; i < g.V; i++)
            {
                if (peeled[i])
                    continue;
                if (outDegree[i] == kc1 or inDegree[i] == kc2)
                {
                    degenOrder.push_back(i);
                    peeled[i] = 1;
                }
            }
            for (ui i = tail; i < degenOrder.size(); i++)
            {
                ui v = degenOrder[i];
                for (ui j = 0; j < g.nsIn[v].size(); j++)
                {
                    ui u = g.nsIn[v][j];
                    if (peeled[u])
                        continue;

                    outDegree[u]--;
                    if (outDegree[u] == kc1)
                    {
                        degenOrder.push_back(u);
                        peeled[u] = 1;
                    }
                }
                for (ui j = 0; j < g.nsOut[v].size(); j++)
                {
                    ui u = g.nsOut[v][j];
                    if (peeled[u])
                        continue;

                    inDegree[u]--;
                    if (inDegree[u] == kc2)
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
    }

    void recodeIDs()
    {
        for (ui i = 0; i < degenOrder.size(); i++)
        {
            rid[degenOrder[i]] = i;
        }
    }

    void getNeighbors(auto &neigh)
    {
        // u the first vertex in P is the one for which we are building two-hop neighbors graph
        ui u = P.front();
        for (ui i = 0; i < neigh.size(); i++)
        {
            ui v = neigh[i];
            if (pruned[v] or inTwoHopG[v])
                continue;
            if (rid[v] < rid[u])
            {
                Xinit.push_back(v);
            }
            else
                Cinit.push_back(v);
            inTwoHopG[v] = 1;
        }
    }
    void getSecondHop(auto &neigh, ui c1h, ui x1h)
    {
        for (size_t i = 0; i < c1h; i++)
        {
            getNeighbors(neigh[Cinit[i]]);
        }
        for (size_t i = 0; i < x1h; i++)
        {
            getNeighbors(neigh[Xinit[i]]);
        }
    }

    void getTwoHopG(ui u)
    {
        // reset inTwoHopG to 0

        inTwoHopG[u] = 1;
        // get 1st hop neighbors, they are appended to C and X
        getNeighbors(g.nsOut[u]);
        getNeighbors(g.nsIn[u]);

        // get 2nd hop neighbors, first hop neighbors are found until c1h, x1h
        // 2nd hop neighbors are appended afterwards
        ui c1h = Cinit.size(), x1h = Xinit.size();
        getSecondHop(g.nsOut, c1h, x1h);
        getSecondHop(g.nsIn, c1h, x1h);
    }

    void reset()
    {
        P.clear();
        Cinit.clear();
        Xinit.clear();

        fill(neighPIn.begin(), neighPIn.end(), 0);
        fill(neighPOut.begin(), neighPOut.end(), 0);
        fill(inTwoHopG.begin(), inTwoHopG.end(), 0);

        vBoundaryIn.clear();
        vBoundaryOut.clear();
    }

    void k1k2CorePrune()
    {
        vector<ui> buffer;
        buffer.reserve(g.V);
        for (ui i = 0; i < g.V; i++)
        {
            inDegree[i] = g.nsIn[i].size();
            outDegree[i] = g.nsOut[i].size();
            if (outDegree[i] + k1 < q or inDegree[i] + k2 < q)
            {
                pruned[i] = 1;
                buffer.push_back(i);
            }
        }
        for (ui i = 0; i < buffer.size(); i++)
        {
            ui v = buffer[i];
            for (ui j = 0; j < g.nsIn[v].size(); j++)
            {
                ui u = g.nsIn[v][j];
                if (pruned[u])
                    continue;
                outDegree[u]--;
                if (outDegree[u] + k1 < q)
                {
                    buffer.push_back(u);
                    pruned[u] = 1;
                }
            }
            for (ui j = 0; j < g.nsOut[v].size(); j++)
            {
                ui u = g.nsOut[v][j];
                if (pruned[u])
                    continue;
                inDegree[u]--;
                if (inDegree[u] + k2 < q)
                {
                    buffer.push_back(u);
                    pruned[u] = 1;
                }
            }
        }
    }

    void addToP(ui v)
    {
        P.push_back(v);
        // update out-neighbors
        for (auto u : g.nsOut[v])
        {
            neighPIn[u]++;
        }
        // update in-neighbors
        for (auto u : g.nsIn[v])
        {
            neighPOut[u]++;
        }
        // update in/out-boundary-vertices
        vBoundaryIn.clear();
        vBoundaryOut.clear();
        for (auto u : P)
        {
            if (neighPIn[u] + k2 == P.size())
                vBoundaryIn.push_back(u);
            if (neighPOut[u] + k1 == P.size())
                vBoundaryOut.push_back(u);
        }
        // print("Pin", neighPIn);
        // print("Pout", neighPOut);
        // print("boundIn: ", vBoundaryIn);
        // print("boundOut: ", vBoundaryOut);
    }

    void removeFromP(ui v)
    {
        // v is at top of stack
        P.pop_back();
        // update out-neighbors
        for (auto u : g.nsOut[v])
        {
            neighPIn[u]--;
        }
        // update in-neighbors
        for (auto u : g.nsIn[v])
        {
            neighPOut[u]--;
        }
        // update in/out-boundary-vertices

        vBoundaryIn.clear();
        vBoundaryOut.clear();
        for (auto u : P)
        {
            if (neighPIn[u] + k2 == P.size())
                vBoundaryIn.push_back(u);
            if (neighPOut[u] + k1 == P.size())
                vBoundaryOut.push_back(u);
        }
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

    vector<ui> update(auto &B)
    {
        vector<ui> res1, res2;
        res1.reserve(B.size());
        res2.reserve(B.size());
        for (ui i = 0; i < B.size(); i++)
        {
            // every vertex in vBoundary is in neighbors of B[i]
            if (intersectsAll(vBoundaryIn, g.nsOut[B[i]]))
                // if (P.size() < k or P.intersect(g.nsIn, B[i]).size() > P.size() - k1)
                if (neighPIn[B[i]] + k2 > P.size())
                    res1.push_back(B[i]);
        }

        for (ui i = 0; i < res1.size(); i++)
        {
            // every vertex in vBoundary is in neighbors of B[i]
            if (intersectsAll(vBoundaryOut, g.nsIn[res1[i]]))
                // if (P.size() < k or P.intersect(g.nsOut, B[i]).size() > P.size() - k2)
                if (neighPOut[res1[i]] + k1 > P.size())
                    res2.push_back(res1[i]);
        }

        return res2;
    }
    void kplex(auto &C, auto &X)
    {
        if (C.empty() and X.empty())
        {
            if (P.size() < q)
                return;
            print("kplex: ", P);
            kplexes++;
            return;
        }
        while (!C.empty())
        {
            // doing pop because loop condition is dependent on C
            ui v = C.back();
            C.pop_back();

            addToP(v);
            // print("C: ", C);
            // print("X: ", X);

            vector<ui> C1 = update(C);
            vector<ui> X1 = update(X);

            // print("P: ", P);
            // print("C1: ", C1);
            // print("X1: ", X1);
            kplex(C1, X1);
            removeFromP(v);
            X.push_back(v);
        }
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

    cout << "Loading Started" << endl;
    Graph g(data_file);
    cout << "n=" << g.V << endl;
    cout << "Loading Done" << endl;

    cout << "starting kplex" << endl;
    auto tick = chrono::steady_clock::now();

    EnumKPlex kp(g, k1, k2, q);
    kp.enumerate();
    cout << data_file << " Enumeration time: " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - tick).count() << endl;

    cout << endl;
    return 0;
}