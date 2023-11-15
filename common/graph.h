#ifndef CUTS_GRAPH_H
#define CUTS_GRAPH_H

#include <bits/stdc++.h>
using namespace std;
typedef unsigned int ui;

#define DS_LOC string("")
#define OUTPUT_LOC string("")

class Graph
{
public:
    vector<vector<ui>> nsIn, nsOut;
    ui V;
    ui E;
    ui AVG_DEGREE = 0;
    Graph(std::string input_file);
    Graph();
    void readFile(string input_file);
};

bool degComp(const pair<ui, ui> &lhs, const pair<ui, ui> &rhs)
{
    return lhs.second > rhs.second;
}

Graph::Graph()
{
    // default constructor
}
void Graph::readFile(string input_file)
{

    ifstream infile;
    infile.open(DS_LOC + input_file);
    if (!infile)
    {
        cout << "load graph file failed " << endl;
        exit(-1);
    }

    ui s, t;

    /**
     * @brief Dataset format:
     * # Number of nodes
     * source destination
     * source destination
     * source destination
     * source destination
     *
     */
    // read number of nodes...
    string line;
    vector<pair<ui, ui>> lines;

    V = 0;
    E = 0;
    while (std::getline(infile, line))
    {
        if (!isdigit(line[0]))
            continue; // it's a comment
        std::istringstream iss(line);
        iss >> s >> t;
        if (s == t)
            continue; // remove self loops
        V = max(s, V);
        V = max(t, V);
        lines.push_back({s, t});
    }
    infile.close();

    // todo recode graph to avoid sparse nodes

    V++; // vertices index starts from 0, so add 1 to number of vertices.
    nsIn.resize(V);
    nsOut.resize(V);

    for (auto &p : lines)
    {
        nsOut[p.first].push_back(p.second);
        nsIn[p.second].push_back(p.first);
        // trying undirected graphs... 
        // nsOut[p.second].push_back(p.first);
        // nsIn[p.first].push_back(p.second);
    }

    lines.clear();

    for (auto &ns : nsIn)
    {
        sort(ns.begin(), ns.end());
        auto last = std::unique(ns.begin(), ns.end());
        ns.erase(last, ns.end());
        E += ns.size();
    }

    for (auto &ns : nsOut)
    {
        sort(ns.begin(), ns.end());
        auto last = std::unique(ns.begin(), ns.end());
        ns.erase(last, ns.end());
    }
}

Graph::Graph(std::string input_file)
{

    auto start = chrono::steady_clock::now();
    readFile(input_file);
    auto end = chrono::steady_clock::now();

    cout << "File Loaded in: " << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << endl;
}

// void Graph::recode(Graph &g, ui* rec){
//     degrees = new ui[g.V];
//     neighbors = new ui[g.E];
//     neighbors_offset = new ui[g.V+1];
//     V = g.V;
//     E = g.E;
//     auto tick = chrono::steady_clock::now();
//     cout<<"Degrees copied"<<endl;
//     for(int i=0;i<g.V;i++)
//         degrees[i] = g.degrees[rec[i]];
//     map<ui, ui> recMapping;
//     for(int i=0;i<g.V;i++)
//         recMapping[rec[i]] = i;

//     neighbors_offset[0] = 0;
//     std::partial_sum(degrees, degrees+V, neighbors_offset+1);

//     for(int v=0;v<V;v++){
//         ui recv = rec[v];
//         ui start = neighbors_offset[v];
//         ui end = neighbors_offset[v+1];
//         for (int j=g.neighbors_offset[recv], k=start; j<g.neighbors_offset[recv+1]; j++, k++){
//             neighbors[k] = recMapping[g.neighbors[j]];
//         }
//         std::sort(neighbors+start, neighbors+end);

//     }
//     cout<<"Reordering Time: "<<chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-tick).count()<<endl;
// }
// Graph::~Graph(){
//     cout<<"Deallocated... "<<endl;
//     delete [] neighbors;
//     delete [] neighbors_offset;
//     delete [] degrees;
// }
#endif // CUTS_GRAPH_H
