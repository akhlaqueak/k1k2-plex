#ifndef CUTS_GRAPH_H
#define CUTS_GRAPH_H

#include <bits/stdc++.h>
using namespace std;
using namespace chrono;
typedef unsigned int ui;

class Graph
{
public:
    vector<vector<ui>> nsIn, nsOut;
    ui V;
    void readTextFile(string);
    void writeBinFile(string);
    void readBinFile(string);
};
void Graph::readBinFile(string fname)
{
    std::ifstream rdfile(fname, std::ios::binary);
    rdfile.read(reinterpret_cast<char *>(&V), sizeof(ui));
    nsIn.resize(V);
    nsOut.resize(V);
    ui nv=0, ne=0;
    ui maxin = 0, maxout = 0;
    for (ui i = 0; i < V; i++)
    {
        ui m;
        rdfile.read(reinterpret_cast<char *>(&m), sizeof(ui));
        nsOut[i].resize(m);
        if (m > maxout)
            maxout = m;
        rdfile.read(reinterpret_cast<char *>(&nsOut[i][0]), m * sizeof(ui));
        ne+=m;
    }
    for (ui i = 0; i < V; i++)
    {
        ui m;
        rdfile.read(reinterpret_cast<char *>(&m), sizeof(ui));
        nsIn[i].resize(m);
        if (m > maxin)
            maxin = m;
        rdfile.read(reinterpret_cast<char *>(&nsIn[i][0]), m * sizeof(ui));
    }
    cout<<"|V| = "<<V<<endl;
    cout<<"|E| = "<<ne<<endl;
    cout << "Max d+ " << maxout << endl;
    cout << "Max d- " << maxin << endl;
}

void Graph::writeBinFile(string fname)
{
    std::ofstream wfile(fname, std::ios::binary);
    wfile.write(reinterpret_cast<char *>(&V), sizeof(ui));
    for (ui i = 0; i < V; i++)
    {
        ui m = nsOut[i].size();
        wfile.write(reinterpret_cast<char *>(&m), sizeof(ui));
        wfile.write(reinterpret_cast<char *>(&nsOut[i][0]), m * sizeof(ui));
    }
    for (ui i = 0; i < V; i++)
    {
        ui m = nsIn[i].size();
        wfile.write(reinterpret_cast<char *>(&m), sizeof(ui));
        wfile.write(reinterpret_cast<char *>(&nsIn[i][0]), m * sizeof(ui));
    }
}

void Graph::readTextFile(string input_file)
{

    ifstream infile;
    infile.open(input_file);
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
    }

    for (auto &ns : nsOut)
    {
        sort(ns.begin(), ns.end());
        auto last = std::unique(ns.begin(), ns.end());
        ns.erase(last, ns.end());
    }
    ui maxdeg = 0;
    for (auto &ns : nsOut)
        if (ns.size() > maxdeg)
            maxdeg = ns.size();
    cout << "Max Out Degree " << maxdeg << endl;

    maxdeg = 0;
    for (auto &ns : nsIn)
        if (ns.size() > maxdeg)
            maxdeg = ns.size();
    cout << "Max In Degree " << maxdeg << endl;
}

#endif // CUTS_GRAPH_H
