#include <bits/stdc++.h>
using namespace std;
typedef unsigned int ui;

#define DS_LOC string("")
#define OUTPUT_LOC string("")

#ifndef CUTS_GRAPH_H
#define CUTS_GRAPH_H

class Graph{
public:
    // vector<vector<int>>
    vector<vector<ui>> nsIn, nsOut;
    ui V;
    ui AVG_DEGREE = 0;
    Graph(std::string input_file);
    Graph();
    void readFile(string input_file);
};
#endif //CUTS_GRAPH_H
