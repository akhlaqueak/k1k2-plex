#include "graph.h"

void processNode(unsigned int v, Graph &g, unsigned int* degrees, unsigned int* buffer, unsigned int &tail, unsigned int level){

    unsigned int start = g.neighbors_offset[v];
    unsigned int end = g.neighbors_offset[v+1];

    for(unsigned int j = start; j<end; j++){
        
        unsigned int u = g.neighbors[j];

        if(degrees[u] > level){
            degrees[u]--; 
            
            if(degrees[u]==level){
                // add to buffer
                buffer[tail++] = u;
            }
            
        }

    }
           
}



unsigned int* find_kcore(Graph &g){
    unsigned int *degrees = new unsigned int[g.V];
    for(unsigned int i=0;i<g.V;i++){
        degrees[i] = g.degrees[i];
    }

    unsigned int *buffer = new unsigned int[g.V];
    vector<unsigned int> degOrderOffsets;
    degOrderOffsets.push_back(0);

    unsigned int count = 0;
    unsigned int tail = 0;
    auto start = chrono::steady_clock::now();
    for(unsigned int level=0; tail<g.V; level++){

        for(int i=0;i<g.V;i++)
            if(degrees[i] == level)
                buffer[tail++] = i;
    //  cout<<"Total nodes here: "<<tail<<endl;
        for(int i=degOrderOffsets.back();i<tail;i++)
            processNode(buffer[i], g, degrees, buffer, tail, level);
            
        degOrderOffsets.push_back(tail);

        //cout<<"*********Completed level: "<<level<<", global_count: "<<tail<<" *********"<<endl;

    }
    auto end = chrono::steady_clock::now();
    cout << "Elapsed Time: "
    << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
   
   
   
    // degrees sorting after degenracy... 
    unsigned int* rec = new unsigned int[g.V];
    vector<pair<unsigned int, unsigned int>> degOrder(g.V);
    for(int i=0;i<g.V;i++){
        unsigned int v = buffer[i];
        degOrder[i]={v, g.degrees[v]};
        // cout<<g.degrees[v]<<" ";
    }
    // sort each k-shell vertices based on their degrees. 
    auto degComp = [](auto a, auto b){return a.second<b.second;};
    for(int i=0;i<degOrderOffsets.size()-1;i++)
        std::sort(degOrder.begin()+degOrderOffsets[i], degOrder.begin()+degOrderOffsets[i+1], degComp);
    // copy back the sorted vertices to rec array... 
    for(int i=0;i<g.V;i++)
        rec[i] = degOrder[i].first;
    delete[] buffer;
    delete[] degrees;
    return rec;
}