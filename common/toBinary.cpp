
#include <bits/stdc++.h>
#include "command_line.h"
using namespace std;
using namespace chrono;
typedef unsigned int ui;
#include "graph.h"
int main(int argc, char *argv[])
{
    CommandLine cmd(argc, argv);
    std::string data_file = cmd.GetOptionValue("-g", "");
    Graph g;
    g.readTextFile(data_file);
    ui ind = data_file.find_last_of(".");
    string binfile = data_file.substr(0, ind) + ".bin";
    g.writeBinFile(binfile);
    return 0;
}
