# Branch and Bound Algorithm
# Setup
```shell
make
```

# Usage
  ./listPlex \<file\> \<k\> \<lb\>

$lb$ is the lower bound.

e.g.

```bash
./listPlex ./dataset/jazz.bin 4 12
./listPlex ./dataset/Email-EuAll.bin 3 12
./listPlex ./dataset/Email-EuAll.bin 4 20
```

# Format

The input graph should be a binary format.
One can convert an edge list format graph file (SNAP format) into this binary format by a converter `toBin` .

usage:

  ./toBin \<input\> \<output\>

e.g.

```bash
./toBin ./dataset/jazz.txt ./datatset/jazz.bin
```

