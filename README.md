# k1k2-plex

 This version performs:
 * (k1,k2)-core pruning
 * Degeneracy Order
 * FaPlex Algorithm 1 style search of (k1,k2)-plex

## Compilation
```
g++ graph.cpp kplex.cpp -o kplex
```
## Usage:
```
./kplex -g <filename> -q <least size> -k1 <k1 value> -k2 <k2 value>
```
