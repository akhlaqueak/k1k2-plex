# Maximal Directed k-Plex Mining

## Serial Programs

### Compile
```
cd serial
make
```

### Execution
```
./kplex -g <dataset.txt> -k1 <k1 value> -k2 <k2 value> -q <min size of (k1,k2)-plex>
```

### Baselines and Ablation Study Versions
All the versions reported as baselines and in ablation study can be produced with:
```
make versions
```

## Parallel Execution

### Compile
```
cd parallel
make
```

### Execution
```
export OMP_NUM_THREADS=<no. of threads>
./kplex -g <dataset.txt> -k1 <k1 value> -k2 <k2 value> -q <min size of (k1,k2)-plex> -t <tau time threshold. default=0.1ms>
```