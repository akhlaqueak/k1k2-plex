# Maximal Directed k-Plex Mining

Code for the paper

> **Maximal Directed $k$-Plex Mining**
> Akhlaque Ahmad, Da Yan, Jiao Han, Jalal Khalil, Lyuheng Yuan, Saugat Adhikari, Yang Zhou.
> *IEEE International Conference on Big Data (IEEE BigData)*, 2026.
> Online appendix: [`appendix.pdf`](appendix.pdf)

A **$(k_1,k_2)$-plex** generalizes the $k$-plex to directed graphs: in a vertex set $P$, every vertex has at
most $k_1$ non-out-neighbors and at most $k_2$ non-in-neighbors in $P$ (counting itself). This repository
enumerates all **maximal** $(k_1,k_2)$-plexes with at least $q$ vertices, using a Bron–Kerbosch style
search with

- **SOP**: first- and second-order pruning of vertices and edges (Section IV),
- **DTP**: directional two-hop pruning of the candidates of each vertex (Section V, Algorithm 2),
- **BRA** and **LAP**: branch pruning and look-ahead pruning (Section VI, Algorithm 4),
- a **parallel** version with OpenMP and timeout-based task decomposition (Section VII).

## Repository layout

| Path | Content |
|---|---|
| `serial/` | serial program (Ours), baselines and ablation variants (`kplex.cpp`, `Makefile`) |
| `parallel/` | parallel OpenMP program (`kplex.cpp`, `Makefile`) |
| `common/` | graph loader, containers and command-line parsing shared by both programs |
| `tests/` | correctness tests (`make test`) |
| `examples/case_study/` | the author-citation case study of Section VIII, as a script and a Jupyter notebook |
| `datasets/` | small example graphs |
| `*/*.job` | Slurm scripts used for the original experiments |

## Build

Requirements: a C++20 compiler (GCC ≥ 10 or Clang ≥ 13) and, for the parallel program, OpenMP.

```bash
make            # builds serial/kplex and parallel/kplex
make versions   # also builds the baselines and ablation variants (see below)
make test       # runs the correctness tests (about one minute)
```

On macOS, Apple Clang needs `libomp` for the parallel program (e.g., `brew install libomp`):

```bash
make OMPFLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp"
```

## Usage

```bash
# serial
./serial/kplex -g datasets/congress.txt -k1 2 -k2 3 -q 8

# parallel, with the timeout threshold tau_time in microseconds (default 100 = 0.1 ms)
OMP_NUM_THREADS=32 ./parallel/kplex -g datasets/congress.txt -k1 2 -k2 3 -q 8 -t 100

# write all maximal plexes to a file
./serial/kplex -g datasets/congress.txt -k1 2 -k2 3 -q 8 -o plexes.txt
```

| Option | Meaning |
|---|---|
| `-g <file>` | input graph (text edge list or the `.bin` file created from it) |
| `-k1 <int>`, `-k2 <int>` | maximum number of missing out-/in-neighbors per vertex, counting the vertex itself ($k_1,k_2 \ge 2$) |
| `-q <int>` | minimum plex size; requires $q \ge 2k_1-1$ and $q \ge 2k_2-1$ (Section V) |
| `-o <file>` | optional: write each maximal plex as one line of sorted vertex IDs |
| `-t <int>` | parallel only: timeout threshold $\tau_{time}$ in microseconds |

**Input format.** One directed edge `src dst` per line with non-negative integer vertex IDs; additional
columns are ignored, lines that do not start with a digit are treated as comments, and self-loops and
duplicate edges are removed. On first use, a binary copy `<name>.bin` is written next to the input file and
can be passed to `-g` directly to skip parsing.

**Output.** The program prints the number of maximal $(k_1,k_2)$-plexes with at least $q$ vertices, their
size distribution, and the time of the preprocessing (`preprocessing cost` in the serial program, `ctcp cost`
in the parallel one), of the enumeration (`search cost`) and in total. The times exclude graph loading.
With `-o`, the plexes are written using the vertex IDs of the input file.

## Reproducing the experiments

**Datasets.** The 15 graphs of Table I are publicly available:

| Graphs | Source |
|---|---|
| wiki-vote, as-caida (`as-caida20071105`), epinions (`soc-Epinions1`), email-euall, amazon0505, web-google, soc-pokec, wiki-talk, math-overflow (`sx-mathoverflow`) | [SNAP](https://snap.stanford.edu/data/) (edge lists, usable directly) |
| bitcoin (`soc-sign-bitcoinotc`) | [SNAP](https://snap.stanford.edu/data/soc-sign-bitcoin-otc.html) (CSV; replace commas by spaces) |
| arabic-2005, uk-2005, it-2004, webbase-2001 | [LAW](https://law.di.unimi.it/datasets.php) (WebGraph format, see below) |
| clue-web (`web-ClueWeb09-50m`) | [Network Repository](https://networkrepository.com/web-ClueWeb09-50m.php) |

The LAW graphs are distributed in the compressed WebGraph format. With the
[WebGraph](https://webgraph.di.unimi.it/) jars (e.g., `it.unimi.dsi:webgraph` and its dependencies from Maven
Central) on the class path, download `<name>.graph` and `<name>.properties` and convert them into an edge list:

```bash
java -Xmx32g -cp "jars/*" it.unimi.dsi.webgraph.ArcListASCIIGraph -g BVGraph arabic-2005 arabic-2005.txt
```

**Algorithms.** `make versions` builds the programs compared in the paper:

| Paper | Binary | Techniques |
|---|---|---|
| Base-0 | `serial/kplex-v0` | core pruning only, Bron–Kerbosch (Algorithm 1) |
| Base-1 | `serial/kplex-v1` | + SOP |
| Base-2 | `serial/kplex-v2` | + undirected two-hop neighborhood |
| Base-3 | `serial/kplex-bk` | two-hop neighborhood replaced by DTP |
| Ours | `serial/kplex` (= `kplex-full`) | + BRA and LAP (Algorithm 4) |
| Ours\SOP, Ours\DTP, Ours\BRA, Ours\LAP | `serial/kplex-no-ctcp`, `-no-dtp`, `-no-branching`, `-no-lookahead` | Ours with one technique disabled |

Runs are stopped after 6 hours (serial) or 3 hours (parallel); `CUTOFF_TIME` in `kplex.cpp` sets the limit.

## Tests

`make test` runs `tests/run_tests.sh` on `datasets/congress.txt` and a small graph with planted dense
communities (11 configurations, including $k_1 > k_2$). For each configuration it checks that

1. Ours returns exactly the result set of plain Bron–Kerbosch without pruning (Base-0),
2. the parallel program with 1, 4 and 8 threads returns the same result set,
3. an independent Python checker (`tests/check_plexes.py`) confirms that every set is a maximal
   $(k_1,k_2)$-plex with at least $q$ vertices and that no set is reported twice,
4. the number of plexes matches `tests/expected_counts.txt`.

## Case study

[`examples/case_study`](examples/case_study) reproduces the author-citation case study of the paper
(Fig. 7) in about 15 seconds: it downloads the ArnetMiner citation data, builds the author citation network,
mines $(2,2)$-plexes and undirected 2-plexes with `kplex -o`, and draws the communities shown in the paper.
See [`case_study.ipynb`](examples/case_study/case_study.ipynb) for a walkthrough with the figures inline.

## Changes for the camera-ready version

The following fixes were made after the submission; the results in the camera-ready paper were produced with
this version (see the commit history for details).

- **Branch pruning** no longer reports duplicate or non-maximal plexes when a vertex cannot be added to $P$
  (the remaining branches of Eqs. (2)–(3) are now skipped, as in the paper).
- **Second-order pruning** uses the $cn^{--}$ counter for Theorem 6 in the initial pass; previously plexes
  could be missed when $k_1 > k_2$.
- **Plexes larger than 1000 vertices** no longer overflow the size histogram.
- **DTP** is implemented as the single-pass Algorithm 2 of the paper.
- **Parallel**: the branches of branch pruning are now also decomposed into tasks on timeout (Section VII).
- **New**: `-o` output, a portable build (Clang/macOS), correctness tests and the case study example.

## Citation

```bibtex
@inproceedings{ahmad2026maximal,
  title     = {Maximal Directed $k$-Plex Mining},
  author    = {Ahmad, Akhlaque and Yan, Da and Han, Jiao and Khalil, Jalal and Yuan, Lyuheng and
               Adhikari, Saugat and Zhou, Yang},
  booktitle = {IEEE International Conference on Big Data (BigData)},
  year      = {2026}
}
```
