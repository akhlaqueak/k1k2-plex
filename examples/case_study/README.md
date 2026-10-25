# Case study: author citation network

Reproduces the case study of the paper (Section VIII, Fig. 7): $(k_1,k_2)$-plexes reveal cohesive
communities of researchers who cite each other, which are hidden when edge directions are ignored.

## Run

```bash
make                                              # in the repository root (builds serial/kplex)
pip install -r examples/case_study/requirements.txt
python3 examples/case_study/run_case_study.py     # about 15 seconds
```

or open [`case_study.ipynb`](case_study.ipynb), which walks through the same steps and shows the figures inline.

The script

1. downloads the [ArnetMiner citation data set](https://lfs.aminer.cn/lab-datasets/soinf/) into `data/`
   (2,555 papers from 10 topics),
2. builds the **author citation network** with 2,160 authors and 19,185 edges, where $u \to v$ if a paper
   (co-)authored by $u$ cites a paper (co-)authored by $v$,
3. runs `serial/kplex -o` for $(2,2)$-plexes with $q=8$ (118 results) and, on the undirected version of the
   network, for 2-plexes with $q=12$ (65,509 results),
4. prints the plexes of Fig. 7 with their missing citation links and draws them into `output/`.

| Output | Content |
|---|---|
| `output/fig7a_db_plex.pdf` | a $(2,2)$-plex of 10 database researchers |
| `output/fig7b_dm_plex.pdf` | a $(2,2)$-plex of 8 data mining researchers with two one-way citations |
| `output/fig7c_undirected_plex.pdf` | an undirected 2-plex that contains the 8 authors of (b); A. Inkeri Verkamo is cited by 11 members but cites none (red edges) |
| `output/plexes_*.txt` | all mined plexes, one per line (vertex IDs; names in `data/authors.tsv`) |

## Data

The citation data set is provided by the ArnetMiner team for research purposes and is downloaded from its
original location rather than redistributed here. If you use it, please cite

- J. Tang, J. Zhang, L. Yao, J. Li, L. Zhang, and Z. Su. ArnetMiner: Extraction and Mining of Academic Social Networks. KDD 2008.
- J. Tang, J. Sun, C. Wang, and Z. Yang. Social Influence Analysis in Large-scale Networks. KDD 2009.
