"""Helpers for the author-citation case study of the k1k2-plex paper (IEEE BigData 2026).

Pipeline: download the ArnetMiner citation data -> build the author citation graph ->
run the compiled k1k2-plex program with `-o` -> load, inspect and draw the resulting plexes.
"""
import math
import os
import subprocess
import time
import urllib.request

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(HERE, 'data')
OUTPUT = os.path.join(HERE, 'output')
CITATION_URL = 'https://lfs.aminer.cn/lab-datasets/soinf/citation-raw.txt'
SERIAL_BIN = os.path.join(REPO, 'serial', 'kplex')


def download_citation_data(url=CITATION_URL):
    """Downloads the raw ArnetMiner citation data set (about 0.4 MB) unless it is already present."""
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, 'citation-raw.txt')
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    return path


def build_author_graph(raw_path):
    """Builds the directed author citation network: edge (u, v) if a paper (co-)authored by u cites a
    paper (co-)authored by v (u != v). Writes author_graph.txt, author_graph_undirected.txt (every edge
    in both directions) and authors.tsv (vertex ID -> name) into data/, and returns their paths."""
    papers, cites, section = {}, [], None
    with open(raw_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('*Vertices'):
                section = 'V'
            elif line.startswith('*Edges'):
                section = 'E'
            elif line.strip() and not line.startswith('#'):
                if section == 'V':
                    fields = line.split('\t')
                    papers[int(fields[0])] = [a.strip() for a in fields[4].split(',') if a.strip()] if len(fields) > 4 else []
                elif section == 'E':
                    s, t = line.split()[:2]
                    cites.append((int(s), int(t)))
    edges = {(u, v) for s, t in cites for u in papers.get(s, []) for v in papers.get(t, []) if u != v}
    names = sorted({a for e in edges for a in e})
    idx = {a: i for i, a in enumerate(names)}
    paths = {k: os.path.join(DATA, k) for k in ('author_graph.txt', 'author_graph_undirected.txt', 'authors.tsv')}
    with open(paths['authors.tsv'], 'w') as f:
        f.writelines(f'{i}\t{a}\n' for a, i in idx.items())
    with open(paths['author_graph.txt'], 'w') as f:
        f.writelines(f'{idx[u]} {idx[v]}\n' for u, v in sorted(edges))
    both = {(idx[u], idx[v]) for u, v in edges} | {(idx[v], idx[u]) for u, v in edges}
    with open(paths['author_graph_undirected.txt'], 'w') as f:
        f.writelines(f'{u} {v}\n' for u, v in sorted(both))
    authors = {a for authors_of_paper in papers.values() for a in authors_of_paper}
    print(f'{len(papers)} papers, {len(cites)} citations -> {len(authors)} authors, {len(edges)} citation edges')
    return paths


def run_kplex(graph, k1, k2, q, out_name, binary=SERIAL_BIN):
    """Runs the compiled serial program and returns (path of the -o file, enumeration time in ms)."""
    if not os.path.exists(binary):
        raise FileNotFoundError(f'{binary} not found; run `make` in the repository root first')
    os.makedirs(OUTPUT, exist_ok=True)
    out = os.path.join(OUTPUT, out_name)
    start = time.time()
    log = subprocess.run([binary, '-g', graph, '-k1', str(k1), '-k2', str(k2), '-q', str(q), '-o', out],
                         check=True, capture_output=True, text=True).stdout
    elapsed = [l for l in log.splitlines() if 'Enumeration time (ms)' in l]
    ms = int(elapsed[0].split(':')[-1]) if elapsed else round((time.time() - start) * 1000)
    return out, ms


def load_graph(path):
    """Returns the directed edge set {(u, v)} of a graph file written by build_author_graph()."""
    with open(path) as f:
        return {tuple(map(int, l.split())) for l in f}


def load_names(path):
    with open(path) as f:
        return {int(i): n for i, n in (l.rstrip('\n').split('\t') for l in f)}


def load_plexes(path):
    """Reads the -o output: one plex per line (vertex IDs)."""
    with open(path) as f:
        return [frozenset(map(int, l.split())) for l in f if l.strip()]


def describe(plex, edges, names):
    """Returns one line per member: in-/out-degree inside the plex and the missing directed links."""
    rows = []
    for v in sorted(plex, key=lambda x: names[x]):
        out_deg = sum((v, u) in edges for u in plex)
        in_deg = sum((u, v) in edges for u in plex)
        not_citing = sorted(names[u] for u in plex if u != v and (v, u) not in edges)
        not_cited_by = sorted(names[u] for u in plex if u != v and (u, v) not in edges)
        row = f'  {names[v]:28s} in {in_deg:2d}  out {out_deg:2d}'
        if not_citing:
            row += f'  | does not cite: {", ".join(not_citing)}'
        if not_cited_by:
            row += f'  | not cited by: {", ".join(not_cited_by)}'
        rows.append(row)
    return '\n'.join(rows)


def find(plexes, names, members=(), size=None):
    """Returns the plexes that contain all given author names (and have the given size, if any)."""
    ids = {i for i, n in names.items() if n in set(members)}
    assert len(ids) == len(set(members)), 'unknown author name'
    return [p for p in plexes if ids <= p and (size is None or len(p) == size)]


def draw_plex(plex, edges, names, out_path, highlight=None, order=None):
    """Draws a plex on an ellipse (style of Fig. 7 of the paper). Bidirectional citations get two arrow heads.
    If `highlight` names a member, its incoming edges are drawn in red and its outgoing edges in blue."""
    members = order or sorted(plex, key=lambda v: names[v])
    n, r = len(members), 0.95
    R = 4.2 * max(1, n / 8.5)
    rx, ry = R * 1.45, R * 0.9
    pos = {v: (rx * math.cos(math.pi / 2 + math.pi / n - 2 * math.pi * i / n),
               ry * math.sin(math.pi / 2 + math.pi / n - 2 * math.pi * i / n)) for i, v in enumerate(members)}
    hl = next((v for v in members if names[v] == highlight), None)
    w, h = rx + r + 0.8, ry + r + 0.2
    fig, ax = plt.subplots(figsize=(2 * w * 0.68, 2 * h * 0.68))
    for i, a in enumerate(members):
        for b in members[i + 1:]:
            ab, ba = (a, b) in edges, (b, a) in edges
            for src, dst, both in ([(a, b, True)] if ab and ba and hl not in (a, b) else
                                   [(a, b, False)] * ab + [(b, a, False)] * ba):
                color = 'red' if dst == hl else 'blue' if src == hl else 'black'
                (x1, y1), (x2, y2) = pos[src], pos[dst]
                d = math.hypot(x2 - x1, y2 - y1)
                ux, uy = (x2 - x1) / d, (y2 - y1) / d
                ax.add_patch(FancyArrowPatch((x1 + ux * r, y1 + uy * r), (x2 - ux * r, y2 - uy * r),
                                             arrowstyle='<->' if both else '->', mutation_scale=22, shrinkA=0, shrinkB=0,
                                             lw=2.4 if color != 'black' else 1.4, color=color, zorder=2 if color != 'black' else 1))
    for v in members:
        x, y = pos[v]
        ax.add_patch(Circle((x, y), r, color='#cdc8c8', zorder=3))
        first, _, last = names[v].rpartition(' ')
        ax.text(x, y, f'{first}\n{last}', ha='center', va='center', fontsize=15, zorder=4, linespacing=1.0)
    ax.set_xlim(-w, w)
    ax.set_ylim(-h, h)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.savefig(out_path, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)
    return out_path
