#!/usr/bin/env python3
"""Generates a small directed test graph with planted dense communities (deterministic for a given seed).

Each community has mostly bidirectional edges plus some one-way and missing edges, so that the
(k1,k2)-plex constraints, branch pruning and look-ahead pruning are all exercised.

Usage: gen_random_graph.py <out.txt> [seed]
"""
import random
import sys

out, seed = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = random.Random(seed)
n, communities = 400, 12
edges = set()
vertices = list(range(n))
for c in range(communities):
    members = rng.sample(vertices, rng.randint(10, 18))
    for u in members:
        for v in members:
            if u != v and rng.random() < 0.88:   # dense but not complete, and not always symmetric
                edges.add((u, v))
for _ in range(6 * n):                           # sparse background edges
    u, v = rng.randrange(n), rng.randrange(n)
    if u != v:
        edges.add((u, v))
with open(out, 'w') as f:
    f.write(f'# planted-community test graph, n={n}, seed={seed}\n')
    for u, v in sorted(edges):
        f.write(f'{u} {v}\n')
print(f'{out}: {len(edges)} edges')
