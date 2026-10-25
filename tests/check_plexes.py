#!/usr/bin/env python3
"""Independent checker for the output of `kplex -o`.

Verifies against the input graph that every reported set is a (k1,k2)-plex with at least q vertices,
that it is maximal, and that no set is reported twice.

Usage: check_plexes.py <graph.txt> <k1> <k2> <q> <plexes.txt>
"""
import sys


def read_graph(path):
    out, inn = {}, {}
    with open(path) as f:
        for line in f:
            if not line[:1].isdigit():
                continue                      # comment line, same rule as the C++ loader
            s, t = map(int, line.split()[:2])
            if s == t:
                continue                      # self-loops are ignored
            out.setdefault(s, set()).add(t)
            inn.setdefault(t, set()).add(s)
    return out, inn


def main():
    path, k1, k2, q, res = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    out, inn = read_graph(path)
    empty = set()
    O = lambda v: out.get(v, empty)
    I = lambda v: inn.get(v, empty)

    def is_plex(S):
        return all(len(S - O(v)) <= k1 and len(S - I(v)) <= k2 for v in S)

    plexes = [frozenset(map(int, l.split())) for l in open(res) if l.strip()]
    errors = []
    if len(set(plexes)) != len(plexes):
        errors.append(f'{len(plexes) - len(set(plexes))} duplicate sets')
    for S in set(plexes):
        if len(S) < q:
            errors.append(f'set smaller than q: {sorted(S)}')
        elif not is_plex(S):
            errors.append(f'not a ({k1},{k2})-plex: {sorted(S)}')
        else:
            # any vertex that can extend S must be adjacent to S, so only neighbors need to be tried
            cand = set().union(*(O(v) | I(v) for v in S)) - S
            ext = next((w for w in cand if is_plex(S | {w})), None)
            if ext is not None:
                errors.append(f'not maximal (can add {ext}): {sorted(S)}')
    if errors:
        print(f'FAIL: {len(errors)} problem(s); first: {errors[0]}')
        sys.exit(1)
    print(f'OK: {len(plexes)} maximal ({k1},{k2})-plexes with at least {q} vertices')


if __name__ == '__main__':
    main()
