#!/usr/bin/env python3
"""Reproduces the case study (Section VIII, Fig. 7) of the k1k2-plex paper end to end.

    make                                   # in the repository root
    python3 examples/case_study/run_case_study.py

Downloads the ArnetMiner citation data, builds the author citation network, mines (2,2)-plexes with q = 8 on
the directed graph and 2-plexes with q = 12 on its undirected (bidirected) version, prints the plexes shown
in the paper and draws them into examples/case_study/output/.
"""
import os

from plexlib import (OUTPUT, build_author_graph, describe, download_citation_data, draw_plex, find,
                     load_graph, load_names, load_plexes, run_kplex)

DM_PLEX = ['Rakesh Agrawal', 'Jiawei Han', 'Raymond T. Ng', 'Philip S. Yu',
           'Ramakrishnan Srikant', 'Ming-Syan Chen', 'Yongjian Fu', 'Jong Soo Park']


def main():
    raw = download_citation_data()
    paths = build_author_graph(raw)
    names = load_names(paths['authors.tsv'])
    edges = load_graph(paths['author_graph.txt'])
    ids = {n: i for i, n in names.items()}

    directed, ms = run_kplex(paths['author_graph.txt'], 2, 2, 8, 'plexes_directed_k2-2_q8.txt')
    plexes = load_plexes(directed)
    print(f'\n(2,2)-plexes with q = 8 on the directed graph: {len(plexes)} found in {ms} ms')

    undirected, ms = run_kplex(paths['author_graph_undirected.txt'], 2, 2, 12, 'plexes_undirected_k2_q12.txt')
    uplexes = load_plexes(undirected)
    print(f'2-plexes with q = 12 on the undirected graph: {len(uplexes)} found in {ms} ms')

    # Fig. 7(a): a (2,2)-plex of database researchers
    (a,) = find(plexes, names, ['Anthony Tomasic', 'Jennifer Widom'], size=10)
    print('\nFig. 7(a): (2,2)-plex on DB\n' + describe(a, edges, names))
    draw_plex(a, edges, names, os.path.join(OUTPUT, 'fig7a_db_plex.pdf'))

    # Fig. 7(b): a (2,2)-plex of data mining researchers with two one-way citations
    (b,) = find(plexes, names, DM_PLEX, size=len(DM_PLEX))
    print('\nFig. 7(b): (2,2)-plex on DM\n' + describe(b, edges, names))
    draw_plex(b, edges, names, os.path.join(OUTPUT, 'fig7b_dm_plex.pdf'), order=[ids[n] for n in DM_PLEX])

    # Fig. 7(c): the undirected 2-plexes containing the DM plex are larger, so the DM plex is not maximal there
    containing = find(uplexes, names, DM_PLEX)
    sizes = sorted((len(p) for p in containing), reverse=True)
    print(f'\n{len(containing)} undirected 2-plexes contain all {len(DM_PLEX)} authors of Fig. 7(b); sizes: {sizes}')
    (c,) = find(containing, names, ['A. Inkeri Verkamo'], size=13)
    print('Fig. 7(c): undirected 2-plex on DM (Verkamo is cited by 11 members but cites none)\n' + describe(c, edges, names))
    first = ['A. Inkeri Verkamo', 'Ada Wai-Chee Fu']
    order = [ids[n] for n in first] + [v for v in sorted(c, key=lambda v: names[v]) if names[v] not in first]
    draw_plex(c, edges, names, os.path.join(OUTPUT, 'fig7c_undirected_plex.pdf'), highlight='A. Inkeri Verkamo', order=order)

    print(f'\nFigures written to {OUTPUT}')


if __name__ == '__main__':
    main()
