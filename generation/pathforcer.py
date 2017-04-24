"""Forcer code generator.

This script generates all transitions from each top-level four-loop
topology all the way down to a single propagator. If the triangle rule
can be applied in a given topology, it is applied, otherwise it is
assumed that through a different method each line can be removed.

We define the momentum flow as going from the lowest vertex number to
a larger one. The momentum label is defined by the index in the sorted
edge list.
"""

from __future__ import division, print_function  # for Python 2 support

import argparse
import collections
import pickle  # to save quickly save and load the graph, for testing
import re
import textwrap
from copy import deepcopy
from itertools import combinations, permutations, product

from igraph import Graph, IN

import numpy
from numpy.linalg import inv

from sympy import Add, Matrix, Mul, Symbol, ilcm, solve, symbols


# add support for the old format, where
# some dots are allowed to be expressed in terms of the second
# external momentum edge line
_old_format = False

# a cache for solved systems
_momentum_cache = {}
_spcv_cache = {}
_propagator_cache = {}
_dot_cache = {}
_basis_cache = {}
_lin_cache = {}


def createFixedTopologies(loops):
    """Create topologies that have a fixed name, dot product and insertion lines.
    """
    fixed = []

    if loops == 3:
        g = Graph(n=8)
        g['name'] = 't1star5'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 7)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[5]['nonint'] = 1
        fixed.append(g)

    if loops == 4:
        g = Graph(n=10)
        g['name'] = 'no6'
        g['convert'] = True
        g.add_edges([(0, 1), (7, 8), (6, 7), (2, 6), (1, 2), (1, 7), (3, 7), (3, 4), (4, 8), (4, 6), (2, 3), (8, 9)])
        g['dots'] = [(1, 2), (2, 4), (0, 2), (0, 3)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'lala'
        g['convert'] = True
        g.add_edges([(0, 3), (5, 8), (4, 5), (3, 4), (2, 3), (1, 2), (1, 8), (2, 4), (1, 4), (1, 5), (8, 9)])
        g['dots'] = [(0, 5), (0, 2), (1, 4), (1, 5), (2, 4)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'nono'
        g['convert'] = True
        g['mindot'] = 14  # one dot product must be flipped for forcer
        g.add_edges([(0, 1), (4, 6), (2, 4), (1, 2), (1, 3), (3, 5), (5, 6), (4, 5), (3, 4), (2, 5), (6, 9)])
        g['dots'] = [(2, 8), (6, 7), (0, 2), (1, 2), (7, 9)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'cross'
        g['convert'] = True
        g.add_edges([(0, 1), (3, 5), (1, 3), (4, 5), (1, 4), (2, 5), (1, 2), (2, 3), (2, 4), (5, 9)])
        g['dots'] = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (2, 3)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'bebe'
        g['convert'] = True
        g.add_edges([(0, 1), (4, 5), (1, 4), (3, 5), (1, 3), (1, 2), (2, 4), (2, 3), (3, 4), (5, 9)])
        g['dots'] = [(0, 2), (0, 4), (0, 6), (1, 2), (2, 6), (1, 4)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'bubu'
        g['convert'] = True
        g.add_edges([(0, 4), (4, 5), (2, 4), (3, 4), (5, 6), (2, 5), (3, 5), (1, 6), (1, 2), (1, 3), (6, 9)])
        g['dots'] = [(0, 2), (0, 8), (2, 3), (2, 9), (3, 8)]
        g.es['nonint'] = 0
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'lastar5'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 6), (5, 6), (4, 5), (3, 5), (2, 6), (4, 9)])  # Q, p1, p2, p3, p5, p8: same direction as mincer
        g['dots'] = [(1, 3)]
        g.es['nonint'] = 0
        g.es[5]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'nostar0'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 6), (5, 6), (4, 5), (3, 6), (2, 5), (4, 9)])  # Q, p1, p2, p3, p5, p7: same direction as mincer
        g['dots'] = [(0, 2)]
        g.es['nonint'] = 0
        g.es[0]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'nostar5'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 6), (5, 6), (4, 5), (3, 6), (2, 5), (4, 9)])  # Q, p1, p2, p3, p5, p7: same direction as mincer
        g['dots'] = [(0, 2)]
        g.es['nonint'] = 0
        g.es[5]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'nostar6'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 6), (5, 6), (4, 5), (3, 6), (2, 5), (4, 9)])  # Q, p1, p2, p3, p5, p7: same direction as mincer
        g['dots'] = [(0, 2)]
        g.es['nonint'] = 0
        g.es[6]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'bustar5'
        g.add_edges([(0, 2), (3, 5), (1, 5), (4, 5), (1, 3), (1, 4), (2, 3), (2, 4), (5, 9)])  # Q, p1, p4, p5, p6: same direction as mincer
        g['dots'] = [(0, 4), (0, 5)]
        g.es['nonint'] = 0
        g.es[5]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'fastar2'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 5), (4, 5), (3, 5), (2, 5), (4, 9)])  # Q, p1, p2, p3, p5, p6: same direction as mincer
        g['dots'] = [(1, 3), (0, 2)]
        g.es['nonint'] = 0
        g.es[2]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 'fastar3'
        g.add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 5), (4, 5), (3, 5), (2, 5), (4, 9)])  # Q, p1, p2, p3, p5, p6: same direction as mincer
        g['dots'] = [(1, 3), (0, 2)]
        g.es['nonint'] = 0
        g.es[1]['nonint'] = 1  # in forcer.h n1 is non-integer instead of 3, check if correct
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 't1star05'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 9)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[0]['nonint'] = 1
        g.es[5]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 't1star24'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 9)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[2]['nonint'] = 1
        g.es[4]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 't1star34'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 9)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[3]['nonint'] = 1
        g.es[4]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 't1star45'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 9)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[4]['nonint'] = 1
        g.es[5]['nonint'] = 1
        fixed.append(g)

        g = Graph(n=10)
        g['name'] = 't1star55'
        g.add_edges([(0, 1), (2, 3), (1, 2), (1, 4), (3, 4), (2, 4), (3, 9)])  # Q, p1, p2, p4: same direction as mincer
        g['dots'] = []
        g.es['nonint'] = 0
        g.es[5]['nonint'] = 2
        fixed.append(g)

    for d in fixed:
        d['carpet'] = []
        d['triangle'] = []
        d['diamond'] = []
        d['insertion'] = []
        d['auto'] = getAutomorphisms(d)

    return fixed


def createTopLevelDiagrams(loops):
    """Creates all top level diagrams.
    """
    diagrams = []
    if loops == 2:
        diagrams = [Graph(n=6) for _ in range(1)]
        diagrams[0]['name'] = 'T1'
        diagrams[0].add_edges([(0, 1), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5)])
        diagrams[0]['dots'] = []

    if loops == 3:
        diagrams = [Graph(n=8) for _ in range(3)]
        diagrams[0]['name'] = 'benz'
        diagrams[0].add_edges([(0, 1), (1, 2), (1, 5), (2, 3), (2, 4), (3, 4), (3, 6), (4, 5), (5, 6), (6, 7)])
        diagrams[0]['dots'] = [(3, 9)]

        diagrams[1]['name'] = 'la'
        diagrams[1].add_edges([(0, 1), (1, 2), (1, 4), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (5, 6), (6, 7)])
        diagrams[1]['dots'] = [(4, 5)]

        diagrams[2]['name'] = 'no'
        diagrams[2].add_edges([(0, 1), (3, 4), (2, 3), (1, 2), (1, 6), (5, 6), (4, 5), (3, 6), (2, 5), (4, 7)])
        diagrams[2]['dots'] = [(0, 2)]

    if loops == 4:
        diagrams = [Graph(n=10) for _ in range(11)]
        # note that the first and the last vertices are always the external vertices
        diagrams[0]['name'] = 'no2'
        diagrams[0]['convert'] = True
        diagrams[0].add_edges([(0, 1), (5, 8), (4, 5), (2, 4), (1, 2), (1, 6), (3, 6), (3, 7), (7, 8), (5, 6), (3, 4), (2, 7), (8, 9)])
        diagrams[0]['dots'] = [(0, 2), (1, 4), (0, 3)]

        diagrams[1]['name'] = 'no1'
        diagrams[1]['convert'] = True
        diagrams[1].add_edges([(0, 1), (4, 8), (3, 4), (2, 3), (1, 2), (1, 5), (5, 6), (6, 7), (7, 8), (3, 7), (2, 6), (4, 5), (8, 9)])
        diagrams[1]['dots'] = [(2, 4), (0, 2), (0, 3)]

        diagrams[2]['name'] = 'haha'
        diagrams[2]['convert'] = True
        diagrams[2].add_edges([(0, 3), (5, 8), (4, 5), (3, 4), (2, 3), (1, 2), (1, 8), (5, 6), (1, 6), (4, 7), (2, 7), (6, 7), (8, 9)])
        diagrams[2]['dots'] = [(0, 2), (0, 5), (1, 4)]

        # topology 27, diamond
        diagrams[9]['name'] = 'no3'
        diagrams[9].add_edges([(0, 1), (4, 8), (7, 8), (1, 4), (1, 5), (3, 4), (2, 7), (6, 7), (2, 5), (5, 6), (2, 3), (3, 6), (8, 9)])
        diagrams[9]['dots'] = [(0, 6), (0, 8), (3, 6)]

        diagrams[3]['name'] = 'la4'
        diagrams[3].add_edges([(0, 1), (1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (3, 6), (4, 7), (4, 8), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[3]['dots'] = [(3, 7), (4, 5), (4, 7)]

        diagrams[4]['name'] = 'j1'
        diagrams[4].add_edges([(0, 1), (1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (3, 8), (4, 5), (4, 6), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[4]['dots'] = [(0, 3), (0, 7), (1, 7)]

        diagrams[5]['name'] = 'j6'
        diagrams[5].add_edges([(0, 1), (1, 2), (1, 6), (2, 3), (2, 4), (3, 4), (3, 8), (4, 5), (5, 6), (5, 7), (6, 7), (7, 8), (8, 9)])
        diagrams[5]['dots'] = [(0, 3), (0, 8), (3, 8)]

        diagrams[6]['name'] = 'j2'
        diagrams[6].add_edges([(0, 1), (1, 2), (1, 7), (2, 3), (2, 6), (3, 4), (3, 5), (4, 5), (4, 8), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[6]['dots'] = [(0, 3), (0, 5), (1, 5)]

        diagrams[7]['name'] = 'j5'
        diagrams[7].add_edges([(0, 1), (1, 2), (1, 6), (2, 3), (2, 5), (3, 4), (3, 5), (4, 7), (4, 8), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[7]['dots'] = [(0, 3), (1, 7), (3, 7)]

        diagrams[8]['name'] = 'j11'
        diagrams[8].add_edges([(0, 1), (1, 2), (1, 5), (2, 3), (2, 5), (3, 4), (3, 7), (4, 6), (4, 8), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[8]['dots'] = [(0, 5), (4, 5), (4, 7)]

        diagrams[10]['name'] = 'j4'
        diagrams[10].add_edges([(0, 1), (1, 2), (1, 6), (2, 3), (2, 7), (3, 4), (3, 5), (4, 5), (4, 8), (5, 6), (6, 7), (7, 8), (8, 9)])
        diagrams[10]['dots'] = [(0, 3), (0, 5), (1, 5)]

    if loops == 5:
        diagrams = [Graph(n=12) for _ in range(64)]
        diagrams[0].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 9), (4, 10), (5, 7), (5, 9), (6, 8), (6, 10), (7, 10), (8, 9), (2, 11)])
        diagrams[1].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 9), (4, 10), (5, 7), (5, 9), (6, 8), (6, 10), (7, 8), (9, 10), (2, 11)])
        diagrams[2].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 9), (4, 10), (5, 7), (5, 8), (6, 9), (6, 10), (7, 9), (8, 10), (2, 11)])
        diagrams[3].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 9), (5, 8), (5, 10), (6, 9), (6, 10), (7, 10), (8, 9), (2, 11)])
        diagrams[4].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 9), (5, 8), (5, 10), (6, 9), (6, 10), (7, 8), (9, 10), (2, 11)])
        diagrams[5].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 9), (5, 7), (5, 10), (6, 8), (6, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[6].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 9), (5, 7), (5, 10), (6, 8), (6, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[7].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 9), (5, 7), (5, 8), (6, 9), (6, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[8].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 8), (5, 9), (5, 10), (6, 9), (6, 10), (7, 9), (8, 10), (2, 11)])
        diagrams[9].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 7), (4, 8), (5, 7), (5, 9), (6, 9), (6, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[10].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 10), (6, 8), (6, 9), (7, 8), (7, 10), (9, 10), (2, 11)])
        diagrams[11].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 9), (6, 10), (7, 9), (7, 10), (8, 10), (2, 11)])
        diagrams[12].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 9), (6, 10), (7, 8), (7, 10), (9, 10), (2, 11)])
        diagrams[13].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 7), (6, 10), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[14].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 7), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[15].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 7), (6, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[16].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 7), (6, 8), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[17].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 7), (6, 8), (6, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[18].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 9), (6, 10), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[19].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 9), (6, 10), (7, 8), (8, 10), (9, 10), (2, 11)])
        diagrams[20].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 8), (6, 10), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[21].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 8), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[22].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 8), (6, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[23].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 7), (6, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[24].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 8), (6, 9), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[25].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 9), (6, 10), (7, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[26].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 9), (6, 10), (7, 8), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[27].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[28].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 8), (6, 9), (7, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[29].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 7), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[30].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 7), (5, 8), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[31].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 7), (4, 8), (5, 9), (6, 10), (7, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[32].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 7), (4, 8), (5, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[33].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 7), (4, 8), (5, 7), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[34].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 9), (5, 10), (6, 9), (6, 10), (7, 8), (7, 9), (8, 10), (2, 11)])
        diagrams[35].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 9), (6, 10), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[36].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 9), (6, 10), (7, 8), (8, 10), (9, 10), (2, 11)])
        diagrams[37].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 10), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[38].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[39].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[40].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 9), (6, 7), (6, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[41].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 7), (5, 8), (6, 9), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[42].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (5, 9), (6, 10), (7, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[43].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (5, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[44].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (5, 7), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[45].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 5), (4, 7), (5, 8), (6, 9), (6, 10), (7, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[46].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 5), (4, 7), (5, 8), (6, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[47].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 5), (4, 7), (5, 8), (6, 8), (6, 9), (7, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[48].add_edges([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 5), (4, 7), (5, 8), (6, 7), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[49].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (5, 10), (6, 7), (6, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[50].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (5, 10), (6, 7), (6, 9), (7, 8), (8, 10), (9, 10), (2, 11)])
        diagrams[51].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 9), (5, 10), (6, 7), (6, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[52].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 10), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[53].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[54].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[55].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 9), (6, 7), (6, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[56].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 9), (6, 10), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[57].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 8), (5, 9), (6, 10), (7, 8), (7, 10), (8, 9), (9, 10), (2, 11)])
        diagrams[58].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 8), (5, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[59].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 8), (5, 9), (6, 8), (7, 9), (7, 10), (8, 10), (9, 10), (2, 11)])
        diagrams[60].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 7), (5, 8), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[61].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 9), (6, 10), (7, 9), (7, 10), (8, 9), (8, 10), (2, 11)])
        diagrams[62].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 9), (6, 10), (7, 8), (7, 9), (8, 10), (9, 10), (2, 11)])
        diagrams[63].add_edges([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (6, 9), (7, 10), (8, 9), (8, 10), (9, 10), (2, 11)])

    for i, d in enumerate(diagrams):
        if 'name' not in d.attributes():
            d['name'] = 'j' + str(i)

        if 'dots' not in d.attributes():
            d['dots'] = []

        d.es['nonint'] = 0  # all edges are integers for now

        findTopologyActionAndDots(d, d['dots'] == [])  # should be before auto
        d['auto'] = getAutomorphisms(d)

    return diagrams


def mapAfterRemoval(numinds, reminds):
    """Create a mapping of indices with removal of elements.

    This is a helper function that maps indices to new indices after some
    elements have been removed. The indices are 0-based.

    Parameters
    ----------
    numinds : int
        The number of indices.
    reminds : sequence of int
        The indices to be removed.

    Returns
    -------
    list of int
        A list that maps the old indices to the new indices after the specified
        elements have been removed. Removed indices are mapped to `numinds`.

    Examples
    --------
    >>> mapAfterRemoval(6, [1, 3])
    [0, 6, 1, 6, 2, 3]

    The result represents a mapping from the old indices (0, 1, 2, 3, 4, 5) to
    the new indices (0, 2, 4, 5), which are generated by removing the indices
    1 and 3, as
      0 -> 0
      1 -> 6 (removed)
      2 -> 1
      3 -> 6 (removed)
      4 -> 2
      5 -> 3

    """
    a = [i for i in range(numinds) if i not in reminds]
    return [a.index(i) if i in a else numinds for i in range(numinds)]


def solvePropagatorSystem(g, basis=None):
    """Solve the given system of propagators for dependent momenta.

    Solve the system of propagators in the given topology for dependent momenta
    and express them in terms of independent momenta.

    Parameters
    ----------
    g : igraph.Graph
        The topology.
    basis : sequence of int, optional
        The indices for the independent loop momenta. If not given, they will
        be automatically chosen via the spanning tree.

    Returns
    -------
    pp : list of sympy.Symbol
        The momenta for all edges.
    sps : list of sympy.Expr
        The independent set of the scalar products.
    dd_reps : dict of sympy.Symbol -> sympy.Expr
        The solutions for the dependent momenta.
    dpind : list of (int, int)
        The list of pairs of edge IDs.
    dp : list of sympy.Symbol
        The list of symbols for scalar products.
    dpmap : dict of (int, int) -> sympy.Symbol
        The mapping from an element in `dpind` to the corresponding element
        in `dp`.
    cachekey : immutable
        The cache key.
    """
    # comparison by edge list of the graph
    cachekey = tuple(g.get_edgelist() + ([] if basis is None else list(basis)))
    if cachekey in _momentum_cache:
        return _momentum_cache[cachekey]

    # symbols for all edges, [p0, p1, ...]
    pp = symbols(['p' + str(i) for i in range(g.ecount())])

    if basis is None:
        st = g.spanning_tree(return_tree=False)
        ll = [i for i in range(g.ecount()) if i not in st]
    else:
        ll = [b for b in basis if b != 0]
    ees = [e.index for e in g.es if isExternal(g, e)]

    # dependent momenta as edge IDs, could contain external momentum
    dd = [i for i in range(g.ecount())
          if i not in ll and (i != ees[0] or basis is not None and
                              ees[0] not in basis)]

    # as symbols
    internals = [pp[i] for i in ll]
    externals = [pp[ees[0]]]
    dependents = [pp[i] for i in dd]

    # solve dependent momenta by momentum conservation at each vertex
    eqs = [pp[ees[1]] - pp[ees[0]]]
    for i in range(g.vcount()):
        if len(g.neighbors(i)) < 2:
            continue
        eq = 0
        for j in set(g.neighbors(i)):
            sign = 1 if j >= i else -1
            for e in g.es:
                if e.tuple == tuple(sorted((i, j))):
                    eq += sign * pp[e.index]
        eqs.append(eq)

    dd_solved = solve(eqs, dependents,
                      dict=True, simplify=False, manual=True)

    # scalar products to be solved as products of two independent momenta
    sps = [a * b for i, a in enumerate(internals) for b in internals[i:]]
    sps += [a * b for a in internals for b in externals]

    # map from a pair of edge IDs to a symbol for the dot product
    dpind = [(i, j) for i in range(g.ecount()) for j in range(i, g.ecount())]
    dp = symbols(['p' + str(i) + '.p' + str(j) for i, j in dpind])
    dpmap = dict(zip(dpind, dp))

    # cache the result
    result = (pp, sps, dd_solved[0], dpind, dp, dpmap, cachekey)
    _momentum_cache[cachekey] = result
    return result


def SPCoeffVector(expr, reps, sps, cachekey):
    """Decompose the given expression in a vector space spanned by `sps`.

    Decompose the given expression as a linear combination of the scalar
    products and return the coefficients and the inhomogeneous term as a list.

    Parameters
    ----------
    expr : sympy.Expr
        The expression to be decomposed.
    reps : iterable
        The replacements substituted before the decomposition.
    sps : list of sympy.Expr
        The independent set of the scalar products spanning the basis.
    cachekey : immutable
        The cache key corresponding to the operation used for getting `reps`
        and `sps`.

    Returns
    -------
    list of sympy.Expr
        The coefficients and the inhomogeneous term. The latter is in the last
        element.
    """
    # After performing E := expr.subs(reps).expand(), decompose it as
    # E = dot(A, sps) + C, then return A + [C]. The result is cached with
    # `cachekey`.
    cachekey2 = (cachekey, expr)
    if cachekey2 in _spcv_cache:
        return _spcv_cache[cachekey2]

    expr = expr.subs(reps).expand()
    (constant, poly) = expr.as_independent(*sps, as_Add=True)
    result = [poly.coeff(sp) for sp in sps]
    result.append(constant)
    _spcv_cache[cachekey2] = result
    return result


def solveDotProductSystem(g, dots, disregard=None, basis=None):
    """Solve the given system of propagators for scalar products.

    Solve the momentum system for topology `g` with the dot products `dots` for
    scalar products involving loop momenta. Optionally, a propagator can be
    treated as a dependent variable using `disregard`.

    Parameters
    ----------
    g : igraph.Graph
        The topology.
    dots : sequence of (int, int)
        The dot products for irreducible numerators.
    disregard : sequence of int, optional
        The edge IDs corresponding to `pj^2` to be treated as dependent
        variables.
    basis : sequence of int, optional, deprecated
        Currently not used.

    Returns
    -------
    dpmap : dict of (int, int) -> sympy.Symbol
        The mapping from a pair of edge IDs to the corresponding scalar
        product.
    sol : dict of sympy.Symbol -> sympy.Expr
        The solutions for all scalar products.
    """
    if disregard is None:
        disregard = []

    # special case for the final state without vertices
    if g.ecount() == 1:
        name = symbols(['p0'])[0]
        return {(0, 0): name}, {name: name}

    cachekey2 = (tuple(g.get_edgelist()), tuple(dots), tuple(disregard))
    if cachekey2 in _dot_cache:
        return _dot_cache[cachekey2]

    # computation depending only on `g'
    # do not use a custom basis here
    pp, sps, dd_reps, dpind, dp, dpmap, cachekey = solvePropagatorSystem(g)

    # the following stuff depends on `dots' and `disregard' as well as `g'

    # propagators, including dot products in the numerator
    pp_reps = [(pp[i] * pp[i], dpmap[(i, i)])
               for i in range(g.ecount())
               if i not in disregard and
               0 not in g.es[i].tuple and g.vcount() - 1 not in g.es[i].tuple]

    pp_reps += [(pp[e0] * pp[e1], dpmap[(e0, e1)]) for e0, e1 in dots]
    propagators = [x[0] for x in pp_reps]

    # check the number of equations vs. variables
    if len(propagators) < len(sps):
        raise AssertionError('Incomplete momentum system.')
    if len(propagators) > len(sps):
        raise AssertionError('Overdetermined momentum system.')

    # solve the system of propagators for the scalar products:
    #   Y_i = A_{ij} sp_j + C_i  => R := sp = A^{-1} (Y-C)

    # replace lists
    qq_reps = next((pp[i] * pp[i], dpmap[(i, i)])
                   for i in range(g.ecount()) if 0 in g.es[i].tuple)

    v = [SPCoeffVector(p, dd_reps, sps, cachekey) for p in propagators]
    A = numpy.asmatrix([[int(a) for a in x[:-1]] for x in v])
    C = [x[-1].subs(*qq_reps) for x in v]

    r = inv(A) * 2  # numerical inversion!
    r = r.astype(numpy.int).tolist()
    Ainv = numpy.array(Matrix(r) / 2)

    Y = numpy.array([x[1] for x in pp_reps])
    R = numpy.dot(Ainv, Y - C)

    # now all pi * pj can be reduced by matrix/vector operations

    # map from a symbol for the dot product to the solved expression
    vv = [SPCoeffVector(pp[i] * pp[j], dd_reps, sps, cachekey)
          for i, j in dpind]
    AA = [x[:-1] for x in vv]
    CC = [x[-1].subs(*qq_reps) for x in vv]

    # NOTE: somehow AA.(Ainv.(Y-C)) is faster than (AA.Ainv).(Y-C)
    dpr = numpy.dot(AA, R) + CC
    dprmap = dict(zip(dp, dpr))

    _dot_cache[cachekey2] = (dpmap, dprmap)
    return dpmap, dprmap


def findBestMomentumBasis(g, dots, olddots, disregard=[], include=set(), exclude=True):
    """Find the momentum system that results in the fewest terms being generated for
    a transition from olddots. If the include is None, the momenta in the dot products
    will be added to the basis. Otherwise, include is used
    """
    if g.ecount() == 1:
        return []

    # exclude propagators that will be removed
    if exclude:
        exclude = []
        if g['insertion']:
            exclude = [g['insertion'][0][1]]
        if g['triangle']:
            exclude = g['triangle'][:2]
        if g['diamond']:
            exclude = g['diamond'][0]

    # make sure the dot products are expressed in independent momenta
    # to prevent unnecessary rewriting
    if len(include) == 0:
        include = set(e for d in dots for e in d if not exclude or e not in exclude)

    basiskey = tuple(e.tuple for e in g.es)
    if basiskey not in _basis_cache:
        bases = []
        loops = g.ecount() - len(set(v for e in g.es for v in e.tuple)) + 1
        for c in combinations([e.index for e in g.es if g.vcount() - 1 not in e.tuple],
                    loops + 1):
            if linearIndependenceProp(g, c):
                bases.append(c)
        _basis_cache[basiskey] = bases
        print('Generated basis for', g['name'])
    else:
        bases = _basis_cache[basiskey]

    """
    # keep on removing until system is no longer reducible
    bases = [(0, [e.index for e in g.es if g.vcount() - 1 not in e.tuple])]
    while True:
        #print(bases)
        #print('loop')
        newbases = []
        for s, x in bases: # only remove if greater than the last one
            if x[-1] < s:
                continue
            newbase = [(r, [xx for xx in x if xx != r]) for r in x if r >= s]  # slow?
            newbase = [xx for xx in newbase if linearIndependenceProp(g, xx[1])]
            #print('nb', newbase)
            newbases += newbase
        if newbases == []:
            break
        bases = newbases

    print(bases)
    exit(1)
    """

    # filter bases
    fbases = [b for b in bases if include <= set(b) and
        (not exclude or len(set(b) & set(exclude)) == 0)]

    # if Q is not in the basis, remove another variable. Pathforcer supports bases without Q,
    # but apparently this behaviour is unwanted for momsubs.
    fbases = [b for b in fbases if 0 in b]

    if fbases == []:
        includes = [list(set(include) - {0})]
        while True:
            includes = [[j for j in inc if j != i] for inc in includes for i in inc]
            fbases = [b for b in bases if any((set(i) | {0}) <= set(b) for i in includes) and
                (not exclude or len(set(b) & set(exclude)) == 0)]
            print('Relaxing basis suggestion', include, 'to', fbases, 'for', g['name'])
            if fbases != []:
                break
    bases = fbases

    # find the basis that minimizes the number of terms in the rewrite
    if bases == []:
        # should not happen anymore
        # FIXME: dotlines is undefined.
        print('WARNING: no basis with', (include if include is not None else dotlines), 'for', g['name'])
    basis, bestscore = None, 0
    for ba in bases:
        varmap, sol = solveDotProductSystem(g, dots, disregard=disregard, basis=ba)
        dotmap = parseMap([(a, b) for a in range(g.ecount()) for b in range(a, g.ecount())], varmap, sol)
        dm = [(k, len(v)) for k, v in dotmap.items() if len(v) > 1 and k[0] in ba and k[1] in ba and (olddots is None or k in olddots)]
        score = sum(d[1]**2 for d in dm)
        if basis is None or score < bestscore:
            basis = ba
            bestscore = score
    return basis


def linearIndependenceProp(g, indep):
    linkey = tuple(e for e in g.es)
    if linkey not in _lin_cache:
        pp = symbols(['p' + str(i) for i in range(g.ecount())])
        ees = [e.index for e in g.es if isExternal(g, e)]

        eqs = [pp[ees[1]] - pp[ees[0]]]
        for i in range(g.vcount()):
            if len(g.neighbors(i)) < 2:
                continue
            eq = 0
            for j in set(g.neighbors(i)):
                sign = 1 if j >= i else -1
                for e in g.es:
                    if e.tuple == tuple(sorted((i, j))):
                        eq += sign * pp[e.index]
            eqs.append(eq)
        _lin_cache[linkey] = (pp, eqs)
    else:
        pp, eqs = _lin_cache[linkey]

    dd = [pp[i] for i in range(g.ecount()) if i not in indep]

    m = []
    for eq in eqs:
        (_, poly) = eq.as_independent(*dd, as_Add=True)
        m.append([poly.coeff(sp) for sp in dd])

    ma = numpy.matrix([[float(a) for a in mm] for mm in m])
    ddd = numpy.dot(ma.T, ma)
    ans = numpy.linalg.det(ddd)

    return abs(ans) > 0.5


def linearIndependence(g, dots, disregard=None):
    """Return True for a linearly independent system.

    Return True if propagators in the topology and the dot products are all
    linearly independent. Optionally, propagators can be excluded in the check
    by specifying the corresponding edge IDs in `disregard`.

    Parameters
    ----------
    g : igraph.Graph
        The topology.
    dots : sequence of (int, int)
        The dot products for irreducible numerators.
    disregard : sequence of int, optional
        The edge IDs corresponding to `pj^2` to be excluded in the check.

    Returns
    -------
    bool
        True if the system is linearly independent.
    """
    if disregard is None:
        disregard = []

    # computation depending only on `g'
    pp, sps, dd_reps, dpind, dp, dpmap, cachekey = solvePropagatorSystem(g)

    # the following stuff depends on `dots' and `disregard' as well as `g'

    # propagators, including dot products in the numerator
    propagators = [
        pp[i] * pp[i] for i in range(g.ecount())
        if (i not in disregard and 0 not in g.es[i].tuple and
            g.vcount() - 1 not in g.es[i].tuple)
    ]
    for e0, e1 in dots:
        propagators.append(pp[e0] * pp[e1])

    # check the number of equations vs. variables
    if len(propagators) > len(sps):
        return False

    # each propagator is translated into a vector consisting of the
    # coefficients of sps
    for i in range(len(propagators)):
        cachekey2 = (cachekey, propagators[i])
        if cachekey2 in _propagator_cache:
            propagators[i] = _propagator_cache[cachekey2]
        else:
            p = SPCoeffVector(propagators[i], dd_reps, sps, cachekey)[:-1]
            # scale the vector such as it becomes an integer vector
            lcm = 1
            for e in p:
                if not e.is_integer and e.is_rational:
                    lcm = ilcm(lcm, e.q)
            if lcm != 1:
                p = [lcm * e for e in p]
            # expect an integer vector (no other objects, e.g., symbols)
            for e in p:
                if not e.is_integer:
                    raise AssertionError('Non-integer vector found: ', p)
            # converting SymPy integer to float makes numpy.dot faster
            p = [float(e) for e in p]
            propagators[i] = _propagator_cache[cachekey2] = p

    # compute the Gram determinant
    m = numpy.matrix(
        [[numpy.dot(i, j) for j in propagators] for i in propagators])
    det = numpy.linalg.det(m)
    # Note that the determinant must be an integer because in the above we
    # ensured the matrix is an integer matrix. Even if the value is computed as
    # a floating-number, we can safely say wheather it is zero or non-zero.
    return abs(det) > 0.5


def solveDotAssignment(relabeling):
    """Finds out which old dot product maps to which new dot product, given relabeling `relabeling'. This is useful
    to reuse indices. It uses backtracking if a solution is difficult to find.
    """
    newlabels = {v for vs in relabeling.values() for v in vs}
    if len(newlabels) != len(relabeling):
        raise AssertionError('Cannot solve system: ', relabeling, newlabels)

    def backtrack(sol, labeling):
        if len(sol) == len(newlabels):
            return sol

        # try the element with the minimal number of elements
        minitem = min(labeling.items(), key=lambda x: len(x[1]))

        for p in minitem[1]:
            newsol = deepcopy(sol)
            newlabeling = deepcopy(labeling)
            newsol[minitem[0]] = p
            del newlabeling[minitem[0]]

            for k, v in newlabeling.items():
                if p in v:
                    v.remove(p)

            r = backtrack(newsol, newlabeling)
            if r is not None:
                return r
        return None

    res = backtrack({}, relabeling)
    if res is None:
        raise AssertionError('Could not solve dot assignment')
    return res


def parseMap(olddots, varnames, dotsystem):
    """Parses the result of a solved dot product system. It writes the fallout of the rewrite,
    namely the edge indices of all the propagators that need to be raised, in a convenient structure.
    """
    nametoedge = dict((v, k) for k, v in varnames.items())

    # convert the option to indices, parse the expression
    dotfull = {}
    for d in olddots:  # re-express the old vars
        var = varnames[d]
        if var not in dotsystem:
            raise AssertionError('Variable', var, 'not found in dot product relations')

        a = dotsystem[var]

        # small fix if there is just a single term
        terms = a.args if a.func == Add else [a]
        dotmap = []
        for x in terms:
            v = x.as_coeff_mul()
            t = nametoedge[v[1][0]]
            dotmap.append((v[0], t))
        dotfull[d] = dotmap

    return dotfull


def remapdots(olddots, varnames, dotsystem):
    """Maps old dots to new dots, given a solved system of dot product relations.
    """
    dotfull = parseMap(olddots, varnames, dotsystem)
    posdots = {k: [t[1] for t in d if t[1][0] != t[1][1]] for k, d in dotfull.items()}

    dotassign = solveDotAssignment(posdots)
    for k, x in dotassign.items():
        dotfull[k] = (x, dotfull[k])

    return dotfull


def findDots(g, edge_ids, dotsug=[], keep_edge=False, max_edge=None):
    """Re-expresses a topology so that it does not depend on a given propagators anymore:
    the current dot products are rewritten, and optionally prop^2 as well. It also maximizes the
    occurrence of max_edge in the dots.
    """
    print('*** Re-express', [e.tuple for e in g.es], 'with dots',
          g['dots'], 'without edges', edge_ids)

    if g.ecount() == 1:
        return []

    # find a linear independent system
    edge_ids = edge_ids + [next(e.index for e in g.es if g.vcount() - 1 == e.target)]
    all_dots = [(i, j) for i in range(g.ecount()) for j in range(i + 1, g.ecount())
                if i not in edge_ids and j not in edge_ids]

    # maximize the number of times an edge_occur is chosen
    if max_edge is not None:
        all_dots = list(sorted(all_dots, key=lambda x: 0 if max_edge in x else 1))

    # make sure the old dots are placed second, after the max_edge ones
    #all_dots = list(sorted(all_dots, key=lambda x: 3 - 2 * (max_edge in x) - (x in dotsug)))

    l = g.ecount() - len({v for e in g.get_edgelist() for v in e}) + 1
    dotsleft = (l + 1) * l // 2 + l - g.ecount() + 2 + (0 if keep_edge else 1)

    dots = []
    mincoeff = 0

    for _ in range(dotsleft):
        for i, d in enumerate(all_dots[mincoeff:]):
            if linearIndependence(g, dots + [d], [] if keep_edge else [edge_ids[0]]):
                dots.append(d)
                mincoeff += i + 1
                break

    print('Found dots', dots)
    return dots


def fixDotsAfterLineRemoval(g, dotmap, edge_id):
    """
    Some dot products in g have a new index in h due to the removal of a line. Fix this in the dot product mapping.
    Additionally, a sign flip could occur due to the renaming of the vertex.
    """
    emap = [i if i < edge_id else i - 1 for i in range(g.ecount())]
    edge = g.es[edge_id].tuple

    def get_sign(t):
        r = (t[0] if t[0] != edge[1] else edge[0], t[1] if t[1] != edge[1] else edge[0])
        return 1 if r[0] <= r[1] else -1

    def replace_dot(coeff, d):
        return coeff * get_sign(g.es[d[0]].tuple) * get_sign(g.es[d[1]].tuple), (emap[d[0]], emap[d[1]])

    newmap = {}
    for k, v in dotmap.items():
        newmap[k] = [replace_dot(1, v[0])[1], [replace_dot(*x) for x in v[1]]]
    return newmap


def removeEdge(g, edge_id, preserveDirection=False):
    """Removes an edge and relabels the highest vertex of the removed edge to the lowest. Leaves the dot products as is and
    leaves the order of the edges as is. This method is not suitable for insertions.
    """
    h = g.copy()
    h.delete_edges(None)  # clear
    for attr in ('transgraphs', 'convert', 'mindot'):
        if attr in h.attributes():
            del h[attr]
    h['removed_edge'] = edge_id
    edge = g.es[edge_id].tuple

    # rewrite the highest vertex to the lowest
    e = [(a.source if a.source != edge[1] else edge[0], a.target if a.target != edge[1] else edge[0]) for a in g.es]

    for ee, i in zip(e, g.es):
        if ee[0] != ee[1]:  # skip the removed line
            if preserveDirection:
                h.add_edge(*ee, nonint=i['nonint'])
            else:
                h.add_edge(*sorted(ee), nonint=i['nonint'])

    return h


def findTadpole(g, v):
    """Finds out if the diagram has a tadpole at vertex v.
    """
    # does the tadpole have no vertices?
    if any(e[0] == e[1] for e in g.get_edgelist()):
        return True

    def dfs(c):
        visited.add(c)
        if g.degree(c) == 1:
            return True

        for n in g.neighbors(c):
            if n not in visited:
                if dfs(n):
                    return True

        return False

    # starting at v and checking every neighbour: if we cannot find a
    # way to an external line without passing v again, we have found a tadpole
    for n in set(g.neighbors(v)):
        visited = {v}
        if not dfs(n):
            return True
    return False


def findDiamond(g):
    """Finds a diamond in a graph. It may not work completely if there are
    insertions from the top vertex to a central vertex.
    """
    # start at bottom vertex
    # diamond should have at least 2 lines with a 3 point vertex, else it is just insertions
    for b in range(g.vcount()):
        if len(g.neighbors(b)) < 3 or len([len(set(g.neighbors(x))) == 3 for x in g.neighbors(b)]) < 2 \
            or 0 in g.neighbors(b) or g.vcount() - 1 in g.neighbors(b):
            continue

        # second neighbours with itself included (latter is necessary for insertions)
        cc = [set(g.neighbors(x)) | {x} for x in set(g.neighbors(b))]
        topv = cc[0]
        for x in cc[1:]:
            topv &= x
        res = list(topv - {b})  # top vertex options

        for t in res:
            # all the neighbours of this vertex should have 3 neighbours, none of which connected to external vertices
            # or be connected to the top vertex
            if all(len(g.neighbors(x)) == 3 and 0 not in g.neighbors(x) and g.vcount() - 1 not in g.neighbors(x) and
                    (len(set(g.neighbors(x))) == 2 or g.es[g.get_eid(x, list(set(g.neighbors(x)) - {b} - {t})[0])]['nonint'] == 0) and
                    g.es[g.get_eid(b, x)]['nonint'] == 0 for x in set(g.neighbors(b)) if x != t):

                # now collect all the components, first from the top to the center lines
                a = [g.get_eid(n, t) for n in g.neighbors(b) if n != t]
                bo = [g.get_eid(n, b) for n in g.neighbors(b) if n != t]
                c = [g.get_eid(n, k) for n in g.neighbors(b) for k in set(g.neighbors(n)) - {b, t}]
                s = [e.index for e in g.es if e.tuple == tuple(sorted([b, t]))]
                return a, bo, c, s
    return []


def mergePropagators(g, dotmap):
    """Merges propagators (edges) that have the same momentum.
    """
    pp, sps, dd_reps, dpind, dp, dpmap, cachekey = solvePropagatorSystem(g)

    # TODO: make sure that one of the outer lines is the representative of Q!
    ddmap = {}
    ddd = {x: i for i, x in enumerate(pp)}
    for op, nps in dd_reps.items():
        if 0 in g.es[pp.index(op)].tuple or g.vcount() - 1 in g.es[pp.index(op)].tuple:
            continue

        if nps in ddd:
            ddmap[pp.index(op)] = (1, ddd[nps])
        elif -nps in ddd:
            ddmap[pp.index(op)] = (-1, ddd[-nps])
        else:
            ddd[nps] = pp.index(op)

    # merge the outer lines if they are the only ones left
    if g.ecount() - len(ddmap) == 2 and len(set(ddmap.values())) == 1:
        ddmap[next(e.index for e in g.es if g.vcount() - 1 in e.tuple)] = list(ddmap.values())[0]

    # replace all dependent propagators with the independent ones
    def replace_dot(coeff, d):
        sign = (ddmap[d[0]][0] if d[0] in ddmap else 1) * (ddmap[d[1]][0] if d[1] in ddmap else 1)
        newd = tuple(sorted(
            (ddmap[d[0]][1] if d[0] in ddmap else d[0], ddmap[d[1]][1] if d[1] in ddmap else d[1])))
        return coeff * sign, newd

    dotmap = {k: (replace_dot(1, v[0])[1], [replace_dot(*x) for x in v[1]]) for k, v in dotmap.items()}

    h = g.copy()

    # copy the noninteger part
    for op, nps in ddmap.items():
        h.es[nps[1]]['nonint'] += h.es[op]['nonint']

    # now fix the dots for all the vertex merges
    for op in sorted(ddmap, reverse=True):
        dotmap = fixDotsAfterLineRemoval(h, dotmap, op)
        h = removeEdge(h, op)

    # special case: when we have only a single propagator left, the upper vertex should be g.vcount() - 1
    if h.ecount() == 1:
        h.add_edges([(0, h.vcount() - 1)])
        h.es[1]['nonint'] = h.es[0]['nonint']
        h.delete_edges([0])

    ddmap = [ddmap[e][1] if e in ddmap else e for e in range(g.ecount())]
    return h, ddmap, dotmap


def removeInsertions(g):
    """Removes one pair of the first repeated edge (insertion) from a topology. It also finds which dot products
    remain after integration. This method is also suited for carpets.
    """
    h = g.copy()
    h['did_insertion'] = True

    idsind, insdots = g['insertion']
    h.es[idsind[0]]['nonint'] += h.es[idsind[1]]['nonint'] + 1

    # all the dot products that are not involved in the insertion
    # should not have the insertion line at this point
    len_indep = len([d for d in g['dots'] if idsind[0] not in d])

    # map the dot products after the line has been removed
    dotmap = {d: (d, [(1, d)]) for d in g['dots']}  # identity mapping
    dotmap = fixDotsAfterLineRemoval(h, dotmap, idsind[1])
    sortednewdots = [dotmap[d][0] for d in g['dots']]

    # generate all extra dot products that are created from the insertion
    extradot = [d[0] if d[1] == idsind[0] else d[1] for d in sortednewdots if idsind[0] in d]
    extradot = [tuple(sorted(d)) for d in combinations(extradot, 2)]
    #print('Extra dots: ', extradot)

    h.delete_edges(idsind[1:])

    #print('New topology {0} with dots {1}, removed {2}'.format([e.tuple for e in h.es], h['dots'], remdots))

    # save the state after the insertion is resolved
    transgraph = h.copy()
    transgraph['dots'] = sortednewdots + extradot
    transgraph['extradots'] = sortednewdots[len_indep:]
    transgraph['namesuffix'] = 'b'
    transgraph['dotmap'] = dotmap.copy()
    transgraph['mapping'] = mapAfterRemoval(g.ecount(), idsind[1:])
    g['transgraphs'] = [transgraph]

    """
    The insertion has now been resolved. Now check if the topology has to change.
    """

    dotmap = {d: (d, [(1, d)]) for d in transgraph['dots']}  # identity mapping

    # note: merging lines could create double dots
    h, ddmap, dotmap = mergePropagators(h, dotmap)
    g['merge'] = ddmap

    transgraph = h.copy()
    transgraph['namesuffix'] = 'c'
    transgraph['dotmap'] = dotmap.copy()
    transgraph['dots'] = [dotmap[x][0] for x in g['transgraphs'][0]['dots']]
    #transgraph['merge'] = g['merge']
    g['transgraphs'].append(transgraph)

    print('Merged edges to {0} with dots {1}'.format([e.tuple for e in h.es], transgraph['dots']))

    # After the merge, we create a new dotmap to remove all dependent dots
    h['dots'] = transgraph['dots'][:len_indep]  # independent dots # FIXME: is this right?
    dotmap = {d: (d, [(1, d)]) for d in transgraph['dots']}

    print('Independent dots', h['dots'])

    # map the insertion dots to squares
    if any(d not in h['dots'] for d in dotmap):
        varnames, solved = solveDotProductSystem(h, h['dots'])
        dotrewrite = parseMap([v[0] for d, v in dotmap.items() if d not in h['dots']], varnames, solved)

        for k, v in dotrewrite.items():
            dotmap[k] = (k, v)

    g['transgraphs'][-1]['dotrewrite'] = dotmap

    l = h.ecount() - len([v for v in range(h.vcount()) if len(h.neighbors(v)) > 0]) + 1
    dotsleft = (l + 1) * l // 2 + l - h.ecount() + 2

    if l > 0 and len(h['dots']) != dotsleft:
        raise AssertionError('{1}/{0} dots are left: {2}'.format(dotsleft, len(h['dots']), h['dots']))

    print('Final dots:', h['dots'])

    return h


def mapEdge(c, h, mapping, mapExt=True):
    """Maps the edges of c to the edges of h, where c and h are isomorphic. It matches edges
    with similar epsilons. Insertions are permitted. The external lines can be matched to the
    lowest one using mapExt.
    """
    def getni(g, e):
        return 0 if (g.vs[e.source].degree() == 1 or
                     g.vs[e.target].degree() == 1) else e['nonint']

    # for duplicated edges, we keep a tracker of what edge it actually is
    edgecount = collections.defaultdict(int)
    edgemap = [0 for _ in range(c.ecount())]
    for e in c.es:
        edgemap[e.index] = edgecount[(e.tuple, getni(c, e))]
        edgecount[(e.tuple, getni(c, e))] += 1

    hmap = collections.defaultdict(list)
    for e in h.es:
        hmap[(e.tuple, getni(h, e))].append(e.index)

    ctoh = [tuple(sorted((mapping[e.tuple[0]], mapping[e.tuple[1]]))) for e in c.es]
    if mapExt:
        ctoh[next(e.index for e in c.es if c.vcount() - 1 == e.target)] = next(e.tuple for e in h.es if 0 == e.source)
        ctoh[next(e.index for e in c.es if 0 == e.source)] = next(e.tuple for e in h.es if 0 == e.source)
    return [hmap[(n, getni(c, e))][edgemap[e.index]] for n, e in zip(ctoh, c.es)]


def isExternal(g, edge):
    """Return True for external edges.

    Parameters
    ----------
    g : igraph.Graph
        The topology.
    edge : igraph.Edge
        The edge to be checked.

    Returns
    -------
    bool
        True if the given edge is one of external edges of the topology.
    """
    return g.vs[edge.source].degree() == 1 or g.vs[edge.target].degree() == 1


def areEdgesIsomorphic(g1, g2, e1, e2):
    """Define when all the edges between two vertices could be isomorphic.

    Note: this function is only called for one e1 and e2 for repeated edges
    """
    # collect nonints but ignore epsilons on outer lines
    ae1 = [e['nonint'] if not isExternal(g1, e) else 0
           for e in g1.es if e.tuple == g1.es[e1].tuple]  # obtain all edges
    ae2 = [e['nonint'] if not isExternal(g2, e) else 0
           for e in g2.es if e.tuple == g2.es[e2].tuple]

    # see if the distribution of epsilons is identical
    return collections.Counter(ae1) == collections.Counter(ae2)


def getAutomorphisms(graph):
    """Finds all automorphisms, including automorphisms of double edges. For each degenerate edge,
    all permutations are generated, keeping the attributes in place. Filters dot products that changes
    the set of dot products and generates automorphisms for cases where some dot products are 0.
    """
    if graph.ecount() == 1:
        return {}

    edgelist = [(e.tuple[0], e.tuple[1], {'nonint': e['nonint']
            if not isExternal(graph, e) else 0}) for e in graph.es]

    # for each degenerate edge, generate all permutations, keeping the attributes in place.
    # group equal edges, and make an index mapping for them
    groups = collections.defaultdict(list)
    edgecounter = collections.defaultdict(int)
    edgeindex = []
    for i, (a, b, j) in enumerate(edgelist):
        groups[(a, b, j['nonint'])].append(i)
        edgeindex.append(edgecounter[(a, b, j['nonint'])])
        edgecounter[(a, b, j['nonint'])] += 1

    # now map edge ids to groups
    grouplist = list(groups.values())
    groupmap = {}
    for i, v in enumerate(grouplist):
        for vv in v:
            groupmap[vv] = i

    # now generate all possible group orderings
    edgemappings = []
    for p in product(*[list(permutations(v)) for v in grouplist]):
        mapping = []
        for i in range(graph.ecount()):
            mapping.append(p[groupmap[i]][edgeindex[i]])
        edgemappings.append(mapping)

    # give each vertex that has no edges a unique color in order to prevent
    # fake isomorphisms
    vcolor = [0 if len(graph.neighbors(v)) > 0 else v + 1 for v in range(graph.vcount())]

    # for each vertex isomorphism, add each possible double edge mapping
    mappings = collections.defaultdict(list)
    for iso in graph.get_isomorphisms_vf2(graph, color1=vcolor, color2=vcolor, edge_compat_fn=areEdgesIsomorphic):
        isov = tuple(iso[v] if v in iso else 0 for v in range(graph.vcount()))
        ctoh = [tuple(sorted((iso[a], iso[b]))) + (j['nonint'], ) for a, b, j in edgelist]
        edgemap = [groups[n][edgeindex[i]] for i, n in enumerate(ctoh)]
        mappings[isov] += [[edgemap[i] for i in p] for p in edgemappings]

    fe = next(e.index for e in graph.es if 0 in e.tuple)
    fl = next(e.index for e in graph.es if graph.vcount() - 1 in e.tuple)

    def dotpreserve(m, ds):
        return {
            tuple(sorted((m[d[0]] if m[d[0]] != fl else fe,
                          m[d[1]] if m[d[1]] != fl else fe))) for d in ds
        } == set(ds)

    def dotmap(m):
        res = []
        for i, d in enumerate(graph['dots']):
            newdot = tuple(sorted((m[d[0]] if m[d[0]] != fl else fe, m[d[1]] if m[d[1]] != fl else fe)))
            res.append(graph.ecount() + (i if newdot not in graph['dots'] else graph['dots'].index(newdot)))
        return res

    dotpreservemapping = {}
    # generate all combinations of dots that are nonzero and check if we find new automorphisms
    for ds in range(0, len(graph['dots']) + 1):
        for dotcombs in combinations(graph['dots'], ds):
            newmapping = {k: [i for i in m if dotpreserve(i, dotcombs)] for k, m in mappings.items()}
            newmapping = {k: v for k, v in newmapping.items() if len(v) > 0}

            if len([x for y in newmapping.values() for x in y]) > 1:
                dotpreservemapping[dotcombs] = {k: {tuple(m + dotmap(m)) for m in v} for k, v in newmapping.items()}

    # it could be that due to removing dots from the set we actually find fewer isos when
    # we set more dots to 0. Therefore, we merge the mappings that had fewer dots to 0 with the mapping
    """
    completemap = {}


    for k, v in dotpreservemapping.items():
        newv = {}
        for kk, vv in dotpreservemapping.items():
            if set(k) <= set(kk):
                for kkk, vvv in vv.items():
                    if kkk in newv:
                        newv[kkk] |= vvv
                    else:
                        newv[kkk] = vvv.copy()
        completemap[k] = newv
    """

    # filter all combinations that happen in simpler variants
    for k, v in dotpreservemapping.items():
        for kk, vv in dotpreservemapping.items():
            if k != kk and set(k) <= set(kk):
                for kkk, vvv in vv.items():
                    if kkk in v:
                        v[kkk] -= vvv

    dotpreservemapping = {k: v for k, v in dotpreservemapping.items() if
        any(len(x) > 0 for x in v.values())}

    return list(sorted(dotpreservemapping.items(), key=lambda x: len(x[0])))


def isIsomorphic(c, h, Qepsilons=True):
    """Checks if two graphs are isomorphic. It does not check the dot products.
    The number of epsilons in the Q are checked if Qepsilons is true.
    """
    vd = c.vcount() - h.vcount()
    if vd < 0:
        c.add_vertices(-vd)
    elif vd > 0:
        h.add_vertices(vd)

    if len(set(i for e in c.es for i in e.tuple)) != len(set(i for e in h.es for i in e.tuple)):
        return False, []

    if Qepsilons:
        if sum(e['nonint'] for e in c.es if isExternal(c, e)) != \
            sum(e['nonint'] for e in h.es if isExternal(h, e)):
            return False, []

    r = c.isomorphic_vf2(h, return_mapping_12=True, edge_compat_fn=areEdgesIsomorphic)
    if vd < 0:
        c.delete_vertices(range(c.vcount() + vd, c.vcount()))
    elif vd > 0:
        h.delete_vertices(range(h.vcount() - vd, h.vcount()))
    return r


def findTopologyActionAndDots(g, doDots=True, dotsug=[], order=('insertion', 'carpet', 'triangle', 'diamond')):
    """For a topology g, find an appropriate action and appropriate dots.
    An action could be solving an insertion, using the triangle rule or a manual
    solution.
    """
    g['carpet'] = []
    g['triangle'] = []
    g['diamond'] = []
    g['insertion'] = []

    for action in order:
        if action == 'diamond':
            g['diamond'] = findDiamond(g)
            if not g['diamond']:
                continue

            # now find dot products
            if doDots:
                # find a dot product system that does not depend on any of the internal lines of the diamond
                g['dots'] = findDots(g, g['diamond'][0], dotsug, True)
            return
        elif action == 'carpet':
            # find out if we have a connection between external lines
            extv1 = g.neighbors(0)[0]
            extv2 = g.neighbors(g.vcount() - 1)[0]
            if extv1 in g.neighbors(extv2) and g.ecount() > 1:
                r = next((e.index for e in g.es if e.tuple == tuple(sorted((extv1, extv2)))), None)
                if r is None:
                    continue

                # find a third line that could be involved
                # this line has momentum consisting of only Q and r
                pp, sps, dd_reps, dpind, dp, dpmap, cachekey = solvePropagatorSystem(g)
                fullrep = dd_reps.copy()
                fullrep.update({x: x for x in pp if x not in dd_reps})

                third = fullrep[pp[g.get_eid(0, extv1)]] - fullrep[pp[r]] if extv1 < extv2 \
                            else fullrep[pp[g.get_eid(0, extv1)]] + fullrep[pp[r]]

                third = next((pp.index(k) for k, v in fullrep.items() if abs(v) == abs(third)), None)

                # number of epsilons in subloop
                eps = sum(e['nonint'] for e in g.es if e.index not in {r, third, g.get_eid(0, extv1),
                    g.get_eid(extv2, g.vcount() - 1)})

                # number of loops in subloop (one minus total loops)
                lc = g.ecount() - len({v for e in g.get_edgelist() for v in e})

                g['carpet'] = (third, r, eps, lc)

                # we can view a carpet as an insertion between Q and r
                idsind = [g.get_eid(0, extv1), r]
                ignore = [r] + ([] if third is None else [third])

                dots = findDots(g, ignore, dotsug, True, idsind[0]) if doDots else g['dots']
                g['dots'] = list(sorted(dots, key=lambda x: 1 if idsind[0] in x else 0))
                g['insertion'] = (idsind, [d for d in dots if idsind[0] in d])
                return
        elif action == 'triangle':
            # find all triangles
            threeneighb = [x for x in range(g.vcount()) if len(g.neighbors(x)) == 3
                           and len({0, g.vcount() - 1} & set(g.neighbors(x))) == 0]
            edge = [(a, b) for a in threeneighb for b in threeneighb if g.are_connected(a, b) and a < b]

            triangles = []
            for a, b in edge:
                c = set(g.neighbors(a)) & set(g.neighbors(b))
                triangles += [(a, b, x) for x in c]

            fulltri = []
            for a, b, c in triangles:
                # remove the extra line at a in (a,b,c)
                d = list(set(g.neighbors(a)) - {b, c})[0]
                e = list(set(g.neighbors(b)) - {a, c})[0]

                # check if the lines are integer
                if g.es[g.get_eid(a, b)]['nonint'] != 0 or g.es[g.get_eid(a, d)]['nonint'] != 0 or \
                                g.es[g.get_eid(b, e)]['nonint'] != 0:
                    continue

                # store the edges of the triangle
                fulltri.append([g.get_eid(*x) for x in [(a, c), (b, c), (a, b), (a, d), (b, e)]])

            # find a triangle of which the sum of the internal indices are lowest
            if len(fulltri) > 0:
                g['triangle'] = min(fulltri, key=lambda x: x[0] + x[1] + x[2])

                if doDots:
                    # find a dot product system that does not depend on any of the internal lines of the triangle
                    g['dots'] = findDots(g, [g['triangle'][0], g['triangle'][1]], dotsug, True)
                return
        elif action == 'insertion' and not g.is_simple():
            el = list(sorted(g.get_edgelist()))
            edge = next(x.tuple for x in g.es if el.count(x.tuple) > 1)
            ids = g.es.select(_source=edge[0], _target=edge[1])[:2]  # pick the first pair
            idsind = [i.index for i in ids]
            ignore = [idsind[1]]

            # re-express dots so that they do not depend on the insertion line that will be removed
            if doDots:
                dots = findDots(g, ignore, dotsug, True, idsind[0]) if doDots else g['dots']
                # sort the dots so that dots unaffected by the insertion are first
                g['dots'] = list(sorted(dots, key=lambda x: 1 if idsind[0] in x else 0))
            g['insertion'] = (idsind, [d for d in g['dots'] if idsind[0] in d])
            return

    # a manual solution is required
    if doDots:
        g['dots'] = findDots(g, [], True)


def renameTopologies(metagraph):
    """Renames all topologies so that the number indicates the order of execution.
    Also labels the transitional graphs for insertions.
    """
    tl = [v for v in range(metagraph.vcount()) if len(metagraph.neighbors(v, mode=IN)) == 0]  # top level
    queue = [[] for _ in range(metagraph.vs[tl[0]]['graph'].ecount() + 1)]
    queue[-1] = tl

    processed = set()

    namemap = {}
    newcount = metagraph.vcount()  # we give the topologies new names

    while True:
        for edgeloop in reversed(range(len(queue))):
            if len(queue[edgeloop]) > 0:
                c = queue[edgeloop].pop()
                break
        else:
            break  # done

        if metagraph.vs[c]['graph']['name'][0] == 'd':
            namemap[c] = 'd' + str(newcount)
        newcount -= 1

        for e in metagraph.es.select(_source=c):
            n = e.tuple[1]
            if n not in processed:
                queue[metagraph.vs[n]['graph'].ecount()].append(n)
                processed.add(n)

    for v in metagraph.vs:
        if v.index in namemap:
            v['name'] = namemap[v.index]
        v['graph']['name'] = v['name']
        if 'transgraphs' in v['graph'].attributes():
            for t in v['graph']['transgraphs']:
                t['name'] = v['name'] + t['namesuffix']


def createTree(diagrams, destination='reductionflow', order=('insertion', 'carpet', 'triangle')):
    """Creates a tree of how the top level diagrams break down into simpler diagrams, all
    the way to a simple propagator.
    """
    # keep track of how the graph g splits up: build a tree
    metagraph = Graph(directed=True)
    metagraph['name'] = destination
    metagraph.add_vertices(g['name'] for g in diagrams)
    metagraph.vs['graph'] = diagrams

    db = list(diagrams)  # database of all visited diagrams
    while len(diagrams) > 0:
        g = diagrams.pop()

        # for visualisation
        if g['triangle']:
            metagraph.vs.find(g['name'])['style'] = ' '
            metagraph.vs.find(g['name'])['color'] = 'black'
        elif g['carpet']:
            metagraph.vs.find(g['name'])['color'] = 'cyan'
            metagraph.vs.find(g['name'])['style'] = 'filled'
        elif g['insertion']:
            metagraph.vs.find(g['name'])['color'] = 'green'
            metagraph.vs.find(g['name'])['style'] = 'filled'
        elif g['diamond']:
            metagraph.vs.find(g['name'])['color'] = 'yellow'
            metagraph.vs.find(g['name'])['style'] = 'filled'
        else:
            metagraph.vs.find(g['name'])['color'] = 'red'
            metagraph.vs.find(g['name'])['style'] = 'filled'

        children = []  # stores all child topologies

        # add insertion child
        if g['insertion']:
            children.append(removeInsertions(g))

        # remove any non-insertion line and add as a child
        ext = {0, g.vcount() - 1}
        for e in g.es:
            if e.tuple[0] in ext or e.tuple[1] in ext or e['nonint'] != 0 \
                    or g.get_edgelist().count(e.tuple) > 1:
                continue

            h = removeEdge(g, e.index)
            h['did_insertion'] = False

            # skip tadpoles
            if not findTadpole(h, e.tuple[0]):
                children.append(h)

        for h in children:
            # check for isomorphisms
            for c in db:
                if c.ecount() != h.ecount():
                    continue

                m = isIsomorphic(c, h)

                if m[0]:
                    # h is isomorphic to c, which means that it can use the same dot products. We need to find the
                    # mapping from the old dot products of g to the new dot products of h/c.
                    edgemap = mapEdge(c, h, m[1])
                    ctoh = [tuple(sorted((edgemap[d[0]], edgemap[d[1]]))) for d in c['dots']]

                    # do we have an insertion?
                    if h['did_insertion']:
                        # for insertions, we create a dotmap inside h
                        # this is less efficient than the line removal case, where we find a map from g
                        # to h
                        h['old_dots'] = h['dots']  # back up the old dots
                        b = findBestMomentumBasis(h, ctoh, h['dots'])  # FIXME: is this correct?
                        varmap, sol = solveDotProductSystem(h, ctoh, basis=b)
                        dotmap = remapdots(h['dots'], varmap, sol)
                    else:
                        # we map h back to g to find how the dots are linked by undoing the line removal
                        gtoh = mapAfterRemoval(g.ecount(), [h['removed_edge']])
                        htog = [tuple(sorted((gtoh.index(d[0]), gtoh.index(d[1])))) for d in ctoh]
                        b = findBestMomentumBasis(g, htog, g['dots'] + [(h['removed_edge'], h['removed_edge'])],
                                disregard=[h['removed_edge']])
                        varmap, sol = solveDotProductSystem(g, htog, [h['removed_edge']], b)
                        dotmap = remapdots(g['dots'] + [(h['removed_edge'], h['removed_edge'])], varmap, sol)
                        dotmap = fixDotsAfterLineRemoval(g, dotmap, h['removed_edge'])

                    h['dots'] = ctoh

                    metagraph.add_edge(g['name'], c['name'], mapping=m[1], dotmap=dotmap, transgraph=h)
                    break
            else:
                # we have a new topology
                h['name'] = 'd' + str(len(db))
                diagrams.append(h)
                db.append(h)  # TODO: add canonical form instead?

                if not h['did_insertion'] and g.ecount() > 1:
                    # map the dots from g to dots in h
                    # this is a suggestion for the dot generation
                    gtoh = mapAfterRemoval(g.ecount(), [h['removed_edge']])
                    dotsug = [(gtoh[e1], gtoh[e2]) for e1, e2 in g['dots']]

                    # this assigns a topology action and dots for h
                    findTopologyActionAndDots(h, dotsug=dotsug, order=order)

                    # find a dotmap from g to h
                    htog = [tuple(sorted((gtoh.index(d[0]), gtoh.index(d[1])))) for d in h['dots']]
                    b = findBestMomentumBasis(g, htog, g['dots'] + [(h['removed_edge'], h['removed_edge'])],
                        disregard=[h['removed_edge']])
                    varmap, sol = solveDotProductSystem(g, htog, [h['removed_edge']], b)
                    dotmap = remapdots(g['dots'] + [(h['removed_edge'], h['removed_edge'])], varmap, sol)
                    dotmap = fixDotsAfterLineRemoval(g, dotmap, h['removed_edge'])
                else:
                    h['old_dots'] = h['dots']
                    findTopologyActionAndDots(h, dotsug=h['old_dots'], order=order)  # overwrites dots
                    b = findBestMomentumBasis(h, h['dots'], h['old_dots'])
                    varmap, sol = solveDotProductSystem(h, h['dots'], basis=b)
                    dotmap = remapdots(h['old_dots'], varmap, sol)

                h['auto'] = getAutomorphisms(h)  # add all automorphisms
                metagraph.add_vertex(h['name'], graph=h)
                metagraph.add_edge(g['name'], h['name'], dotmap=dotmap, mapping=list(range(g.vcount())), transgraph=h)

    renameTopologies(metagraph)

    metagraph.vs['label'] = metagraph.vs['name']
    metagraph.write_dot(metagraph['name'] + '.dot')
    return metagraph


def findTadpoles(metagraph, c, f):
    """Identify what lines have to be missing in order to create a tadpole.
    """
    if 'tadp' in metagraph.vs[c].attributes() and metagraph.vs[c]['tadp'] is not None:
        return metagraph.vs[c]['tadp']

    curgraph = metagraph.vs[c]['graph']
    cgi = mapAfterRemoval(curgraph.ecount(), [e.index for e in curgraph.es
                                          if 0 in e.tuple or curgraph.vcount() - 1 in e.tuple])
    curt = set()
    for e in metagraph.es.select(_source=c):
        n = e.tuple[1]
        transgraph = e['transgraph']  # transitional graph
        findTadpoles(metagraph, n, f)

        if transgraph['did_insertion']:
            continue

        gtoh = mapAfterRemoval(curgraph.ecount(), [transgraph['removed_edge']])
        ctot = mapEdge(metagraph.vs[n]['graph'], transgraph, e['mapping'])
        # translate the tadpoles back to current graph
        curt |= {tuple(sorted([gtoh.index(ctot[t]) for t in tl] + [transgraph['removed_edge']])) for tl in metagraph.vs[n]['tadp']}

    # find tadpoles in current topology
    rem = {e['transgraph']['removed_edge'] for e in metagraph.es.select(_source=c)
        if 'removed_edge' in e['transgraph'].attributes() and not e['transgraph']['did_insertion']}
    curt |= {(e.index,) for e in curgraph.es if e.index not in rem and e['nonint'] == 0
                and not isExternal(curgraph, e)}

    # filter subsets from other rules:
    curt = list(sorted(curt, key=lambda x: len(x)))
    for i, x in enumerate(curt):
        for j, y in reversed(list(enumerate(curt))):
            if j == i:
                break
            if set(x) <= set(y):
                del curt[j]

    def printtadpoles(curt):
        if len(curt) == 1:
            if len(curt[0]) == 1:
                f.write("if (count(`i'{},1) >= 0) Discard;\n".format(curt[0][0]))
            else:
                f.write('if (' + '&&'.join("(count(`i'{},1) >= 0)".format(t) for t in curt[0]) + ') Discard;\n')
        else:
            f.write("if (" +
                '||'.join('(' + '&&'.join("(count(`i'{},1) >= 0)".format(t) for t in tl) + ')' for tl in curt) + ') Discard;\n')

    if len(curt) > 0:
        f.write('#case M{}\n'.format(metagraph.vs[c]['name']))
        for e in metagraph.es.select(_target=c):
            if e['transgraph']['did_insertion']:
                f.write('#case M{}d\n'.format(metagraph.vs[e.source]['graph']['name']))
            else:
                f.write('#case M{}M{}\n'.format(metagraph.vs[e.source]['graph']['name'], cgi[e['transgraph']['removed_edge']] + 1))

        printtadpoles(curt)

        # new code: too slow when minterms is large
        """
        # TODO: make use of ignore?
        minterms = [[[1] if cgi.index(x) in tl else [0,1] for x in range(curgraph.ecount() - 2)] for tl in curt]
        minterms = [t for tl in minterms for t in list(product(*tl))]

        # TODO: change the order in and and or?
        terms = SOPform(['(count(i{},1) >=0)'.format(i + 1) for i in range(curgraph.ecount() - 2)], minterms)
        print(terms)
        if terms.func == Or:
            a = '||'.join('(' + (str(term) if term.func != And else
                    '&&'.join('({})'.format(t) for t in term.args)) + ')' for term in terms.args)
        else:
            if terms.func == And:
                a = '&&'.join('({})'.format(t) for t in terms.args)
            else:
                a = terms


        # better, but even slower...
        #terms = POSform(['(count(i{},1) >=0)'.format(i + 1) for i in range(curgraph.ecount() - 2)], minterms)
        #b = '||'.join('(' + '&&'.join('({})'.format(t) for t in term.args) + ')' for term in terms.args)

        f.write('if ({}) Discard;\n'.format(a))
        print(curgraph['name'])
        """
        f.write('#break\n')

    # for insertions, we map the dots back to the isograph
    for e in metagraph.es.select(_source=c):
        if e['transgraph']['did_insertion'] and len(metagraph.vs[e.target]['tadp']) > 0:
            f.write('#case M{}c\n'.format(curgraph['name']))
            # the target should be computed at this point
            ctot = mapEdge(metagraph.vs[e.target]['graph'], e['transgraph'], e['mapping'])
            curtf = [[ctot[t] for t in tl] for tl in metagraph.vs[e.target]['tadp']]
            printtadpoles(curtf)
            f.write('#break\n')

    metagraph.vs[c]['tadp'] = curt


def createTadpoleFilter(metagraph, filename='tadpole.prc'):
    with open(filename, 'w') as f:
        f.write("#procedure tadpole(TOPO,i)\n")
        f.write("if (count(`TOPO',1));\n")
        f.write("#switch `TOPO'\n")

        tl = [v for v in range(metagraph.vcount()) if len(metagraph.neighbors(v, mode=IN)) == 0]

        for t in tl:
            findTadpoles(metagraph, t, f)

        f.write('#endswitch\n')
        f.write('endif;\n')
        f.write("#endprocedure\n")


def generateMomentumConversion(graphs, filename='momsubs.prc'):
    """Generate a conversion from the physics results expressed in momenta p1.p1,p1.p2,... to forcer notation.
    """
    print('Generating momentum conversion')

    with open(filename, 'w') as f:
        # gs = []  # never used
        f.write("#procedure momsubs(TOPO,MOMS,DOTS)\n")
        f.write("#ifndef `SUBTOPO'\n\t#define SUBTOPO \"\"\n#endif\n")
        f.write("#switch `TOPO'\n")

        for gc, g in enumerate(graphs):
            if gc % 10 == 0:
                print('Status: {}/{}'.format(gc, len(graphs)))
            if g.ecount() == 1:
                continue

            f.write("*--#[ {0} :\n#case M{0}\n".format(g['name']))

            ext = [e.index for e in g.es if 0 in e.tuple or g.vcount() - 1 in e.tuple]
            cgi = mapAfterRemoval(g.ecount(), ext)

            # do a shared case
            f.write("#procedure SimpleDot{}()\n".format(g['name']))
            f.write("if (count(M{},1));\n".format(g['name']))
            f.write("\tid Q.Q^n? = `Q2'^n;\n")
            for i in range(g.ecount() - 2):
                f.write("\tid p{0}.p{0}^n? = i{0}^n;\n".format(i + 1))
            for i, (e1, e2) in enumerate(g['dots']):
                f.write("\tid {0}.{1}^n? = i{2}^n;\n".format(
                    'Q' if cgi[e1] == len(cgi) else 'p' + str(cgi[e1] + 1),
                    'Q' if cgi[e2] == len(cgi) else 'p' + str(cgi[e2] + 1),
                    g.ecount() - 1 + i))
            f.write("endif;\n")
            f.write("#endprocedure\n")

            """
            First, we find out all P flows through the diagram
            """
            # find shortest distance from all points to all points
            pathmap = [[0 for x in range(g.ecount())] for x in range(g.ecount())]

            for x in range(g.ecount()):
                queue = [x]
                while len(queue) > 0:
                    c = queue.pop(0)
                    for a in [e.index for e in g.es if
                            len(set(g.es[c].tuple) & set(e.tuple)) == 1
                            and pathmap[x][e.index] == 0 and e.index != x]:
                        pathmap[x][a] = pathmap[x][c] + 1
                        queue.append(a)

            # generate shortest paths
            paths = []
            for source in range(g.ecount()):
                for dest in range(source, g.ecount()):
                    p = [source]
                    while p[-1] != dest:
                        # convdia should always pick the first option
                        minv = min(pathmap[x][dest] for x in range(g.ecount()) if pathmap[p[-1]][x] == 1)
                        opt = [x for x in range(g.ecount()) if pathmap[p[-1]][x] == 1 and pathmap[x][dest] == minv]
                        p.append(opt[0])
                    paths.append(tuple(p))

            last = next(e.index for e in g.es if g.vcount() - 1 in e.tuple)

            # now filter double paths due to symmetry
            emaps = set(em for d, ms in g['auto'] for vm, ems in ms.items() for em in ems) | \
                {tuple(range(g.ecount()))}
            newpaths = set()
            for p in paths:
                isos = [tuple(iso[a] if iso[a] != last else 0 for a in p) for iso in emaps]
                isos += [tuple(reversed(iso)) for iso in isos]
                isos = list(sorted(tuple(iso) for iso in isos))

                if len(isos) > 0 and isos[0] not in newpaths and isos[0].count(0) <= 1 and isos[0] != (0,):
                    newpaths.add(isos[0])

            """
            Now find a suitable momentum basis, aligned with the P flow and with the dot products
            """
            basisdict = collections.defaultdict(list)
            for path in list(sorted(newpaths)) + [[]]:  # empty path is used for default
                basis = findBestMomentumBasis(g, g['dots'], None, disregard=[], include=set(path), exclude=False)
                if basis is not None:
                    basisdict[basis].append(tuple(path))

            f.write("#switch `SUBTOPO'\n")
            for basis, paths in sorted(basisdict.items(), key=lambda x: 100 if () in x[1] else 0):
                for path in sorted(paths):  # empty path is used for default
                    if path != ():
                        # find maximum distance to find start and end of path
                        f.write("#case {}{}\n".format('Q' if path[0] == 0 else path[0], path[-1]))
                        f.write("* path: {}\n".format(path))
                    else:
                        # TODO: use dot momenta as include for default case?
                        f.write('#default\n')
                        #f.write("#message Using default for subtopo `SUBTOPO' in {}\n".format(g['name']))
                f.write("#if ( `MOMS' > 0 )\n")

                pp, _, dd_reps, _, _, _, _ = solvePropagatorSystem(g, basis=basis)
                independent = [pp[i] for i in range(len(pp)) if pp[i] not in dd_reps.keys()]
                depp = [pp.index(x) for x in dd_reps.keys()]
                dependent = [pp[i] for i in range(len(pp)) if pp[i] in dd_reps.keys() and
                            i != next(e.index for e in g.es if g.vcount() - 1 in e.tuple)]

                # generate momentum conservation
                verts = [v for v in sorted(range(g.vcount()), key=lambda x: len(g.neighbors(x)))
                            if len(g.neighbors(v)) > 2]
                eqv = []
                for i in verts:
                    eqs = []
                    for j in set(g.neighbors(i)):
                        sign = 1 if j >= i else -1
                        for e in g.es:
                            if e.tuple == tuple(sorted((i, j))):
                                eqs.append((sign, pp[ext[0]] if e.index == ext[1] else pp[e.index]))
                    eqv.append(eqs)

                # find all systems that can solve all dependent variables
                queue = [([], [])]
                schemes = []
                while len(queue) > 0:
                    curi, curp = queue.pop()

                    if len(curp) == len(dependent):
                        schemes.append((curi, curp))
                        continue

                    upperindex = 0 if curi == [] else curi[-1] + 1
                    for c, p in eqv[upperindex]:
                        if p not in independent + curp:
                            queue.append((curi + [upperindex], curp + [p]))

                bestscheme = min(schemes, key=lambda x: sum(len(eqv[a])**2 for a in x[0]))

                # reorganize the vertex rules and see if you can make replacements for the
                # dependent variables straight into the independent ones, if this does not generate more terms
                depchain = []
                for i, p in zip(*bestscheme):
                    dep = next(aa for aa, a in enumerate(eqv[i]) if a[1] == p)
                    eqv[i][0], eqv[i][dep] = eqv[i][dep], eqv[i][0]

                    if len(dd_reps[p].args) < len(eqv[i]):  # eqv[i] also has lhs
                        depchain.append([(-1, p)] + [(1, t) if t in pp else (-1, t.args[1]) for t in dd_reps[p].args])
                    else:
                        depchain.append(eqv[i])

                # calculate maximum distance from each rewrite to the original independent basis
                # then print in order of decreasing distance
                dist = {}
                queue = [(0, x) for x in independent]

                while len(queue) > 0:
                    d, cur = queue.pop(0)

                    if cur not in dist or d > dist[cur]:
                        dist[cur] = d

                    # find neighbours
                    for x in depchain:
                        if any(cur == a[1] for a in x[1:]):
                            queue.append((d + 1, x[0][1]))

                # make sorted dependency chain
                # TODO: split off in disjoint groups
                dcc = {a[0][1]: len(a) - 1 for a in depchain}
                di = list(sorted(dist.items(), key=lambda x: str(x[0])))
                depchainsorted = list(sorted(di, key=lambda x: x[1] * 20 -
                        (0 if x[1] == 0 else dcc[x[0]]), reverse=True))

                momrewlist = [next(a for a in depchain if a[0][1] == x[0]) for x in
                                depchainsorted if x[1] > 0]
                for i, x in enumerate(momrewlist):
                    # do not rewrite the independent dot products
                    dp = ','.join(sorted(set(('Q' if d[0] == ext[0] else str(pp[d[0]])) if pp.index(x[0][1]) == d[1] else
                        ('Q' if d[1] == ext[1] else str(pp[d[1]])) for d in g['dots'] if pp.index(x[0][1]) in d)))
                    f.write('\tif (count(M{3},1)) id {0}.p?!{{{0},[P+{0}],[P-{0}]{2}}} = {1};\n'.format(
                        '{}'.format('Q' if x[0][1] == pp[ext[0]] else x[0][1]), ' '.join(
                        '- {}.p'.format('Q' if a == pp[ext[0]] else a) if s * -x[0][0] < 0 else
                        '+ {}.p'.format('Q' if a == pp[ext[0]] else a) for s, a in x[1:]),
                        ',' * (len(dp) > 0) + dp, g['name']))
                    f.write('\t.sort:momrewrite-M{}-{};\n'.format(g['name'], len(momrewlist) - i))

                f.write("#endif\n#if ( `DOTS' > 0 )\n")
                f.write("\t#call SimpleDot{}()\n".format(g['name']))

                # rewrite all momentum combinations to the basis
                varmap, sol = solveDotProductSystem(g, g['dots'], basis=basis)
                dotmap = parseMap([(a, b) for a in range(g.ecount()) for b in range(a, g.ecount())], varmap, sol)

                symbs = symbols(['i' + str(i + 1) for i in range(g.ecount() - 2 + len(g['dots']))] + ["`Q2'"])
                # strtot = []  # never used
                sortcount = 0
                # sorttable = {1: 100, 2: 3, 3: 2, 4: 1}  # never used
                for d, dm in sorted(dotmap.items(), key=lambda x: len(x[1]) * 1000 + x[0][1] * 100 + x[0][0]):
                    if ext[1] in d or len(dm) == 1 or (len(dm) > 1 and (d[0] in depp or d[1] in depp)):
                        continue

                    exp = 0
                    for coeffd, dot in dm:
                        if dot[0] == dot[1]:  # propagator increase
                            if 0 in g.es[dot[0]].tuple or g.vcount() - 1 in g.es[dot[0]].tuple:
                                exp += coeffd * symbs[-1]
                            else:
                                exp += coeffd * symbs[cgi[dot[0]]]
                        else:
                            exp += coeffd * symbs[g.ecount() - 2 + g['dots'].index(dot)]

                    if len(dm) > 2:
                        prop = ','.join(list(sorted(str(s) for s in exp.atoms(Symbol) if str(s) != '`Q2\'')))
                        f.write('\t#call replace({0}.{1},{2},i,{5},M{3},dotrewrite-{3}-{4})\n'.format(
                            'Q' if d[0] in ext else 'p' + str(cgi[d[0]] + 1),
                            'Q' if d[1] in ext else 'p' + str(cgi[d[1]] + 1),
                            str(exp).replace(' ', ''),
                            g['name'],
                            sortcount,
                            prop))
                    else:
                        f.write('\tif (count(M{},1)) id {}.{}{} = {}{};\n'.format(
                            g['name'],
                            'Q' if d[0] in ext else 'p' + str(cgi[d[0]] + 1),
                            'Q' if d[1] in ext else 'p' + str(cgi[d[1]] + 1),
                            '^n?' * (len(dm) == 1),
                            exp,
                            '^n' * (len(dm) == 1)))

                    sortcount += 1
                f.write('#endif\n#break\n')
            f.write('#endswitch\n#break\n*--#] {} : \n'.format(g['name']))
        f.write('#endswitch\n')
        f.write('#endprocedure\n')


def convertTreeToFORM(metagraph, destination='pathforcer.h'):
    """Converts a reduction graph to a FORM application.
    """
    print('Creating FORM output')

    f = open(destination, 'w')

    # print header
    f.write(textwrap.dedent("""\
        #ifndef `FORCERMAINHFILE'
        #define FORCERMAINHFILE
        *--#[ setup :

        AutoDeclare CFunction Z,Y;
        AutoDeclare Symbol M,i,j;
        Symbol Q2;
        #ifndef `SETQONE'
            #define SETQONE "1"
        #endif
        #if `SETQONE'
            #define Q2 "1"
        #else
            #define Q2 "Q2"
        #endif

        #ifndef `SORTLIMIT'
            #define SORTLIMIT "0"
        #endif

        #include- forcer-aux.h
        #include- forcer-rules.h
        #include- turbo.h

        *--#] setup : 
        *--#[ DoForcer :

        #procedure DoForcer()
    """))

    # do an initial count of all topologies
    f.write('*--#[ topocount :\n')
    special_names = ','.join([g['graph']['name'] for g in metagraph.vs
                             if g['graph']['name'][0] != 'd'])
    special_names = wrap_text(special_names, 80 - 7, 80, ',', '\\', '')

    f.write(textwrap.dedent("""\
        #do i={{{},\\
        d1,...,d{}}}
        \t#define HASTOPO`i' "0"
        \tif (count(M`i',1)) redefine HASTOPO`i' "1";
        #enddo
    """).format(special_names, metagraph.vcount()))

    f.write('.sort:topocount;\n')
    f.write('*--#] topocount : \n')

    # we output the FORM file sorted by the edge count, so that a .sort can be placed in between
    tl = [v for v in range(metagraph.vcount()) if len(metagraph.neighbors(v, mode=IN)) == 0]  # top level
    queue = [[] for _ in range(metagraph.vs[tl[0]]['graph'].ecount() + 1)]
    queue[-1] = tl

    processed = set()

    current_edge_loop = None
    lastnames = []
    while True:
        for edgeloop in reversed(range(len(queue))):
            if len(queue[edgeloop]) > 0:
                c = queue[edgeloop].pop()
                lastnames.append(metagraph.vs[c]['graph']['name'])
                break
        else:
            break  # done

        # for every edge going from c to n
        for e in metagraph.es.select(_source=c):
            n = e.tuple[1]

            if n not in processed:
                queue[metagraph.vs[n]['graph'].ecount()].append(n)
                processed.add(n)

        if current_edge_loop is None:
            current_edge_loop = max(edgeloop - 2, 0)
            f.write('*--#[ {0} edges :\n'.format(current_edge_loop))

        curgraph = metagraph.vs[c]['graph']

        # map edge id to index by filtering external lines
        cgi = mapAfterRemoval(curgraph.ecount(), [e.index for e in curgraph.es
                                                  if 0 in e.tuple or curgraph.vcount() - 1 in e.tuple])

        f.write("*--#[ reduce-{} :\n".format(curgraph['name']))
        # first, we create a function for a possible insertion
        if curgraph['insertion']:
            generateInsertion(metagraph, c, curgraph, f)

        f.write('#procedure Prepare{}()\n'.format(curgraph['name']))
        f.write('* Topology: ' + ' '.join('{}-{}{}'.format(e.source, e.target, '*' * e['nonint']) for e in curgraph.es) + '\n')
        f.write('* Dots: ' + ','.join('p{}.p{}'.format(*e) for e in curgraph['dots']) + '\n')

        if curgraph.ecount() > 1:
            generateAutomorphism(curgraph, f)
            f.write('\n')
            generateRelabeling(metagraph, c, curgraph, cgi, f)
        else:
            f.write("\tid M{} = 1;\n".format(curgraph['name']))
            f.write('.sort:reduce {};\n'.format(curgraph['name']))
        f.write("#endprocedure\n\n")

        # topohash
        """
        topohash = ['{}-{}{}'.format(e.source, e.target, '*' * e['nonint'])
                    for e in curgraph.es]
        topohash += ['p{}.p{}'.format(*e) for e in curgraph['dots']]
        if 'convert' in curgraph.attributes():
            topohash += ['convert']
        if 'mindot' in curgraph.attributes():
            topohash += ['mindot: ' + str(curgraph['mindot'])]
        topohash = ' '.join(topohash).encode(encoding='utf-8')
        topohash = hashlib.sha1(topohash).hexdigest()
        """

        # do the topology action here
        f.write("#procedure Reduce{}()\n".format(curgraph['name']))
#       f.write('\t#define topohash "{}"\n'.format(topohash))
        f.write("\t#call Prepare{}()\n".format(curgraph['name']))
        if curgraph.ecount() > 1:
            f.write("\t#if `$DOACTION'\n")
            generateAction(curgraph, cgi, f)
            f.write("\t#endif\n")
        f.write("#endprocedure\n")
        f.write("*--#] reduce-{} : \n".format(curgraph['name']))

        if len(queue[edgeloop]) == 0:
            f.write(textwrap.dedent("""\
                #call ReduceEdge({0},\\
                {1}\\
                )
                *--#] {0} edges : 
            """).format(
                current_edge_loop,
                wrap_text(','.join(n for n in lastnames),
                          80, 80, ',', '\\', '')
            ))
            current_edge_loop = None
            lastnames = []

    f.write("#call EnablePolyratfun()\n")
    f.write(".sort:DoForcer end;\n")
    f.write("#endprocedure\n")
    f.write("\n")
    f.write("*--#] DoForcer : \n")
    f.write("#endif\n")
    f.write("* vim: ft=form et ts=8 sts=2 sw=2 fdm=marker fmr=#[,#] fdc=4\n")
    f.close()


def generateAutomorphism(curgraph, f):
    """Write FORM code for applying automorphism of the graph."""
    assert curgraph.ecount() > 1

#   f.write("#$MAXNUMSUM = 0;\n")
    f.write('* Apply automorphisms\n')
    fullcgi = mapAfterRemoval(curgraph.ecount() + len(curgraph['dots']), [e.index for e in curgraph.es
                                              if 0 in e.tuple or curgraph.vcount() - 1 in e.tuple])

    f.write('if (match(M{0}*<1/i1^n1?$n1>*...*<1/i{1}^n{1}?$n{1}>));\n'.format(
        curgraph['name'], curgraph.ecount() - 2 + len(curgraph['dots'])))

    # calculate the sum of numerators
#   if len(curgraph['dots']) == 0:
#       f.write("\t$NUMSUM = 0;\n")
#   else:
#       f.write("\t$NUMSUM = {};\n".format('+'.join('abs_($n{})'.format(i) for i in range(curgraph.ecount() - 1,
#           curgraph.ecount() - 1 + len(curgraph['dots'])))))
#   f.write("\tif ($NUMSUM > $MAXNUMSUM) $MAXNUMSUM = $NUMSUM;\n")

    # count number of terms and intialize $iso
    # find the case where all dots are present: this is the baseline that can be
    # removed from all others
    # FIXME: ugly
    ignoredots = set()
    autocountmap = {}
    counter = 0
    for autoc, (dotnotzero, automap) in enumerate(sorted(curgraph['auto'])):
        for vm, ems in sorted(automap.items()):
            for em in sorted(ems):
                tr = tuple(fullcgi[em[fullcgi.index(i)]] + 1 for i in range(len(fullcgi) - 2))

                if len(dotnotzero) == len(curgraph['dots']):
                    ignoredots.add(tr)
                if tr in autocountmap:
                    continue
                autocountmap[tr] = counter
                counter += 1

    # if there is only the identity, skip
    if counter > 1:
        f.write('#do i=0,{}\n'.format(counter - 1))
        f.write("\t$iso`i'=`HASHBASE'^{};\n".format(len(fullcgi)))
        f.write('#enddo\n')

        printelse = False
        for autoc, (dotnotzero, automap) in enumerate(sorted(curgraph['auto'])):
            dotzero = set(range(curgraph.ecount() - 1, curgraph.ecount() - 1 + len(curgraph['dots']))) - \
                set(curgraph.ecount() - 1 + curgraph['dots'].index(d) for d in dotnotzero)
            dotzerocheck = '&&'.join(['($n{}==0)'.format(x) for x in dotzero])

            # handle the general case where all dots are present
            if printelse and len(dotzerocheck) == 0:
                #print('endif;', file=f)
                printelse = False

            if len(dotzerocheck) > 0:
                f.write('{}if ({});\n'.format('else' if printelse else '', dotzerocheck))
                #printelse = True

            autoset = set()  # set(ignoredots) if len(dotzerocheck) > 0 else set()
            for vm, ems in sorted(automap.items()):
                for em in sorted(ems):
                    tr = tuple(fullcgi[em[fullcgi.index(i)]] + 1 for i in range(len(fullcgi) - 2))
                    if tr in autoset:
                        continue
                    autoset.add(tr)
                    # for the triangle, the main sorting criterium is the sum of the lines that will be decreased
                    # so we add it as a first sorting index
                    firstterm = ''
                    if curgraph['triangle']:
                        #firstterm = 'abs_({})*10^15+'.format('+'.join('$n' + str(fullcgi[em[i]] + 1) for i in curgraph['triangle'][2:]))
                        firstterm = '+'.join('$n' + str(fullcgi[em[i]] + 1) for i in curgraph['triangle'][2:])

                    #print("\t$iso{}={}{};".format(autocountmap[tr],firstterm,
                    #    '+'.join('abs_($n{})*10^{}'.format(fullcgi[em[fullcgi.index(i)]] + 1, len(fullcgi) - 2 - i) for i in range(len(fullcgi) - 2))), file=f)

                    #f.write("\t$iso{}={}{}{};\n".format(autocountmap[tr],firstterm,
                    #    '+'.join('abs_($n{})*10^{}'.format(tr[i], len(fullcgi) - 2 - i) for i in range(curgraph.ecount() - 2)),
                    #    '-'*(len(curgraph['dots']) > 0) + '-'.join('$n{}*10^{}'.format(tr[i], len(fullcgi) - 2 - i) for i in range(curgraph.ecount() - 2, len(fullcgi) - 2))))

                    liststr = [firstterm] if len(firstterm) > 0 else []
                    liststr += ['$n{}'.format(tr[i]) for i in range(len(fullcgi) - 2)]

                    # collect all terms in reverse order
                    f.write("\t#call hash({},{})\n".format(autocountmap[tr], ','.join(reversed(liststr))))

            if len(dotzerocheck) > 0:
                f.write('endif;\n')
        if printelse:
            f.write('endif;\n')

        f.write('$bestiso=-1;\n')
        f.write("$bestisov=`HASHBASE'^{};\n".format(len(fullcgi)))
        f.write("#do i=0,{0}\nif ($iso`i' < $bestisov);\n\t$bestiso=`i';\n\t$bestisov=$iso`i';\nendif;\n#enddo\n".format(counter - 1))

        printelse = False
        autoset = set()
        for autoc, (dotnotzero, automap) in enumerate(sorted(curgraph['auto'])):
            for vm, ems in sorted(automap.items()):
                signchange = []
                for i, di in enumerate(curgraph['dots']):
                    d = (curgraph.es[di[0]].tuple, curgraph.es[di[1]].tuple)

                    if (vm[d[0][0]] > vm[d[0][1]]) != (vm[d[1][0]] > vm[d[1][1]]):  # sign change
                        signchange.append('$n' + str(curgraph.ecount() - 1 + i))
                for em in sorted(ems):
                    tr = tuple(fullcgi[em[fullcgi.index(i)]] + 1 for i in range(len(fullcgi) - 2))
                    if tr in autoset:
                        continue
                    autoset.add(tr)
                    coeff = '' if not signchange else 'sign_(' + '+'.join(signchange) + ')*'
                    r = ','.join('i{},i{}'.format(i + 1, fullcgi[em[fullcgi.index(i)]] + 1) for i in range(len(fullcgi) - 2)
                            if i != fullcgi[em[fullcgi.index(i)]])

                    if len(r) > 0:
                        f.write('{}if ($bestiso == {}); Multiply {}replace_({});\n'.format('else' * printelse, autocountmap[tr], coeff, r))
                        printelse = True
        f.write('endif;\n')

    f.write('endif;\n')
#   f.write('ModuleOption maximum,$MAXNUMSUM;\n')
#   f.write('ModuleOption local,<$n1>,...,<$n{}>{},$NUMSUM{};\n'.format(
    f.write('ModuleOption local,<$n1>,...,<$n{}>{}{};\n'.format(
        curgraph.ecount() - 2 + len(curgraph['dots']),
        ',<$iso0>,...,<$iso{}>'.format(counter - 1) if counter > 0 else '',
        ',$bestiso,$bestisov' if counter > 0 else ''))
    f.write('.sort:Z{} automorphism;\n'.format(curgraph['name']))


def generateRelabeling(metagraph, c, curgraph, cgi, f):
    """Write FORM code for handling missing lines."""
    assert curgraph.ecount() > 1

    #print("#do i=1,{}\n\t#$M`i'=0;\n#enddo\n\nif (count(M{},1));".format(curgraph.ecount() - 2, curgraph['name']), file=f)

    # filter tadpoles. all topologies that reach this point and still have a propagator
    # with a negative power are tadpoles
    # get all edges that have transitions if removed:
    rem = {e['transgraph']['removed_edge'] for e in metagraph.es.select(_source=c)
        if 'removed_edge' in e['transgraph'].attributes() and not e['transgraph']['did_insertion']}

    # make an exception for propagators with an epsilon in the power and edges that are
    # already checked. If the topology action is an insertion, ignore ALL edges.
    tadp = [str(cgi[e.index] + 1) for e in curgraph.es if e.index not in rem and e['nonint'] == 0
        and 0 not in e.tuple and curgraph.vcount() - 1 not in e.tuple]
    skip = [str(cgi[e.index] + 1) for e in curgraph.es if ((e.index not in rem and e['nonint'] > 0) or
            (curgraph['insertion'] and (e.index in rem or e['nonint'] > 0))) and 0 not in e.tuple and curgraph.vcount() - 1 not in e.tuple]
    specials = ','.join(skip + tadp)
    f.write('#call CheckRemovedPropagators({},{},{}{})\n'.format(curgraph['name'], curgraph.ecount() - 2, len(skip),
        ',' * (len(specials) > 0) + specials))

    """
    printelse = ''
    if len(checkind) > 0:
        print("if ({});\n\tDiscard;".format(" || ".join('(count(i{},1) >= 0)'.format(x) for x in checkind)), file=f)
        printelse = 'else'

    for e in metagraph.es.select(_source=c):
        if not e['transgraph']['did_insertion']:
            rid = cgi[e['transgraph']['removed_edge']] + 1
            f.write("{2}if (count(i{0},1) >= 0);\n\tMultiply replace_(M{1},M{1}M{0});\n\t$M{0}=1;\n".format(rid, curgraph['name'], printelse))
            printelse = 'else'

    print('{}$DOACTION = 1;{}'.format(printelse + ";\n\t" if printelse != '' else '', "\nendif;" if printelse != '' else ''), file=f)
    print("endif;\nModuleOption maximum,<$M1>,...,<$M{}>;\n.sort:relabel;".format(curgraph.ecount() - 2), file=f)
    """

    f.write("#if `$DOACTION' == 0\n")

    # for every edge going from c to n, write the topology change
    for e in metagraph.es.select(_source=c):
        oldname = curgraph['name']  # name could be changed

        n = e.tuple[1]

        transgraph = e['transgraph']  # transitional graph

        # skip all insertions, we always do them
        if curgraph['insertion'] or transgraph['did_insertion']:
            continue

        f.write('\n')

        rid = cgi[transgraph['removed_edge']] + 1
        curgraph['name'] = curgraph['name'] + 'M' + str(rid)  # change the name
        gtoh = mapAfterRemoval(curgraph.ecount(), [transgraph['removed_edge']])

        f.write("#if `$M{}'\n".format(rid))

        # re-express the dot products and edges from curgraph to transgraph to newgraph
        if len(e['dotmap']) > 0:
            f.write(generateDotmap(e['dotmap'], curgraph, transgraph,
                    gtoh, metagraph.vs[n]['graph'], e['mapping']))

        f.write("#endif\n")

        curgraph['name'] = oldname  # restore name

    f.write('.sort:reduce {};\n'.format(curgraph['name']))
    f.write("#endif\n")


def generateAction(curgraph, cgi, f):
    """Write FORM code for the topology action."""
    if curgraph['triangle'] or curgraph['diamond'] or curgraph['insertion']:
        f.write("#ifdef `HASSCHEMEZ{}'\n".format(curgraph['name']))

    if curgraph.ecount() > 1:
        if 'mindot' in curgraph.attributes():
                print('\t\tif (count(M{0},1)) Multiply replace_(i{1},-i{1});'.format(curgraph['name'], curgraph['mindot']), file=f)

        print("\t\t#call MtoZ(M{0},Z{0},i,{1})".format(curgraph['name'], curgraph.ecount() - 2 + len(curgraph['dots'])), file=f)

        if 'convert' in curgraph.attributes():
            print('\t\tid Z{0}(n1?,...,n14?) = Z{0}(n1,...,n14) * {1}2^({2});\n\t\t.sort:convert {0};'.format(curgraph['name'],
                '', '+'.join(['n{}'.format(i) for i in range(curgraph.ecount() - 1, curgraph.ecount() - 1 + len(curgraph['dots']))])), file=f)

        print("\t\t#call Reduce({})".format(curgraph['name']), file=f)

        if 'convert' in curgraph.attributes():
            print('\t\tid Z{0}(n1?,...,n14?) = Z{0}(n1,...,n14) * {1}2^({2});'.format(curgraph['name'],
                '1/', '+'.join(['n{}'.format(i) for i in range(curgraph.ecount() - 1, curgraph.ecount() - 1 + len(curgraph['dots']))])), file=f)

        print("\t\t#call ZtoM(Z{0},M{0},i,{1})".format(curgraph['name'], curgraph.ecount() - 2 + len(curgraph['dots'])), file=f)

        if 'mindot' in curgraph.attributes():
                print('\t\tif (count(M{0},1)) Multiply replace_(i{1},-i{1});'.format(curgraph['name'], curgraph['mindot']), file=f)

    if curgraph['insertion']:
        f.write("\t\t#call Prepare{}()\n".format(curgraph['name']))

    if curgraph['triangle'] or curgraph['diamond'] or curgraph['insertion']:
        f.write("#else\n")
        f.write("\t\t#if `USETURBO'\n")
        f.write("\t\t\t#call TurboM{}()\n".format(curgraph['name']))
        f.write("\t\t#endif\n")

    if curgraph['triangle']:
        trinames = ['i' + str(cgi[i] + 1) for i in curgraph['triangle']]
        dot = ['i' + str(curgraph.ecount() - 1 + i) for i, d in enumerate(curgraph['dots']) if curgraph['triangle'][2] in d]
        f.write("\t\t#call triangleS(M{},{}{},{},{})\n".format(
            curgraph['name'],
            ','.join(trinames),
            ',' * (len(dot) > 0) + ','.join(dot),
            curgraph.es[curgraph['triangle'][0]]['nonint'],
            curgraph.es[curgraph['triangle'][1]]['nonint']))

    if curgraph['diamond']:
        dianames = [['i' + str(cgi[i] + 1) for i in k] for k in curgraph['diamond']]
        dotopt = curgraph['diamond'][1] + curgraph['diamond'][3]  # dots can be from B or S
        dot = ['i' + str(curgraph.ecount() - 1 + i) for i, d in enumerate(curgraph['dots']) if d[0] in dotopt or
            d[1] in dotopt]
        f.write("\t\t#call diamondS(M{},{},{},{}{},{})\n".format(
            curgraph['name'],
            len(curgraph['diamond'][1]),
            len(curgraph['diamond'][3]),
            ','.join(','.join(x) for x in dianames if x != []),
            ',' * (len(dot) > 0) + ','.join(dot),
            ','.join(str(curgraph.es[e]['nonint']) for e in curgraph['diamond'][1])))

    # handle insertion and carpet
    if curgraph['insertion']:
        # propagator location in the old system
        prop = [d[0] if d[0] != curgraph['insertion'][0][0] else d[1] for d in curgraph['insertion'][1]]
        prop = ["`Q2'" if cgi[m] == len(cgi) else 'i' + str(cgi[m] + 1) for m in prop]

        # extra dots can be inferred since the number of indices is known
        extradotstart = curgraph.ecount() - 1 + len(curgraph['dots']) - len(curgraph['insertion'][1])
        numnewdots = len(curgraph['insertion'][1]) * (len(curgraph['insertion'][1]) - 1) // 2
        olddots = ','.join(['i{}'.format(x) for x in range(extradotstart, curgraph.ecount() - 1 + len(curgraph['dots']))])

        if curgraph['carpet']:
            prop = [d[0] if d[0] != curgraph['insertion'][0][0] else d[1] for d in curgraph['insertion'][1]]
            prop = [0 if cgi[m] == len(cgi) else str(cgi[m] + 1) for m in prop]

            f.write("\t\t#call MtoZ(M{0},Z{0},i,{1})\n".format(curgraph['name'], curgraph.ecount() - 2 + len(curgraph['dots'])))
            f.write("\t\t#call carpet(Z{},{},{},{},{},{},{},{}{})\n".format(curgraph['name'],
                            0 if not curgraph['carpet'][0] else cgi[curgraph['carpet'][0]] + 1,
                            cgi[curgraph['carpet'][1]] + 1,
                            0 if not curgraph['carpet'][0] else curgraph.es[curgraph['carpet'][0]]['nonint'],
                            curgraph.es[curgraph['carpet'][1]]['nonint'],
                            curgraph['carpet'][2],  # epsilons
                            curgraph['carpet'][3],  # loop count
                            extradotstart,
                            ',' * (len(prop) > 0) + ','.join(prop)))
            f.write("\t\t#call ZtoM(Z{0}b,M{0}b,i,{1})\n".format(curgraph['name'], curgraph.ecount() - 3 + len(curgraph['dots']) + numnewdots))
        else:
            newdots = ','.join(['i{}'.format(x) for x in range(extradotstart, curgraph.ecount() - 1 + len(curgraph['dots']) + numnewdots)])

            # the function should rename Z'name' to Z'name'b
            f.write("\t\t#call insertionS(M{},i{},i{}{},i{}{}{},{},{})\n".format(curgraph['name'],
                            cgi[curgraph['insertion'][0][0]] + 1,  # p1.p1
                            cgi[curgraph['insertion'][0][1]] + 1,  # p2.p2
                            ',' * (len(olddots) > 0) + olddots,
                            cgi[curgraph['insertion'][0][0]] + 1,  # output: p1.p1
                            ',' * (len(prop) > 0) + ','.join(prop),  # propagators
                            ',' * (len(newdots) > 0) + newdots,
                            curgraph.es[curgraph['insertion'][0][0]]['nonint'],  # epsilons
                            curgraph.es[curgraph['insertion'][0][1]]['nonint']))

            if cgi[curgraph['insertion'][0][1]] + 1 < curgraph.ecount() - 2 + len(curgraph['dots']) + numnewdots:
                f.write("\t\tif (count(M{}b,1)) Multiply replace_(<i{},i{}>,...,<i{},i{}>);\n".format(
                curgraph['name'],
                cgi[curgraph['insertion'][0][1]] + 2, cgi[curgraph['insertion'][0][1]] + 1,
                curgraph.ecount() - 2 + len(curgraph['dots']) + numnewdots,
                curgraph.ecount() - 2 + len(curgraph['dots']) + numnewdots - 1))

        f.write("\t\t#call Insertion{}()\n".format(curgraph['name']))

    if curgraph['triangle'] or curgraph['diamond'] or curgraph['insertion']:
        f.write("#endif\n")

    if curgraph.ecount() > 1 and not curgraph['insertion']:
        f.write("\t\t#call Prepare{}()\n".format(curgraph['name']))


def generateInsertion(metagraph, c, curgraph, f):
    """Write FORM procedure for cleaning up insertion integrals."""
    assert curgraph['insertion']
    assert curgraph.ecount() > 1

    e = next(e for e in metagraph.es.select(_source=c) if e['transgraph']['did_insertion'])

    f.write('#procedure Insertion{}()\n'.format(curgraph['name']))

    n = e.tuple[1]
    transgraph = e['transgraph']  # transitional graph

    # it could happen that at this point, we have a repetition of dot products due to the extra
    # dots the insertion generates. We map them here.
    dotmapd, dotmapc, dotmapv = {}, curgraph['transgraphs'][0].ecount(), []
    for d in curgraph['transgraphs'][-1]['dots']:
        if d not in dotmapd:
            dotmapd[d] = dotmapc
            dotmapc += 1
        dotmapv.append(dotmapd[d])

    # merge similar edges and dots
    if curgraph['merge'] != [i for i in range(curgraph['transgraphs'][0].ecount())] or dotmapc != len(dotmapv):
        # get the signs
        signmap = [0 for _ in range(curgraph['transgraphs'][0].ecount())] + \
            [1 if curgraph['transgraphs'][1]['dotmap'][d][1][0][0] == -1 else 0 for d in curgraph['transgraphs'][0]['dots']]

        f.write('* merge propagators with same momentum, and repeated dots\n')
        f.write('\tif (count(M{}b,1)) '.format(curgraph['name']))
        f.write(mergeLines(curgraph['transgraphs'][0], curgraph['transgraphs'][1],
            curgraph['merge'] + dotmapv, signmap) + '\n')

    # remove double dots from list
    curgraph['transgraphs'][-1]['dots'] = [x for i, x in enumerate(curgraph['transgraphs'][-1]['dots']) if
        x not in curgraph['transgraphs'][-1]['dots'][:i]]

    # write the dot product mapping of all the dot products that are removed in the insertion
    f.write('* rewrite insertion dots\n')
    t = transgraph.copy()
    if 'old_dots' in transgraph.attributes():
        t['dots'] = transgraph['old_dots']  # fix for the dots

    # dr = str(curgraph['transgraphs'][-1]['dotrewrite'])  # never used

    t['name'] = curgraph['name'] + 'd'
    f.write(generateDotmap(curgraph['transgraphs'][-1]['dotrewrite'], curgraph['transgraphs'][-1],
        t, [i for i in range(transgraph.ecount())]))

    gtoh = [i for i in range(curgraph.ecount())]  # identity mapping

    f.write(generateDotmap(e['dotmap'], t, transgraph,
        gtoh, metagraph.vs[n]['graph'], e['mapping']))

    f.write('.sort:end insertion;\n')
    f.write("#endprocedure\n\n")


# TODO: merge with mergeLines. the edgemap could map to the same edges: this means we sum them up
def generateDotmap(dotmap, oldgraph, newgraph, edgemap, isograph=None, isovm=None):
    """Generates for code for a dot product mapping and index mapping (from old to new).
    """
    # we cannot map if there are no lines
    if newgraph.ecount() == 1:
        if isograph is not None:
            return "\tid M{0} = M{1};\n\t#redefine HASTOPO{1} \"1\"\n".format(oldgraph['name'], isograph['name'])
        else:
            return "\tid M{0} = M{1};\n".format(oldgraph['name'], newgraph['name'])

    rewritestr = ''
    if isograph is not None:
        rewritestr += '* map dots to {} and change topology to {}\n'.format(
            ','.join('p{}.p{}'.format(*e) for e in isograph['dots']),
            isograph['name']
            )

    # Mappings to remove the external lines.
    oldnoext = mapAfterRemoval(oldgraph.ecount(), [e.index for e in oldgraph.es
                                                   if 0 in e.tuple or oldgraph.vcount() - 1 in e.tuple])
    newnoext = mapAfterRemoval(newgraph.ecount(), [e.index for e in newgraph.es
                                                   if 0 in e.tuple or newgraph.vcount() - 1 in e.tuple])

    nmap = [oldnoext[edgemap.index(newnoext.index(i))] + 1 for i in range(newgraph.ecount() - 2)] + \
            [i + newgraph.ecount() - 1 for i in range(len(newgraph['dots']))]

    coeff = ''
    isotonew = [x for x in range(oldgraph.ecount() - 2 + len(oldgraph['dots']))]
    if isograph is not None:
        coeff, isotonew = createEdgeDotMap(newgraph, isograph, isovm)
        nmap = [nmap[i] for i in isotonew]

    # map dots from old to ISO
    symbs = symbols(['j' + str(i + 1) for i in range(newgraph.ecount() - 2 + len(newgraph['dots']))] + ["`Q2'"])
    strtot = []
    for d, dm in sorted(dotmap.items(), key=lambda x: len(x[1][1])):
        if d[0] == d[1]:  # old position is propagator
            oldid = oldnoext[d[0]] + 1
        else:
            oldid = oldgraph.ecount() - 1 + oldgraph['dots'].index(d)
        exp = 0
        for coeffd, dot in dm[1]:
            if dot[0] == dot[1]:  # propagator increase
                if 0 in newgraph.es[dot[0]].tuple or newgraph.vcount() - 1 in newgraph.es[dot[0]].tuple:
                    exp += coeffd * symbs[-1]
                else:
                    exp += coeffd * symbs[isotonew.index(newnoext[dot[0]])]
            else:
                exp += coeffd * symbs[isotonew.index(newgraph.ecount() - 2 + newgraph['dots'].index(dot))]
        strtot.append((oldid, exp))

    # filter out simple mappings such as i4 -> i5 , i14 -> -i13
    simplemaps = []
    for i in reversed(range(len(strtot))):
        if strtot[i][1].func == Mul or len(strtot[i][1].args) == 0:
            simplemaps.append('i{},{}'.format(strtot[i][0], strtot[i][1]))
            del strtot[i]

    # find the part that has a constant mapping
    cons = sorted(set(range(1, oldgraph.ecount() - 1)) - {x[0] for x in strtot})
    rewritestr += "\tif (count(M{},1)) Multiply replace_({}{});\n".format(oldgraph['name'],
        ','.join('i{0},j{1}'.format(c, nmap.index(c) + 1) for c in cons),
        ',' * (len(simplemaps) > 0) + ','.join(sorted(simplemaps)))

    for i, x in enumerate(sorted(strtot, key=lambda x: x[0])):
        # filter out Q2 for now, because it could be 0
        prop = list(sorted(symbs.index(s) + 1 for s in x[1].atoms(Symbol) if str(s) != '`Q2\''))
        rewritestr += "\t#call replace(i{0},{1},j,{2},M{3},dotmapM{3}R{4})\n".format(x[0], str(x[1]).replace(' ', ''),
            ','.join('j' + str(a) for a in prop),
            oldgraph['name'], i + 1)

    #rewritestr += "\t#call tadpoleFilter(M{},j)\n".format(oldgraph['name'])

    if isograph is not None:

        if len(coeff) > 0:
            dr = ','.join('j{0},{1}i{0}'.format(x, '-' * (x in coeff)) for x in range(newgraph.ecount() - 1,
                newgraph.ecount() - 1 + len(newgraph['dots'])))

            rewritestr += "\tif (count(M{0},1)) Multiply replace_(M{0},M{1},<j1,i1>,...,<j{2},i{2}>{3});\n".format(oldgraph['name'],
                isograph['name'],
                newgraph.ecount() - 2, ',' * (len(dr) > 0) + dr)

        else:
            rewritestr += "\tif (count(M{0},1)) Multiply replace_(M{0},M{2},<j1,i1>,...,<j{1},i{1}>);\n".format(
                oldgraph['name'], newgraph.ecount() - 2 + len(newgraph['dots']), isograph['name'])
        rewritestr += "\t#redefine HASTOPO{} \"1\"\n".format(isograph['name'])
    else:
        # WARNING: we do not do an edgemap here
        rewritestr += "\tif (count(M{0},1)) Multiply replace_(M{0},M{2},<j1,i1>,...,<j{1},i{1}>);\n".format(
            oldgraph['name'], newgraph.ecount() - 2 + len(newgraph['dots']), newgraph['name'])

    return rewritestr


def createEdgeDotMap(oldgraph, newgraph, vm):
    """Maps edges and dots from oldgraph to new graph, given a vertex mapping.
    """
    ngi = mapAfterRemoval(newgraph.ecount(), [e.index for e in newgraph.es
                                if 0 in e.tuple or newgraph.vcount() - 1 in e.tuple])

    # map the edges and dot products back from h to c, preserving the index
    ctot = mapEdge(newgraph, oldgraph, vm, not _old_format)
    cdott = [tuple(sorted((ctot[d[0]], ctot[d[1]]))) for d in newgraph['dots']]

    # find the correct sign. The sign of the vertex merge was absorbed in the dotmap earlier.
    signchange = []
    for i, di in enumerate(newgraph['dots']):
        d = (newgraph.es[di[0]].tuple, newgraph.es[di[1]].tuple)
        if (vm[d[0][0]] > vm[d[0][1]]) != (vm[d[1][0]] > vm[d[1][1]]):  # sign change
            signchange.append(oldgraph['dots'].index(cdott[i]) + oldgraph.ecount() - 1)

    # map the edges from reduced newgraph back to reduced oldgraph
    tgi = mapAfterRemoval(oldgraph.ecount(), [e.index for e in oldgraph.es
                                                if 0 in e.tuple or oldgraph.vcount() - 1 in e.tuple])
    ntog = [tgi[ctot[ngi.index(e)]] for e in range(oldgraph.ecount() - 2)]
    doti = [oldgraph.ecount() - 2 + oldgraph['dots'].index(i) for i in cdott]

    return signchange, ntog + doti


def mergeLines(oldgraph, newgraph, mergemap, signmap):
    """Writes FORM code that sums up indices of lines and dot products that will be merged.
    """
    def isext(g, i):
        return (i < g.ecount() and
                (0 in g.es[i].tuple or g.vcount() - 1 in g.es[i].tuple))

    mmi = mapAfterRemoval(oldgraph.ecount() + len(oldgraph['dots']), [e for e in range(oldgraph.ecount()) if isext(oldgraph, e)])
    nmi = mapAfterRemoval(newgraph.ecount() + len(newgraph['dots']), [e for e in range(newgraph.ecount()) if isext(newgraph, e)])

    mergemap = [[j for j in range(len(mergemap)) if mergemap[j] == i] for i in sorted(set(mergemap))]
    indices = [[('i' + str(mmi[i] + 1), '-' * signmap[i]) for i in e if not isext(oldgraph, i)] for e in mergemap]

    kk = [','.join(['{},{}i{}'.format(j, js, nmi[i] + 1) if not isext(newgraph, i) else j + ",{}`Q2'".format(js)
            for j, js in x if ('i' + str(nmi[i] + 1) != j) or (signmap[i] == 1)]) for i, x in enumerate(indices) if len(x) > 0]
    k = ','.join(z for z in kk if len(z) > 0)

    return "Multiply replace_(M{},M{}{});\n".format(oldgraph['name'], newgraph['name'],
                ',' * (len(k) > 0) + k)


def wrap_text(text, first_width, width=None, separator=' ', pre_newline='',
              post_newline=''):
    """Return the text wrapped by the given width.

    Parameters
    ----------
    text : str
        The text to be processed.
    first_width : int
        The line width for the first line.
    width : int, optional
        The line width for other lines.
    separator : str or sequence of str, optional
        The word(s) after which a line can break.
    pre_newline : str, optional
        The word inserted before every line break.
    post_newline : str, optional
        The word inserted after every line break.

    Returns
    -------
    str
        The modified text to fit in the line width.

    """
    if width is None:
        width = first_width

    if not hasattr(separator, '__contains__'):
        separator = [separator]
    separator = '|'.join(re.escape(str(s)) for s in separator)

    first_width -= len(pre_newline)
    if first_width <= 0:
        raise ValueError('text width of the first line is too small')

    width -= len(pre_newline) + len(post_newline)
    if width <= 0:
        raise ValueError('text width is too small')

    if len(text) <= first_width:
        return text

    lines = []
    words = re.compile(r'\s*({0})\s*'.format(separator)).split(text)
    if words:
        words.reverse()
        next_line = ''
        while words:
            next_word = words[-1]
            if (len(next_line) + len(next_word) <=
                    (width if lines else first_width)):
                next_line += words.pop()
            else:
                if next_line:
                    lines.append(next_line)
                    next_line = ''
                else:
                    raise ValueError('word "{0}" is too long'.format(
                        next_word))
        if next_line:
            lines.append(next_line)
    return (pre_newline + '\n' + post_newline).join(lines)


def generateAll(fromfile=False, loops=4):
    if not fromfile:
        d = createTopLevelDiagrams(loops) + createFixedTopologies(loops)
        m = createTree(d, order=('insertion', 'carpet', 'triangle', 'diamond'))

        with open('metagraph{}.pickle'.format(loops), 'wb') as f:
            pickle.dump(m, f, protocol=2)
    else:
        with open('metagraph{}.pickle'.format(loops), 'rb') as f:
            m = pickle.load(f)
            if _old_format:
                # Fix 'convert' and 'mindot' in old pickle files.
                for v in m.vs:
                    g = v['graph']
                    if g['triangle'] or g['diamond'] or g['insertion']:
                        for attr in ('convert', 'mindot'):
                            if attr in g.attributes():
                                del g[attr]

            m.write_dot(m['name'] + '.dot')

    # notation file should be generated from notation.py
    generateMomentumConversion([x['graph'] for x in m.vs])
#    createTadpoleFilter(m)
#    convertTreeToFORM(m)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--loops', default=4, type=int, help='Number of loops')
    parser.add_argument('--restore', action='store_true', default=False, help='Restore from pickle file')
    parser.add_argument('--old', action='store_true', default=False, help='Support old format')
    args = parser.parse_args()
    _old_format = args.old
    generateAll(args.restore, args.loops)
