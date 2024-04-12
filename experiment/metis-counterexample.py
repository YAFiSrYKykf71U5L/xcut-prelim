#!/usr/bin/env python3

"""This file contains a few routines that (hopefully) generate examples on which METIS performs poorly."""

import argparse
import math
import sys
import numpy as np

def output(args, g, n, m):
    print(f'{n} {m // 2}')
    #print(f'{n} {m // 2} 1')
    for u, vs in g.items():
        for (v, w) in vs:
            print(f'{v+1} ', end='')
#            print(f'{v+1} {w} ', end='')
        print('')

# create a clique with hairs sticking out from each vertex
# this is a hard example for metis
def clique_with_fur(args):
    g = {}
    n = args.n
    w = args.w
    W = args.W
    m = 0

    # generate a clique
    for u in range(n):
        g[u] = []
        for v in range(n):
            if u == v:
                continue

            g[u].append((v, w))
            m += 1

    for u in range(n, 2 *n):
        g[u] = []

        g[u].append((u - n, W))
        g[u - n].append((u, W))
        m += 2

    output(args, g, 2 * n, m)

# same as above but use a margulis expander instead of a complete graph
def margulis(args):
    g = {}
    n = args.n
    w = args.w
    W = args.W
    m = 0

    def f(a, b):
        return a * n + b

    def k(a, b):
        return ((a + n) % n) * n + ((b + n) % n)

    def S(a, b):
        return k(a, a + b)

    def T(a, b):
        return k(a + b, b)

    def s(a, b):
        return k(a + 1, b)

    def t(a, b):
        return k(a, b + 1)

    def Sinv(a, b):
        return k(a, b - a)

    def Tinv(a, b):
        return k(a - b, b)

    def sinv(a, b):
        return k(a - 1, b)

    def tinv(a, b):
        return k(a, b - 1)


    for i in range(n):
        for j in range(n):
            g[f(i ,j)] = [
                (S(i, j), w),
                (T(i, j), w),
                (s(i, j), w),
                (t(i, j), w),
                (Sinv(i, j), w),
                (Tinv(i, j), w),
                (sinv(i, j), w),
                (tinv(i, j), w),
            ]
            m += 8

    for i in range(n):
        for j in range(n):
            g[(n + i) * n + j] = [(i * n + j, W)]
            g[i * n + j].append(((n+i) * n + j, W))
            m += 2

    output(args, g, 2 * n * n, m)

def pair_iter(list):
    i = 0
    while i < len(list):
        yield (list[i], list[(i+1) % len(list)])
        i += 1

# same as above but use a random regular graph as an expander
# the algorithm used to generate expanders follows the approach of the networkx library
def random_expander(args):
    g = {}
    n = args.n
    w = args.w
    W = args.W

    edges = set()

    perm = np.arange(n)

    print(f"Iterations: {(int(math.log(n)) // 2) + 1}", file=sys.stderr)

    for i in range((int(math.log(n)) // 2) + 1):
        print(i, file=sys.stderr)

        iterations = n

        while len(edges) != (i + 1) * n:
            iterations -= 1
            print(iterations, file = sys.stderr)

            np.random.shuffle(perm)

            new_edges = {
                (u, v)
                for u, v in pair_iter(perm) # zip(perm, np.concatenate((perm[1:], [perm[0]])))
                if (u, v) not in edges and (v, u) not in edges
            }

            if len(new_edges) == n:
                edges.update(new_edges)
            
            if iterations == 0:
                print("No more iterations left. Exiting.", file=sys.stderr)
                exit(1)

    for i in range(n):
        g[i] = []

    m = 0

    for (u, v) in edges:
        g[u].append((v, w))
        g[v].append((u, w))
        m += 2

    for u in range(n, 2 *n):
        g[u] = []

        g[u].append((u - n, W))
        g[u - n].append((u, W))
        m += 2

    output(args, g, 2 * n, m)

def silly(args):
    g = {}
    n = args.n
    w = args.w
    W = args.W
    m = 0

    for u in range(n):
        g[u] = []

    # generate a clique
    for u in range(n):
        while len(g[u]) < math.log(n):
            v = np.random.randint(n)

            if u == v: 
                continue
                
            g[u].append((v, w))
            g[v].append((u, w))

            m+= 2

    for u in range(n, 2 *n):
        g[u] = []

        g[u].append((u - n, W))
        g[u - n].append((u, W))
        m += 2

    output(args, g, 2 * n, m)

parser = argparse.ArgumentParser(description = 'METIS counterexample utility')
subparsers = parser.add_subparsers()

fur_parser = subparsers.add_parser('fur', help='clique of size n with "fur"')
fur_parser.set_defaults(func=clique_with_fur)
fur_parser.add_argument('-n', type=int, default = 100, help='number of vertices in the clique')
fur_parser.add_argument('-w', type=int, default = 1, help='weight of inner edges')
fur_parser.add_argument('-W', type=int, default = 2, help='weight of spike edges')

margulis_parser = subparsers.add_parser('margulis', help='margulis expander of size n^2 with "fur"')
margulis_parser.set_defaults(func=margulis)
margulis_parser.add_argument('-n', type=int, default = 10, help='number of vertices in the clique')
margulis_parser.add_argument('-w', type=int, default = 1, help='weight of inner edges')
margulis_parser.add_argument('-W', type=int, default = 2, help='weight of spike edges')

random_parser = subparsers.add_parser('random', help='random regular graph of size n with "fur"')
random_parser.set_defaults(func=random_expander)
random_parser.add_argument('-n', type=int, default = 100, help='number of vertices in the clique')
random_parser.add_argument('-w', type=int, default = 100, help='weight of inner edges')
random_parser.add_argument('-W', type=int, default = 101, help='weight of spike edges')

silly_parser = subparsers.add_parser('silly', help='silly regular graph of size n with "fur"')
silly_parser.set_defaults(func=silly)
silly_parser.add_argument('-n', type=int, default = 100, help='number of vertices in the clique')
silly_parser.add_argument('-w', type=int, default = 100, help='weight of inner edges')
silly_parser.add_argument('-W', type=int, default = 101, help='weight of spike edges')

args = parser.parse_args()
args.func(args)
