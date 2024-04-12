#!/usr/bin/env python3

"""Given a sparse matrix in mtx format, convert it to the format used by metis.
"""

import argparse
import os
import sys

from fnmatch import fnmatch
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument('filename', help = 'Filename or directory of the file(s) that need(s) to be converted.')
parser.add_argument('output_directory', help = 'Directory the chaco file is going to be written to.')

args = parser.parse_args()

if args.output_directory[-1] != '/':
   args.output_directory = args.output_directory + '/'

# First we want to get all files, by looking for mtx files in the input directory
FILETYPE = '*.mtx'
BAD_WORDS = ['_coord', '_area', '_population', '_Gcoord', '_b.']
graph_files = []

if os.path.isfile(args.filename):
  graph_files = [ Path(args.filename) ]
else:
  for path, subdirs, files in os.walk(args.filename):
    for name in files:
      if fnmatch(name, FILETYPE) and all(bad_word not in name for bad_word in BAD_WORDS):
        graph_files.append(Path(os.path.join(path, name)))


for filename in graph_files:
    lines = []

    print("Converting " + filename.stem)

    infile = open(filename, "r")
    outfile = open(args.output_directory +  filename.stem + ".chaco", "w")

    print("Reading graph.")

    found_header = False
    n = 0
    m = 0

    edges = []

    for line in infile:
        line = line.strip()
        if line[0] == "%" or len(line) == 0:
            continue
        if not found_header:
           print("Found header")
           found_header = True
           n = int(line.split()[0])
           edges = [set() for i in range(n)]

        u, v, *_ = list(map(int, line.strip().split()[0:2]))

        if u == v or v in edges[u-1]:
            continue

        m += 1

        edges[u - 1].add(v)
        edges[v - 1].add(u)


    print("Writing to file " + args.output_directory + filename.stem + ".chaco")

    outfile.write(f"{n} {m}")
    for adjacency_list in edges:
        outfile.write('\n')
        outfile.write(' '.join(str(v) for v in adjacency_list))

    infile.close()
    outfile.close()
