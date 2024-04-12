#!/usr/bin/env python3

import argparse
import os

from fnmatch import fnmatch
from pathlib import Path
from subprocess import Popen, PIPE, TimeoutExpired

parser = argparse.ArgumentParser()
parser.add_argument('type', help='Type of experiment to run. (0 = sparse cut, 1 = low conductance cut, 2 = normalized cut, 3 = greedy normalized cut, 4 = dynamic programming normalized cut)')
parser.add_argument('filename', help='Filename or directory containing the graph(s) used.')
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--logtostderr', action='store_true', help='Log to standard error')
parser.add_argument('-p', '--phi', help='Parameter phi used in the experiment')
parser.add_argument('-t', '--potential', help='Threshold of the random walk to certify expansion.')
parser.add_argument('-c', '--chaco', action='store_true', help='Parse graphs in chaco format.')
parser.add_argument('-k', '--num_parts', help='Parameter for multicut problems.')
parser.add_argument('-n', '--num_tries', help='Number of times to run the experiments.')

args = parser.parse_args()

print(args)

MTX_PARSER_PATH = './experiment/graphs/parse_mtx_graph.py'
EXECUTABLE_PATH = './bazel-bin/main/rw'

# First we want to get all files, by looking for mtx files in the input directory
FILETYPE = '*.mtx'
BAD_WORDS = ['_coord', '_area', '_population', '_Gcoord', '_b.']

# Parse the experiment type to run and build the command to execute.
process_options = [ EXECUTABLE_PATH, '-hierarchy', '-writeToFile' ]

seeds = [
    2044263138,
    26474971,
    305290509,
    3840931955,
    2533745830,
    2960824372,
    3479911025,
    1020351391,
    110446324,
    103896342,
]

num_tries = 10
if args.num_tries:
    num_tries = int(args.num_tries)

if args.logtostderr:
    process_options.append('-v')
    process_options.append('1')
    process_options.append('-logtostderr')

if args.chaco:
    process_options.append('-chaco')
    FILETYPE = '*.chaco'

if args.phi:
    process_options.append('-phi')
    process_options.append(str(args.phi))

if args.potential:
    process_options.append('-potential')
    process_options.append(str(args.potential))

if args.num_parts:
    process_options.append('-numParts')
    process_options.append(str(args.num_parts))

process_options.append('-experimentType')
experiment_type = int(args.type)
experiment_name = ''
match experiment_type:
    case 0:
        process_options.append('0')
        experiment_name = 'sparse cut'
    case 1:
        process_options.append('1')
        experiment_name = 'low conductance cut'
    case 2:
        process_options.append('2')
        experiment_name = 'normalized cut'
    case 3:
        process_options.append('3')
        experiment_name = 'greedy normalized cut'
    case 4:
        process_options.append('4')
        experiment_name = 'dynamic programming normalized cut'
    case _:
        print('Please specify an experimental type in the range [0, 3].')
        exit(1)

TIMEOUT = 1200

# phis = [ 0.3, 0.2, 0.1 ] #, 0.05, 0.025 ] #, 0.01 ]
phis = [ 0.0 ]
# thresholds = [ 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001 ]
thresholds = [ 0.0001 ]

# find the files 
graph_files = []

if os.path.isfile(args.filename):
  graph_files = [ Path(args.filename) ]
else:
  for path, subdirs, files in os.walk(args.filename):
    for name in files:
      if fnmatch(name, FILETYPE) and all(bad_word not in name for bad_word in BAD_WORDS):
        graph_files.append(Path(os.path.join(path, name)))

print(f"Found {len(graph_files)} files")

def experiment_run_with_phi(phi, thresholds, num_tries):
    for potential in thresholds:
        print(f'phi: {phi}, potential: {potential}')
        for t in range(num_tries):
            print(f'Run {t + 1} of {num_tries}')
            main = Popen(process_options
                         + [path]
                         + ['-name', path.stem]
                         + ['-seed', str(seeds[t])]
                         + ['-phi', str(phi)])
                   #      + ['-potential', str(potential)], stdout=PIPE)
            try:
                main.wait(timeout=TIMEOUT)
                if main.returncode != 0:
                    print('Process terminated improperly')
                    print(main.returncode)
                if main.returncode == 42:
                    print("The graph hierarchy does not terminate with this phi. Skipping...")
            except TimeoutExpired:
                main.terminate()
                print(f'The graph was not decomposed in the alotted time of {TIMEOUT} seconds.')
                break

for path in graph_files:
    print(f'Finding {experiment_name} on {path.stem}')
    for phi in phis:
        experiment_run_with_phi(phi, thresholds, num_tries)
