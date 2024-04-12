#!/usr/bin/env python3

import argparse
import os

from fnmatch import fnmatch
from pathlib import Path
from subprocess import Popen, PIPE, TimeoutExpired

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Filename or directory containing the graph(s) used.')
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--logtostderr', action='store_true', help='Log to standard error')
parser.add_argument('-p', '--phi', help='Parameter phi used in the experiment')
parser.add_argument('-t', '--potential', help='Threshold of the random walk to certify expansion.')
parser.add_argument('-c', '--chaco', help='Parse graphs in chaco format.')
parser.add_argument('-k', '--num_parts', help='Parameter for multicut problems.')
parser.add_argument('-n', '--num_tries', help='Number of times to run the experiments.')
parser.add_argument('--subsetType', help='Whether to sample subsets greedily or randomly.')

args = parser.parse_args()
print(args)

EXECUTABLE_PATH = './bazel-bin/main/quality'

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

# Parse the experiment type to run and build the command to execute.
process_options = [ EXECUTABLE_PATH, '-writeToFile', ]

num_tries = 5
if args.num_tries:
    num_tries = args.num_tries

if args.logtostderr:
    process_options.append('-v')
    process_options.append('1')
    process_options.append('-logtostderr')

if args.chaco:
    process_options.append('-chaco')

if args.phi:
    process_options.append('-phi')
    process_options.append(str(args.phi))

if args.potential:
    process_options.append('-potential')
    process_options.append(str(args.potential))

if args.subsetType:
    process_options.append('-subsetType')
    process_options.append(str(args.subsetType))

TIMEOUT = 300

phis = [ 0.4, 0.3, 0.2, 0.1, 0.05, 0.025, 0.01 ]
thresholds = [ 0.00001 ]

def experiment_run_with_phi(path, phi, thresholds, num_tries):
    for potential in thresholds:
        print(f'phi: {phi}, potential: {potential}')
        for t in range(num_tries):
            print(f'Run {t + 1} of {num_tries}')
            main = Popen(process_options
                         + [ path ]
                         + ['-name', path.stem]
                         + ['-phi', str(phi)]
                         + ['-potential', str(potential)], stdout=PIPE)
            try:
                main.wait(timeout=TIMEOUT)
                if main.returncode != 0:
                    print('Process terminated improperly')
                    print(main.returncode)
                if main.returncode == 42:
                    print("The graph hierarchy does not terminate with this phi. Skipping...")
                    break
            except TimeoutExpired:
                main.terminate()
                print(f'The graph was not decomposed in the alotted time of {TIMEOUT} seconds.')
                return

for path in graph_files:
    print(f'Evaluating sparsifier quality on {path.stem}')
    for phi in phis:
        experiment_run_with_phi(path, phi, thresholds, num_tries)
