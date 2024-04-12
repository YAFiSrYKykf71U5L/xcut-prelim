#!/usr/bin/env python3

import argparse
import os

from fnmatch import fnmatch
from pathlib import Path
from subprocess import Popen, PIPE, TimeoutExpired
import resource

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Filename or directory containing the graph(s) used.')
#parser.add_argument('output_directory', help = 'Directory the chaco file is going to be written to.')

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-k', '--num_parts', help='Parameter for multicut problems.')
parser.add_argument('-n', '--num_tries', help='Number of times to run the experiments.')

args = parser.parse_args()
print(args)

# if args.output_directory[-1] != '/':
#    args.output_directory = args.output_directory + '/'


#METIS_PATH = '/home/voetschm94cs/local/bin/gpmetis'
METIS_PATH = '../../local/bin/gpmetis'
if not os.path.exists(METIS_PATH):
    print("Could not find the metis executable, please correct the path")
    exit(1)

# might have to fix this if you want to run it
GRACLUS_PATH = '../Graclus/graclus'
GRACLUS_OPTIONS = [ "-l", "20", "-b" ]
if not os.path.exists(METIS_PATH):
    print("Could not find the Graclus executable, please correct the path")
    exit(1)

KAHIP_PATH = '../KaHIP/deploy/kaffpa'
KAHIP_OPTIONS = ["--preconfiguration=fsocial"]
if not os.path.exists(METIS_PATH):
    print("Could not find the KaHiP executable, please correct the path")
    exit(1)

EVAL_PATH = './bazel-bin/main/eval_ncut'
if not os.path.exists(METIS_PATH):
    print("Could not find the evaluation executable, please correct the path")
    exit(1)

print("Found all the executables. Proceeding...")

# First we want to get all files, by looking for chaco files in the input directory
FILETYPE = '*.chaco'
BAD_WORDS = ['_coord', '_area', '_population', '_Gcoord', '_b.']
graph_files = []

if os.path.isfile(args.filename):
  graph_files = [ Path(args.filename) ]
else:
  for path, subdirs, files in os.walk(args.filename):
    for name in files:
      if fnmatch(name, FILETYPE) and all(bad_word not in name for bad_word in BAD_WORDS):
        graph_files.append(Path(os.path.join(path, name)))

num_tries = 1
if args.num_tries != None:
    num_tries = int(args.num_tries)


TIMEOUT = 1200

parts = [2, 4, 8, 16, 32, 64, 128]

print("Parsed args.")

def run_tool(process_options, process_name):
    #outfile = open(args.output_directory + path.stem + f".{process_name}.log", "w")
    usage_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    main = Popen(process_options)#, stdout=outfile)

    try:
        main.wait(timeout=TIMEOUT)

        if main.returncode != 0:
            print('Process terminated improperly')
            print(main.returncode)
            exit(1)

    except TimeoutExpired:
        main.terminate()
        print(f'The graph was not decomposed in the alotted time of {TIMEOUT} seconds.')
        return

    usage_end = resource.getrusage(resource.RUSAGE_CHILDREN)

    return usage_end.ru_utime - usage_start.ru_utime 
    
    #outfile.close()

def run_eval(path, alg_name, partition_file_path, time):
    main = Popen([ EVAL_PATH, 
                 '-logtostderr', 
                 '-name', path.stem, 
                 '-algo', alg_name,
                 '-time', str(time),
                 #'-numParts', numParts,
                 path,
                 partition_file_path  ], )
    try:
        main.wait(timeout=TIMEOUT)

        if main.returncode != 0:
            print('Process terminated improperly')
            print(main.returncode)
            exit(1)
    except TimeoutExpired:
        main.terminate()
        print(f'The graph was not decomposed in the alotted time of {TIMEOUT} seconds.')
        return

def run_METIS(path, num_tries, num_parts):
    options =[METIS_PATH, path, str(num_parts)]
    for i in range(num_tries):
        print(f'Run {i + 1} of {num_tries}')
        time = run_tool(options, 'metis')
        print(time)
        partition_file = Path(str(path) + f'.part.{num_parts}')
        run_eval(path, 'metis', partition_file, time)
        partition_file.unlink()

def run_graclus(path, num_tries, num_parts):
    options = [GRACLUS_PATH, path, str(num_parts)] + GRACLUS_OPTIONS
    for i in range(num_tries):
        print(f'Run {i + 1} of {num_tries}')
        time = run_tool(options, 'graclus')
        print(time)
        partition_file = Path(path.stem + f'.chaco.part.{num_parts}')
        run_eval(path, 'graclus', partition_file, time)
        partition_file.unlink()

def run_KaHIP(path, num_tries, num_parts):
    options = [KAHIP_PATH, path, '--k', str(num_parts)] + KAHIP_OPTIONS
    for i in range(num_tries):
        print(f'Run {i + 1} of {num_tries}')
        time = run_tool(options, 'kahip')
        print(time)
        partition_file = Path(f'tmppartition{num_parts}')
        run_eval(path, 'kahip', partition_file, time)
        partition_file.unlink()

print("Getting started.")

for path in graph_files:
    for size in parts:
        print(f'Evaluating sparsifier quality on {path.stem}')
        run_METIS(path, num_tries, size)
        run_KaHIP(path, num_tries, size)
        run_graclus(path, num_tries, size)
