import argparse
import os
import numpy as np

from fnmatch import fnmatch
from pathlib import Path

def process_chaco_file(file_path, output_file):
    # Initialize variables to store statistics
    degrees = []
    
    # Read the file
    with open(file_path, 'r') as file:
        # Read header
        header = file.readline().split()
        num_vertices, num_edges = map(int, header)

        # Process lines
        for line in file:
            # Split adjacency list
            adjacency = list(map(int, line.split()))

            # Add degree of vertex to degrees list
            degree = len(adjacency)
            degrees.append(degree)

    # Write statistics to the outfile
    if degrees:
        max_degree = max(degrees)
        avg_degree = np.mean(degrees)
        variance = np.std(degrees)
        percentile_values = list(map(int, np.percentile(degrees, [25, 50, 75, 90])))
        
        # Write the statistics to the file
        output_file.write(
            f'{file_path.stem},{num_vertices},{num_edges},{max_degree},{avg_degree:.2f},{percentile_values[0]},{percentile_values[1]},{percentile_values[2]},{percentile_values[3]},{variance}\n'
        )
    else:
        print("No lines found in the file.")


# Create the command line parser
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Filename or directory containing the graph(s) used.')

# Parse Arguments
args = parser.parse_args()

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

# Create the csv file for writing outputs
outfile_name = 'graph-statistics.csv'

if os.path.exists(outfile_name):
   append_write = 'a'
   header = ''
else:
   append_write = 'w'
   header = 'Graph Name,Vertices,Edges,Max Degree,Avg Degree,25th,50th,75th,90th,STD\n'

# Open the file and write the header if applicable
outfile = open(outfile_name, append_write)
outfile.write(header)

for graph_file in graph_files:
   print(f'Processing graph {graph_file.stem}')
   process_chaco_file(graph_file, output_file=outfile)
