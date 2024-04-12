#ifndef GRAPH_IO_H_
#define GRAPH_IO_H_

#include <glog/logging.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "lib/datastructures/array_graph.hpp"
#include "lib/datastructures/dynamic_graph.hpp"
#include "lib/datastructures/sparsifier.hpp"

namespace fs = std::filesystem;

/**
 * Read graph from standard in.
 */
ArrayGraph::Graph *readGraphMTX();

/**
 * Read graph from the specified file.
 */
ArrayGraph::Graph *readGraphMTX(fs::path const &filename);

/**
 * Read graph from standard in.
 */
DynamicGraph::Graph *readDynamicGraphMTX();

/**
 * Read graph from the specified file.
 */
DynamicGraph::Graph *readDynamicGraphMTX(fs::path const &filename);

ArrayGraph::Graph *readGraphChaco();
ArrayGraph::Graph *readGraphChaco(fs::path const &filename);

void writeGraphMTX(ArrayGraph::Graph *graph, fs::path const &filename);
void writeGraphChaco(ArrayGraph::Graph *graph, fs::path const &filename);

template <typename W>
CutSparsifier::CutSparsifier<W> readSparsifier(fs::path const &filename) {
  LOG(INFO) << "Reading sparsifier from file " << filename;

  std::ifstream f{filename};
  if (!f) {
    LOG(ERROR) << "Could not read file " << filename;
    throw 1;
  }

  // read parent and weight arrays
  std::string line;
  std::vector<std::vector<int>> parents;
  std::vector<std::vector<W>> weights;

  while (std::getline(f, line)) {
    // create new layer
    parents.emplace_back();
    weights.emplace_back();

    std::stringstream ss{line};

    int p;
    W w;
    ss >> p >> w;

    parents.back().push_back(p);
    weights.back().push_back(w);
  }

  // construct the sparsifier
}

template <typename W>
void writeSparsifier(CutSparsifier::CutSparsifier<W> const &sparsifier,
                     fs::path const &filename) {
  std::ofstream f{filename};

  for (int i = 0; i < sparsifier.size(); i++) {
    for (int j = 0; j < sparsifier.getLayerSize(i); j++) {
      auto [p, w] = sparsifier.getParentEdge(i, j);
      f << p << " " << w << " ";
    }
    f << std::endl;
  }
}

std::vector<int> readPartition(fs::path const &filename);
void writePartition(std::vector<int> const &cut, fs::path const &filename);

/**
 * Read Partition and set the cut array.
 */
void readPartition(ArrayGraph::Graph *graph, fs::path const &filename);

void writePartition(ArrayGraph::Graph *graph, fs::path const &filename);

#endif  // GRAPH_IO_H_
