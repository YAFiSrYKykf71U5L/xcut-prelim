#include "graph_io.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

/**
 * Write a vector to a file, separated by newlines. These are not exported.
 */
template <typename T>
void writeVector(std::vector<T> const &v, fs::path const &filename) {
  std::ofstream f(filename);

  for (auto &e : v) {
    f << e << std::endl;
  }
}

template <typename T>
std::vector<T> readVector(fs::path const &filename) {
  std::string line;
  std::vector<T> vec;

  std::ifstream f(filename);
  if (!f) {
    LOG(ERROR) << "Error opening file " << filename;
    return {};
  }

  while (std::getline(f, line)) {
    if (line[0] == '%' || line.empty()) continue;

    std::stringstream ss(line);
    T val;

    ss >> val;

    vec.push_back(val);
  }

  return vec;
}

/**
 * Read an unweighted graph from a matrix-market file.
 */
ArrayGraph::Graph *readGraphMTX(std::istream &f) {
  std::string line;
  // get first non-comment line
  while (std::getline(f, line) && line[0] == '%') continue;

  int n, m;
  std::stringstream ss(line);
  ss >> n >> n >> m;

  VLOG(2) << "Graph file has " << n << " vertices and " << m << " non-zeros.";

  std::vector<ArrayGraph::Edge> es;

  int u, v;

  int i = 0;
  while (std::getline(f, line)) {
    if (line[0] == '%') continue;

    std::stringstream ss(line);
    ss >> u >> v;

    LOG_EVERY_N(INFO, 100000) << "Read " << i++ * 100000 << " edges.";

    if (u == v) {
      VLOG(4) << "Edge is a self-loop. Skipping. This may be implemented at a "
                 "later time.";
      continue;
    }

    if (v < u) {
      std::swap(u, v);
    }

    es.emplace_back(u - 1, v - 1);
  }

  std::sort(es.begin(), es.end(),
            [](ArrayGraph::Edge const &e, ArrayGraph::Edge const &f) {
              return e.from < f.from || (e.from == f.from && e.to < f.to);
            });

  auto last =
      std::unique(es.begin(), es.end(),
                  [](ArrayGraph::Edge const &e, ArrayGraph::Edge const &f) {
                    return e.from == f.from && e.to == f.to;
                  });
  es.erase(last, es.end());

  VLOG(0) << "Finished reading graph.";

  return new ArrayGraph::Graph(n, es);
}

ArrayGraph::Graph *readGraphMTX(fs::path const &filename) {
  VLOG(0) << "Reading graph from file " << filename;

  std::ifstream f(filename);
  if (!f) {
    VLOG(0) << "Could not read file " << filename;
    return new ArrayGraph::Graph();
  }

  return readGraphMTX(f);
}

ArrayGraph::Graph *readGraphMTX() {
  LOG(INFO) << "Reading graph from standard input.";
  return readGraphMTX(std::cin);
}

/**
 * Read an unweighted graph from a matrix-market file.
 */
DynamicGraph::Graph *readDynamicGraphMTX(std::istream &f) {
  std::string line;
  // get first non-comment line
  while (std::getline(f, line) && line[0] == '%') continue;

  int n, m;
  std::stringstream ss(line);
  ss >> n >> n >> m;

  VLOG(2) << "Graph file has " << n << " vertices and " << m << " non-zeros.";

  std::vector<DynamicGraph::Edge> es;

  int u, v;

  int i = 1;
  while (std::getline(f, line)) {
    if (line[0] == '%') continue;

    std::stringstream ss(line);
    ss >> u >> v;

    LOG_EVERY_N(INFO, 100000) << "Read " << i++ * 1000000 << " edges.";

    if (u == v) {
      VLOG(4) << "Edge is a self-loop. Skipping. This may be implemented at a "
                 "later time.";
      continue;
    }

    if (v < u) {
      std::swap(u, v);
    }

    es.emplace_back(u - 1, v - 1);
  }

  std::sort(es.begin(), es.end(),
            [](DynamicGraph::Edge const &e, DynamicGraph::Edge const &f) {
              return e.from < f.from || (e.from == f.from && e.to < f.to);
            });

  auto last =
      std::unique(es.begin(), es.end(),
                  [](DynamicGraph::Edge const &e, DynamicGraph::Edge const &f) {
                    return e.from == f.from && e.to == f.to;
                  });
  es.erase(last, es.end());

  VLOG(0) << "Finished reading graph.";

  return new DynamicGraph::Graph(n, es);
}

DynamicGraph::Graph *readDynamicGraphMTX(fs::path const &filename) {
  if (filename.extension() != ".mtx") {
    LOG(WARNING)
        << "The file extension does not match .mtx, you might be reading "
           "an incompatible file.";
  }

  VLOG(0) << "Reading graph from file " << filename;

  std::ifstream f(filename);
  if (!f) {
    VLOG(0) << "Could not read file " << filename;
    return new DynamicGraph::Graph();
  }

  return readDynamicGraphMTX(f);
}

DynamicGraph::Graph *readDynamicGraphMTX() {
  LOG(INFO) << "Reading graph from standard input.";
  return readDynamicGraphMTX(std::cin);
}

/**
 * Write the graph to matrix market format.
 */
void writeGraphMTX(ArrayGraph::Graph *graph, fs::path const &filename) {
  std::ofstream f(filename);

  f << graph->size() << " " << graph->size() << " " << graph->getVolume()
    << std::endl;

  for (int u = 0; u < graph->size(); u++) {
    for (auto e = graph->cbeginEdges(u); e != graph->cendEdges(u); e++) {
      f << u + 1 << " " << *e + 1 << std::endl;
    }
  }
}

/**
 * Read an unweighted graph from a chaco/metis file.
 */
ArrayGraph::Graph *readGraphChaco(std::istream &f) {
  std::string line;

  // get first non-comment line
  while (std::getline(f, line) && line[0] == '%') continue;

  int n, m;
  std::stringstream ss(line);
  ss >> n >> m;

  std::vector<int> edges(m * 2);
  std::vector<int> degrees(n);
  std::vector<int> start(n + 1, 0);

  int i = 0;
  int u;
  int node = 0;
  while (std::getline(f, line)) {
    if (line[0] == '%') continue;

    std::stringstream ss(line);

    start[node] = i;

    while (ss >> u) {
      edges[i] = u - 1;
      i++;
    }

    degrees[node] = i - start[node];

    node++;
  }

  start[n] = 2 * m;

  VLOG(0) << "Finished reading graph.";

  return new ArrayGraph::Graph(n, m, std::move(edges), std::move(degrees),
                               std::move(start));
}

ArrayGraph::Graph *readGraphChaco(fs::path const &filename) {
  VLOG(0) << "Reading graph from file " << filename;

  std::ifstream f(filename);
  if (!f) {
    VLOG(0) << "Could not read file " << filename;
    return new ArrayGraph::Graph();
  }

  return readGraphChaco(f);
}

ArrayGraph::Graph *readGraphChaco() {
  LOG(INFO) << "Reading graph from standard input.";
  return readGraphChaco(std::cin);
}

/**
 * Write the the graph to chaco/metis file format.
 */
void writeGraphChaco(ArrayGraph::Graph *graph, fs::path const &filename) {
  std::ofstream f(filename);

  f << graph->size() << " " << graph->getNumEdges() << std::endl;

  for (int u = 0; u < graph->size(); u++) {
    for (auto e = graph->cbeginEdges(u); e != graph->cendEdges(u); e++) {
      f << *e + 1 << " ";
    }
    f << std::endl;
  }
}

std::vector<int> readPartition(fs::path const &filename) {
  return readVector<int>(filename);
}

void writePartition(std::vector<int> const &cut, fs::path const &filename) {
  writeVector(cut, filename);
}

/**
 * Read a partition vector.
 */
void readPartition(ArrayGraph::Graph *graph, fs::path const &filename) {
  auto cut = readVector<int>(filename);

  int n = 0;

  for (auto i : cut) {
    n = std::max(n, i);
  }

  graph->setCutArray(n + 1, std::move(cut));
}

/**
 * Write a partition vector.
 */
void writePartition(ArrayGraph::Graph *graph, fs::path const &filename) {
  writeVector(graph->getCutArray(), filename);
}
