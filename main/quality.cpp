// These first two are to get it to build with bazel on Arch, may be redundant
// on other systems.
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <glog/stl_logging.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "lib/algorithms/ford_fulkerson.hpp"
#include "lib/algorithms/push_relabel.hpp"
#include "lib/io/graph_io.hpp"
#include "lib/random_walk.hpp"
#include "util.hpp"

DEFINE_uint32(seed, 0,
              "Seed randomness with any positive integer. Default value '0' "
              "means a random seed will be chosen based on system time.");
DEFINE_uint32(numRepetitions, 50, "Number of times to repeat the experiment.");
DEFINE_uint32(subsetType, 0,
              "Choose what algorithm is used to choose the subset for the "
              "quality experiment.");
DEFINE_uint32(subsetSize, 0, "Size of the sampled subsets.");
DEFINE_bool(edmondskarp, false,
            "Use Edmonds-Karp algorithm instead of Push-Relabel.");
DEFINE_double(
    phi, 0.2,
    "Value of \\phi such that expansion of each cluster is at least \\phi");
DEFINE_uint32(boost, 1,
              "Value of the boundary-linkedness of each cluster. In other "
              "words, this is how much more we weigh outgoing edges versus "
              "internal ones.");
DEFINE_double(potential, 0.000001,
              "Potential of the random walk at which we accept convergence and "
              "certify a component is an expander.");
DEFINE_string(name, "unknown", "Name of the graph that is being processed.");
DEFINE_bool(stdin, false, "Read the graph from standard input.");
DEFINE_bool(chaco, false,
            "Input graph is given in the Chaco graph file format");
DEFINE_bool(writeToFile, false,
            "Write the output to a csv file, append if it already exists.");

enum SamplingAlgo {
  Random,
  Greedy,
};

std::vector<int> defaultSizes{10, 100, 1000, 10000, 100000};

class RandomSubset {
 private:
  std::vector<int> perm;
  std::mt19937 randomSource;

 public:
  RandomSubset(int n) : perm(n, 0) {
    std::random_device rd;
    randomSource = std::mt19937(rd());
    std::iota(perm.begin(), perm.end(), 0);
  }

  std::pair<std::vector<int>, std::vector<int>> randomSubset(int size) {
    // ensure that ss is the smaller of the two sets
    if (static_cast<size_t>(size) > perm.size() / 2) {
      size = perm.size() - size;
    }

    std::shuffle(perm.begin(), perm.end(), randomSource);

    auto mid = perm.begin() + size;

    std::vector<int> ss(perm.begin(), mid), ts(mid, perm.end());

    return std::make_pair(ss, ts);
  }
};

class GreedySubset {
 private:
  DynamicGraph::Graph *graph;
  std::mt19937 randomSource;
  std::uniform_int_distribution<> dist;

 public:
  GreedySubset(DynamicGraph::Graph *graph)
      : graph(graph), dist(0, graph->size() - 1) {
    std::random_device rd;
    randomSource = std::mt19937(rd());
  }

  std::pair<std::vector<int>, std::vector<int>> randomSubset(int size) {
    // ensure that ss is the smaller of the two sets
    if (size > graph->size() / 2) {
      size = graph->size() - size;
    }

    std::vector<int> ss;
    std::vector<bool> visited(graph->size(), false);
    std::queue<int> queue;

    auto start = dist(randomSource);

    queue.push(start);
    visited[start] = true;

    int i = 0;

    while (!queue.empty() && ss.size() < static_cast<size_t>(size)) {
      auto u = queue.front();
      queue.pop();

      ss.push_back(u);
      i++;

      for (auto i = 0; i < graph->getNumNeighbors(u); i++) {
        auto v = graph->getEdge(u, i);

        if (!visited[v]) {
          queue.push(v);
          visited[v] = true;
        }
      }
    }

    std::sort(ss.begin(), ss.end());

    std::vector<int> ts;
    i = 0;

    for (int u = 0; u < graph->size(); u++) {
      if (i >= size || u < ss[i]) {
        ts.emplace_back(u);
      } else {
        i++;
      }
    }

    return std::make_pair(ss, ts);
  }
};

template <typename G>
int minCutPartition(std::shared_ptr<G> graph, std::vector<int> const &ss) {
  std::vector<int> cut(graph->size(), 0);

  for (auto s : ss) {
    cut[s] = 1;
  }

  int minCut = 0;

  for (auto u = 0; u < graph->size(); u++) {
    for (auto it = graph->cbeginEdges(u); it != graph->cendEdges(u); it++) {
      auto v = *it;

      if (cut[u] != cut[v]) {
        minCut++;
      }
    }
  }

  // edges are counted twice
  return minCut / 2;
}

template <typename G, typename S, typename Sampler>
void qualityComparison(G *g, std::shared_ptr<S> sparsifier, Sampler &sampler,
                       std::ofstream &file, unsigned long time) {
  // auto prg = PushRelabel(g);
  // auto prs = PushRelabel(sparsifier);
  // auto prs = FordFulkerson(s);

  if (FLAGS_subsetSize != 0) {
    if (FLAGS_subsetSize >= static_cast<size_t>(g->size())) {
      LOG(ERROR) << "Can not sample more than n nodes.";
    }

    for (size_t i = 0; i < FLAGS_numRepetitions; i++) {
      auto [ss, ts] = sampler.randomSubset(FLAGS_subsetSize);

      // auto minCutPR = prg.maxFlow(ss, ts);
      auto minCut = minCutPartition(g, ss);

      auto [s, t] = sparsifier->addSourceSink(ss, ts);
      auto prs = PushRelabel(sparsifier);
      auto sparsifierMinCut = prs.maxFlow(s, t);
      sparsifier->removeSourceSink();

      LOG(INFO) << "Cut values: " << sparsifierMinCut << ", " << minCut;

      file << FLAGS_name << ",";
      file << FLAGS_phi << ",";
      file << FLAGS_potential << ",";
      file << time << ",";
      file << FLAGS_subsetSize << ",";
      file << sparsifierMinCut << ",";
      file << minCut << ",";
      file << (double)sparsifierMinCut / (double)minCut << std::endl;
    }
  } else {
    for (auto size : defaultSizes) {
      if (size >= g->size()) {
        continue;
      }

      for (size_t i = 0; i < FLAGS_numRepetitions; i++) {
        LOG(INFO) << size << ": " << i;
        auto [ss, ts] = sampler.randomSubset(size);

        auto minCut = minCutPartition(g, ss);

        auto [s, t] = sparsifier->addSourceSink(ss, ts);
        auto prs = PushRelabel(sparsifier);
        auto sparsifierMinCut = prs.maxFlow(s, t);
        sparsifier->removeSourceSink();

        LOG(INFO) << "Cut values: " << sparsifierMinCut << ", " << minCut;

        file << FLAGS_name << ",";
        file << FLAGS_phi << ",";
        file << FLAGS_potential << ",";
        file << time << ",";
        file << size << ",";
        file << sparsifierMinCut << ",";
        file << minCut << ",";
        file << (double)sparsifierMinCut / (double)minCut << std::endl;
      }
    }
  }
}

/**
 * TODO: fix push relabel
 */

template <typename G, typename S, typename Sampler>
void qualityComparison(G *g, std::shared_ptr<S> sparsifier, Sampler &sampler) {
  //   auto prg = PushRelabel(g);

  //   auto prgf = FordFulkerson(g);
  //   auto prf = FordFulkerson(sparsifier);

  if (FLAGS_subsetSize != 0) {
    if (FLAGS_subsetSize >= g->size()) {
      LOG(ERROR) << "Can not sample more than n nodes.";
    }

    for (int i = 0; i < FLAGS_numRepetitions; i++) {
      LOG(INFO) << "Sampling.";
      auto [ss, ts] = sampler.randomSubset(FLAGS_subsetSize);
      LOG(INFO) << "Sampling done.";

      auto minCut = minCutPartition(g, ss);
      // auto minCutPR = prg.maxFlow(ss, ts);
      // auto minCutFF = prgf.maxFlow(ss, ts);

      VLOG(0) << "min cut done. " << minCut;

      //   auto prf = FordFulkerson(sparsifier);
      //   auto sparsifierMinCut_ = prf.maxFlow(ss, ts);
      //   VLOG(0) << "ff sparsifire done";

      auto [s, t] = sparsifier->addSourceSink(ss, ts);
      auto prs = PushRelabel(sparsifier);
      auto sparsifierMinCut = prs.maxFlow(s, t);
      sparsifier->removeSourceSink();

      VLOG(0) << "sparsifier done.";

      LOG(INFO) << "Cut value: " << minCut;
      // LOG(INFO) << "Cut value (PR): " << minCutPR;
      // LOG(INFO) << "Cut value (FF): " << minCutFF;
      LOG(INFO) << "Sparsifier Cut value (PR): " << sparsifierMinCut;
      //   LOG(INFO) << "Sparsifier Cut value (FF): " << sparsifierMinCut_;
      LOG(INFO) << "Quality: " << (double)sparsifierMinCut / (double)minCut;
      // assert(minCutPR == minCutFF && "min cuts must match.");
      assert(sparsifierMinCut == sparsifierMinCut_ && "min cuts must match.");

      // LOG(INFO) << ss;
    }
  } else {
    for (auto size : defaultSizes) {
      if (size >= g->size()) {
        continue;
      }

      for (int i = 0; i < FLAGS_numRepetitions; i++) {
        LOG(INFO) << i;
        auto [ss, ts] = sampler.randomSubset(size);
        LOG(INFO) << ss.size() << " " << ts.size();

        auto minCut = minCutPartition(g, ss);
        // auto minCutPR = prg.maxFlow(ss, ts);
        // auto minCutFF = prgf.maxFlow(ss, ts);
        VLOG(0) << "min cut done.";

        // auto sparsifierMinCut_ = prf.maxFlow(ss, ts);

        auto [s, t] = sparsifier->addSourceSink(ss, ts);
        auto prs = PushRelabel(sparsifier);
        auto sparsifierMinCut = prs.maxFlow(s, t);
        sparsifier->removeSourceSink();

        VLOG(0) << "sparsifier done.";

        LOG(INFO) << "Cut value: " << minCut;
        // LOG(INFO) << "Cut value (PR): " << minCutPR;
        // LOG(INFO) << "Cut value (FF): " << minCutFF;
        LOG(INFO) << "Sparsifier Cut value (PR): " << sparsifierMinCut;
        // LOG(INFO) << "Sparsifier Cut value (FF): " << sparsifierMinCut_;
        LOG(INFO) << "Quality: " << (double)sparsifierMinCut / (double)minCut;
        // assert(minCutPR == minCutFF && "min cuts must match.");
        assert(sparsifierMinCut == sparsifierMinCut_ && "min cuts must match.");
      }
    }
  }
}

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);

  gflags::SetUsageMessage("Expander Hierarchy");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (argc != 2 && !FLAGS_stdin) {
    VLOG(0) << "Please provide a filename holding a graph file.";
    return 1;
  }

  auto randomGen = configureRandomness(FLAGS_seed);

  DynamicGraph::Graph *g;

  if (FLAGS_stdin) {
    VLOG(0) << "Reading input.";
    g = readDynamicGraphMTX();
    VLOG(0) << "Finished reading input.";
  } else if (FLAGS_chaco) {
    // TODO
    // g = readGraphChaco(argv[1]);
  } else {
    g = readDynamicGraphMTX(argv[1]);
  }

  LOG(INFO) << "Graph has " << g->size() << " vertices and " << g->getNumEdges()
            << " edges.";

  // compute the expander hierarchy and get the sparsifier
  const auto start = std::chrono::high_resolution_clock::now();
  auto samplingType = static_cast<SamplingAlgo>(FLAGS_subsetType);

  auto solver = RandomWalk::ExpanderHierarchy<double, DynamicGraph::Graph,
                                              DynamicGraph::WeightedGraph>(
      g, randomGen.get(), FLAGS_potential, FLAGS_boost);
  solver.computeExpanderHierarchy(FLAGS_phi);

  LOG(INFO) << "Computed the expander hierarchy. Getting sparsifier.";

  auto sparsifier = solver.getSparsifier();

  const auto end = std::chrono::high_resolution_clock::now();

  LOG(INFO) << "Got sparsifier. Converting to graph.";

  auto sg = sparsifier.toGraph();

  LOG(INFO) << "Converted to graph.";

  // now, we sample subsets and calculate the cuts for these subsets

  if (FLAGS_writeToFile) {
    std::filesystem::path f("quality-experiments.csv");
    if (!std::filesystem::exists(f)) {
      std::ofstream file(f, std::ios::out | std::ios::app);

      LOG(INFO) << "Experiment file does not exist, writing header.";

      file << "Graph Name,Phi,Potential,Running "
              "Time,Size,Value,Actual,Ratio"
           << std::endl;
    }

    std::ofstream file(f, std::ios::out | std::ios::app);

    if (!file) {
      LOG(FATAL) << "There was a problem opening the file.";
    } else {
      unsigned long time =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();

      if (samplingType == SamplingAlgo::Random) {
        LOG(INFO) << "Sampling " << FLAGS_numRepetitions << " random subsets.";

        RandomSubset sampler(g->size());
        qualityComparison(g, sg, sampler, file, time);
      }

      if (samplingType == SamplingAlgo::Greedy) {
        LOG(INFO) << "Sampling " << FLAGS_numRepetitions << " greedy subsets.";

        GreedySubset sampler(g);
        qualityComparison(g, sg, sampler, file, time);
      }
    }
  } else {
    if (samplingType == SamplingAlgo::Random) {
      LOG(INFO) << "Sampling " << FLAGS_numRepetitions << " random subsets.";

      RandomSubset sampler(g->size());
      qualityComparison(g, sg, sampler);
    }

    if (samplingType == SamplingAlgo::Greedy) {
      LOG(INFO) << "Sampling " << FLAGS_numRepetitions << " greedy subsets.";

      GreedySubset sampler(g);
      qualityComparison(g, sg, sampler);
    }
  }

  // sparsifier.printDebug();

  return 0;
}
