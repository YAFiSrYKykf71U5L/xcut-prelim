// These first two are to get it to build with bazel on Arch, may be redundant
// on other systems.
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <glog/stl_logging.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "lib/io/graph_io.hpp"
#include "lib/util/cut_evaluators.hpp"
#include "util.hpp"

DEFINE_string(name, "unknown", "Name of the graph that is being processed.");
DEFINE_string(algo, "unknown",
              "Name of the algorithm that produced the partition.");
DEFINE_double(time, 0.0, "Execution time of the partitioner.");
// DEFINE_uint32(numParts, 0, "Number of parts in the partition.");
DEFINE_bool(writeToFile, false,
            "Write the output to a csv file, append if it already exists.");

constexpr std::string_view EVALVERSION = "0.0.1";

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);

  gflags::SetUsageMessage(
      "Utility for computing the normalized cut value of a partition of a "
      "graph.");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (argc != 3) {
    LOG(FATAL)
        << "Please provide a filename holding a graph file and a partition.";
    return 1;
  }

  // std::shared_ptr<ArrayGraph::Graph> g;

  // g = readGraphMTX(argv[1]);
  auto g = readGraphChaco(argv[1]);
  LOG(INFO) << "Graph has " << g->size() << " vertices and " << g->getNumEdges()
            << " edges.";

  LOG(INFO) << "Reading partition from file " << argv[2];
  auto partition = readPartition(argv[2]);
  auto numParts = *std::max_element(partition.begin(), partition.end()) + 1;
  LOG(INFO) << "Partition has " << numParts << " parts.";

  LOG(INFO) << "Computing normalized cut.";
  auto value = normalizedCut(g, partition, numParts);
  LOG(INFO) << "Normalized cut is of value " << value << ". Writing to file.";

  std::filesystem::path f("comparison-experiments.csv");
  if (!std::filesystem::exists(f)) {
    std::ofstream file("comparison-experiments.csv",
                       std::ios::out | std::ios::app);

    VLOG(0) << "Experiment file does not exist, writing header.";

    file << "Graph Name,Experiment Type,Parts,Time,Value" << std::endl;
  }

  std::ofstream file("comparison-experiments.csv",
                     std::ios::out | std::ios::app);

  if (!file) {
    VLOG(0) << "There was a problem opening the file.";
    delete g;
    return 1;
  } else {
    VLOG(0) << "Appending to experiment file.";
    file << FLAGS_name << ",";
    file << FLAGS_algo << ",";
    file << numParts << ",";
    file << FLAGS_time << ",";
    file << value << std::endl;
    delete g;
    return 0;
  }
}
