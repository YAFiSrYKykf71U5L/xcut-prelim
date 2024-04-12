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
#include "lib/random_walk.hpp"
#include "lib/util/cut_evaluators.hpp"
#include "util.hpp"

DEFINE_uint32(seed, 0,
              "Seed randomness with any positive integer. Default value '0' "
              "means a random seed will be chosen based on system time.");
DEFINE_uint32(experimentType, 0, "Type of the experiment we are running.");
DEFINE_uint32(numParts, 2, "Number of parts in multicut objective.");
DEFINE_double(phi, 0.0,
              "Value of \\phi such that expansion of each cluster is at least "
              "\\phi. If set to zero, the expander decomposition tries to find "
              "a good value itself.");
DEFINE_uint32(boost, 1,
              "Value of the boundary-linkedness of each cluster. In other "
              "words, this is how much more we weigh outgoing edges versus "
              "internal ones.");
DEFINE_double(potential, 0.0001,
              "Potential of the random walk at which we accept convergence and "
              "certify a component is an expander.");
DEFINE_string(name, "unknown", "Name of the graph that is being processed.");
DEFINE_bool(hierarchy, false, "Compute an expander hierarchy.");
DEFINE_bool(stdin, false, "Read the graph from standard input.");
DEFINE_bool(chaco, false,
            "Input graph is given in the Chaco graph file format");
DEFINE_bool(writeToFile, false,
            "Write the output to a csv file, append if it already exists.");

constexpr std::string_view RWVERSION = "0.0.1";

enum ExperimentType {
  SparseCut,
  LowConductanceCut,
  NormalizedCut,
  GreedyNormalizedCut,
  DynamicNormalizedCut,
};

std::vector<int> numClusterSizes{2, 4, 8, 16, 32, 64, 128};
// std::vector<int> numClusterSizes{30};

/**
 * Function to validate cut produced.
 */
double validator(std::shared_ptr<ArrayGraph::Graph> graph,
                 std::vector<int> const &cut, int numComponents,
                 ExperimentType type) {
  std::vector<int> cutEdges(numComponents, 0), volumes(numComponents, 0),
      sizes(numComponents, 0);

  // iterate all edges
  for (int u = 0; u < graph->size(); u++) {
    sizes[cut[u]]++;

    volumes[cut[u]] += graph->getInitialDegree(u);

    for (auto e = graph->cbeginEdges(u); e != graph->cendEdges(u); e++) {
      auto v = *e;

      if (cut[u] != cut[v]) {
        cutEdges[cut[u]]++;
      }
    }
  }

  double target;

  switch (type) {
    case SparseCut:
      assert(numComponents == 2);
      target = (double)cutEdges[0] / (double)std::min(sizes[0], sizes[1]);
      break;

    case LowConductanceCut:
      assert(numComponents == 2);
      target = (double)cutEdges[0] / (double)std::min(volumes[0], volumes[1]);
      break;

    case NormalizedCut:
      assert(numComponents == 2);
    case GreedyNormalizedCut:
    case DynamicNormalizedCut:
      target = 0.0;

      for (int i = 0; i < numComponents; i++) {
        VLOG(1) << "Edges: " << cutEdges[i] << " Volume: " << volumes[i];
        target += (double)cutEdges[i] / (double)volumes[i];
      }
      break;

    default:
      VLOG(0) << "Unknown Experiment Type. Please use a valid type.";
  }

  return target;
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

  ArrayGraph::Graph *g;

  if (FLAGS_stdin) {
    VLOG(0) << "Reading input.";
    g = readGraphMTX();
    VLOG(0) << "Finished reading input.";
  } else if (FLAGS_chaco) {
    g = readGraphChaco(argv[1]);
  } else {
    g = readGraphMTX(argv[1]);
  }

  VLOG(0) << "Graph has " << g->size() << " vertices and " << g->getNumEdges()
          << " edges.";

  auto dfsResult = g->dfs();

  VLOG(0) << "Graph has " << dfsResult.numComponents
          << " connected components.";

  const auto start = std::chrono::high_resolution_clock::now();
  auto experimentType = static_cast<ExperimentType>(FLAGS_experimentType);

  if (FLAGS_hierarchy) {
    auto solver = RandomWalk::ExpanderHierarchy<double, ArrayGraph::Graph,
                                                ArrayGraph::WeightedGraph>(
        g, randomGen.get(), FLAGS_potential, FLAGS_boost);
    solver.computeExpanderHierarchy(FLAGS_phi);
    const auto endExpander = std::chrono::high_resolution_clock::now();

    auto sparsifier = solver.getSparsifier();

    // writeSparsifier(sparsifier, "sparsifier.txt");

    CutSparsifier::CutResult result;

    for (auto clusters : numClusterSizes) {
      const auto startRefine = std::chrono::high_resolution_clock::now();

      switch (experimentType) {
        case ExperimentType::SparseCut:
          result = sparsifier.sparseCut();
          VLOG(0) << "Sparsifier's sparsest cut is of sparsity " << result.cost;
          break;
        case ExperimentType::LowConductanceCut:
          result = sparsifier.lowConductanceCut();
          VLOG(0) << "Sparsifier's lowest conductance cut is of conductance "
                  << result.cost;
          break;
        case ExperimentType::NormalizedCut:
          result = sparsifier.normalizedCut();
          VLOG(0) << "Sparsifier's best normalized cut is of value "
                  << result.cost;
          break;
        case ExperimentType::GreedyNormalizedCut:
          // result = sparsifier.greedyNormalizedCut(FLAGS_numParts);
          result = sparsifier.greedyNormalizedCut(clusters);
          VLOG(0) << "Sparsifier's best normalized cut for k = " << clusters
                  << " is of value " << result.cost << " (greedy)";
          break;
        case ExperimentType::DynamicNormalizedCut:
          result = sparsifier.normalizedCutDynamic(clusters);
          VLOG(0) << "Sparsifier's best normalized cut for k = "
                  << FLAGS_numParts << " is of value " << result.cost
                  << " (dynamic)";
          break;
        default:
          VLOG(0) << "Experiment type not known.";
          exit(1);
      }

      auto graphHierarchy = solver.getHierarchy();
      graphHierarchy.setLeaders(result.edges);
      graphHierarchy.setTaboo(result.tabooNodes);

      while (graphHierarchy.size() > 1) {
        auto gp = graphHierarchy.popAndProject();
        // delete gp;
        graphHierarchy.refine();
      }
      const auto endRefine = std::chrono::high_resolution_clock::now();

      assert(g == graphHierarchy.getBaseGraph());

      double actual;

      switch (experimentType) {
        case ExperimentType::SparseCut:
          actual = sparseCut(g, g->getClusterArray());
          break;
        case ExperimentType::LowConductanceCut:
          actual = lowConductanceCut(g, g->getClusterArray());
          break;
        case NormalizedCut:
        case GreedyNormalizedCut:
        case DynamicNormalizedCut:
          actual = normalizedCut(g, g->getClusterArray(), clusters);
          break;
        default:
          LOG(ERROR)
              << "Experiment type not recognized. This should not happpen.";
          assert(false);
      }

      VLOG(0) << "Validator: " << actual;

      if (FLAGS_writeToFile) {
        std::filesystem::path f("hierarchy-experiments.csv");
        if (!std::filesystem::exists(f)) {
          std::ofstream file("hierarchy-experiments.csv",
                             std::ios::out | std::ios::app);

          VLOG(0) << "Experiment file does not exist, writing header.";

          file
              << "Graph Name,Experiment Type,Phi,Potential,Parts,Decomposition "
                 "Time,Solving Time,Value,Actual"
              << std::endl;
        }

        std::ofstream file("hierarchy-experiments.csv",
                           std::ios::out | std::ios::app);

        if (!file) {
          VLOG(0) << "There was a problem opening the file.";
        } else {
          VLOG(0) << "Appending to experiment file.";
          file << FLAGS_name << ",";
          file << FLAGS_experimentType << ",";
          file << FLAGS_phi << ",";
          file << FLAGS_potential << ",";
          file << clusters << ",";
          file << std::chrono::duration_cast<std::chrono::milliseconds>(
                      endExpander - start)
                      .count()
               << ",";
          file << std::chrono::duration_cast<std::chrono::milliseconds>(
                      endRefine - startRefine)
                      .count()
               << ",";
          file << result.cost << ",";
          file << actual << std::endl;
        }
      }
    }

  } else {
    auto solver = RandomWalk::ExpanderDecomposition<ArrayGraph::Graph, double>(
        g, randomGen.get(), FLAGS_potential, FLAGS_boost);
    solver.computeExpanderDecomp(FLAGS_phi);

    const auto end = std::chrono::high_resolution_clock::now();

    if (FLAGS_writeToFile) {
      std::ofstream file("decomposition-experiments.csv",
                         std::ios::out | std::ios::app);

      std::filesystem::path f("hierarchy-experiments.csv");
      if (!std::filesystem::exists(f)) {
        VLOG(0) << "Experiment file does not exist, writing header.";

        file << "Graph Name,Experiment Type,Phi,Potential,Running Time,Value"
             << std::endl;
      }

      if (!file) {
        VLOG(0) << "There was a problem opening the file.";
      } else {
        file << FLAGS_name << ",";
        file << FLAGS_experimentType << ",";
        file << FLAGS_phi << ",";
        file << FLAGS_potential << ",";
        file << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                     start)
                    .count()
             << std::endl;
        // file << value << std::endl;
      }
    }
  }

  delete g;
  return 0;
}
