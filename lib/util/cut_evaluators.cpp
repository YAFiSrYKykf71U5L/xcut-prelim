#include "cut_evaluators.hpp"

void computePartitionParameters(ArrayGraph::Graph *graph,
                                std::vector<int> const &cut, int numParts,
                                std::vector<int> &cutEdges,
                                std::vector<int> &sizes,
                                std::vector<int> &volumes) {
  cutEdges.resize(numParts);
  sizes.resize(numParts);
  volumes.resize(numParts);

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
}

int cutSize(ArrayGraph::Graph *graph, std::vector<int> const &cut) {
  std::vector<int> cutEdges, sizes, volumes;

  computePartitionParameters(graph, cut, 2, cutEdges, sizes, volumes);

  assert(cutEdges[0] == cutEdges[1]);

  return cutEdges[0];
}

double sparseCut(ArrayGraph::Graph *graph, std::vector<int> const &cut) {
  std::vector<int> cutEdges, sizes, volumes;

  computePartitionParameters(graph, cut, 2, cutEdges, sizes, volumes);

  return (double)cutEdges[0] / (double)std::min(sizes[0], sizes[1]);
}

double lowConductanceCut(ArrayGraph::Graph *graph,
                         std::vector<int> const &cut) {
  std::vector<int> cutEdges, sizes, volumes;

  computePartitionParameters(graph, cut, 2, cutEdges, sizes, volumes);

  return (double)cutEdges[0] / (double)std::min(volumes[0], volumes[1]);
}

double normalizedCut(ArrayGraph::Graph *graph, std::vector<int> const &cut,
                     int numParts) {
  std::vector<int> cutEdges, sizes, volumes;

  computePartitionParameters(graph, cut, numParts, cutEdges, sizes, volumes);

  double cutValue = 0.0;

  for (int i = 0; i < numParts; i++) {
    VLOG(4) << "Edges: " << cutEdges[i] << " Volume: " << volumes[i];
    cutValue += (double)cutEdges[i] / (double)volumes[i];
  }

  return cutValue;
}