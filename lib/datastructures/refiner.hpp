#include "graph_hierarchy.hpp"

struct Clustering {
  std::vector<int> c;
  int numClusters;
  std::vector<int> clusterSize;
  std::vector<int> clusterVolume;
  std::vector<int> clusterCut;

  template <typename G>
  void addCluster(G *graph, int u) {
    auto prevIndex = clustering[u];

    auto vertexWeight = graph->getVertexWeight(u);
    auto vertexDegree = graph->getInitialDegree(u);

    clustering[u] = numClusters++;

    /**
     * TODO: cluster updating
     */
  }
};

template <typename G, typename GW>
class Refiner {
 private:
  // the graph hierarchy
  GraphHierarchy<G, GW> *hierarchy;

 public:
  template <typename G>
  void refineClustering(G *graph, Clustering &clustering) {}

  /**
   * Project the coarser clustering onto a finer graph.
   */
  template <typename G, typename H>
  Clustering projectClusteringFromParent(G *graph, H *coarser, int numClusters,
                                         Clustering const &coarserClustering) {
    Clustering clustering;

    clustering.c.resize(graph->size());

    // project the clusters
    for (int u = 0; u < graph->size(); u++) {
      clustering.c[u] = coarserClustering.c[graph->getCutIndex[u]];
    }

    return clustering;
  }

  template <typename G>
  void setLeaders(G *graph, Clustering &clustering, int level,
                  std::vector<std::pair<int, int>> &leaders) {
    while (!leaders.empty() && leaders.back().first == level) {
      auto [l, u] = leaders.back();

      VLOG(2) << "Setting leader " << l << ", " << u;

      clustering.addCluster(graph, u);
    }
  }

  template <typename G, typename H>
  Clustering setupRefinement(G *graph, H *coarser,
                             int num std::vector < std::pair<int, int> &
                                 leaders){
      auto clustering = projectClusteringFromParent(graph, coarser, )}

  std::vector<int> refine(std::vector<std::pair<int, int>> leaders,
                          std::vector<std::pair<int, int>> taboos) {}
};