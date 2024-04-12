#ifndef GRAPH_HIERARCHY_H_
#define GRAPH_HIERARCHY_H_

#include <algorithm>
#include <memory>
#include <stack>
#include <variant>
#include <vector>

#include "lib/datastructures/array_graph.hpp"

template <typename G, typename GW>
class GraphHierarchy {
 private:
  G *baseGraph;
  std::vector<GW *> graphStack;

  /**
   * Number of clusters in the current graph.
   */
  int numClusters;

  /**
   * Label of parts from the hierarchy that are their own clusters in the k-way
   * partition. First element is level of the partition, second label is the
   * part.
   */
  std::vector<std::pair<int, int>> leaderVertices;

  /**
   * Vertex ids that are taboo, meaning that all their children are cut off,
   * which results in an empty partition set, if the vertex ends up in its own
   * partition.
   */
  std::vector<std::pair<int, int>> tabooNodes;

 public:
  GraphHierarchy(G *g) : baseGraph(g), numClusters(1) {}

  /**
   * All graphs higher up must be weighted, so there is no way to push others.
   */
  void push(GW *g) { graphStack.push_back(g); }

  void setLeaders(std::vector<std::pair<int, int>> leaders) {
    leaderVertices = std::move(leaders);

    std::sort(leaderVertices.begin(), leaderVertices.end(),
              [](std::pair<int, int> a, std::pair<int, int> b) {
                return a.first < b.first;
              });
  }

  void setTaboo(std::vector<std::pair<int, int>> taboo) {
    tabooNodes = std::move(taboo);

    std::sort(tabooNodes.begin(), tabooNodes.end(),
              [](std::pair<int, int> a, std::pair<int, int> b) {
                return a.first < b.first;
              });
  }

  GW *peek() const { return graphStack.back(); }

  GW *popAndProject() {
    VLOG(2) << "Popping graph at level " << graphStack.size();

    auto parent = graphStack.back();

    VLOG(7) << "Parent size: " << parent->size();

    graphStack.pop_back();

    if (size() > 1) {
      auto current = graphStack.back();

      current->setClustering(numClusters, parent->getClusterArray());
    } else {
      baseGraph->setClustering(numClusters, parent->getClusterArray());
    }

    // finally paint all new ones until there is no more to paint
    if (leaderVertices.empty()) {
      return parent;
    }

    if (!leaderVertices.empty()) {
      VLOG(4) << "top leader: " << leaderVertices.back();
    }

    while (!leaderVertices.empty() &&
           leaderVertices.back().first == size() - 1) {
      auto [l, u] = leaderVertices.back();

      VLOG(4) << "Setting element " << l << ", " << u;

      if (size() > 1) {
        auto current = graphStack.back();
        current->newCluster(u);
      } else {
        baseGraph->newCluster(u);
      }
      numClusters++;

      leaderVertices.pop_back();
    }

    while (!tabooNodes.empty() && tabooNodes.back().first > size() - 1) {
      tabooNodes.pop_back();
    }

    // now set taboo vertices
    while (!tabooNodes.empty() && tabooNodes.back().first == size() - 1) {
      auto [l, u] = tabooNodes.back();

      VLOG(4) << "Making (" << l << ", " << u << ") taboo.";

      if (size() > 1) {
        auto current = graphStack.back();
        current->makeTaboo(u);
      } else {
        baseGraph->makeTaboo(u);
      }

      tabooNodes.pop_back();
    }

    return parent;
  }

  void refine() {
    if (numClusters == 1) {
      VLOG(1) << "Nothing to refine, current level contains only one cluster.";
      return;
    }

    if (size() > 1) {
      VLOG(1) << "Refining graph with " << numClusters << " clusters and "
              << graphStack.back()->size() << " vertices.";
      auto current = graphStack.back();
      current->refine();
    } else {
      VLOG(1) << "Refining graph with " << numClusters << " clusters and "
              << baseGraph->size() << " vertices.";
      baseGraph->refine();
    }
  }

  G *getBaseGraph() { return baseGraph; }

  int size() { return graphStack.size() + 1; }

  void del() {
    for (auto graphPointer : graphStack) {
      delete graphPointer;
    }
  }
};

#endif  // GRAPH_HIERARCHY_H_
