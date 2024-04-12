#ifndef MAX_FLOW_H_
#define MAX_FLOW_H_

#include <algorithm>
#include <limits>
#include <memory>
#include <unordered_set>

#include "bfs.hpp"

template <typename G>
class FordFulkerson {
 private:
  std::shared_ptr<G> graph;
  std::vector<std::vector<int>> flow;
  BFS<G> bfs;
  bool fresh;

  void reset() {
    for (auto &f : flow) {
      std::fill(f.begin(), f.end(), 0);
    }
  }

  /**
   * Get the maximum flow along the BFS path between s and t.
   */
  int pathFlow(int s, int t) {
    int routable = std::numeric_limits<int>::max();

    while (s != t) {
      auto i = bfs.getParent(t);

      // VLOG(0) << "Edge (" << t << ", " << graph->getEdge(t, i) << ")";
      // VLOG(0) << "flow: " << flow[t][i] << " out of "
      //         << graph->getEdgeWeight(t, i);

      // need to add since we're going the reverse direction here so the sign is
      // flipped
      routable = std::min(routable, graph->getEdgeWeight(t, i) + flow[t][i]);

      t = graph->getEdge(t, i);
    }

    return routable;
  }

 public:
  FordFulkerson(std::shared_ptr<G> graph)
      : graph(graph), flow(graph->size()), bfs(graph), fresh(true) {
    for (int u = 0; u < graph->size(); u++) {
      flow[u].resize(graph->getNumNeighbors(u), 0);
    }
  }

  /**
   * Compute the maximum flow between vertex s and vertex t.
   */
  int maxFlow(int s, int t) {
    if (!fresh) {
      reset();
    }
    fresh = false;

    int maxFlow = 0;

    // each iteration begins with a fresh BFS query
    while (bfs.connectedResidual(s, t, flow)) {
      // get the amount of routable flow

      VLOG(4) << "Found path in residual graph.";
      int routable = pathFlow(s, t);

      maxFlow += routable;

      VLOG(4) << "Can route " << routable << " flow along this path.";

      // route the flow
      int u = t;

      while (u != s) {
        // i is an edge on the path from t to s
        auto i = bfs.getParent(u);

        auto v = graph->getEdge(u, i);
        auto rev = graph->getReverse(u, i);

        flow[u][i] -= routable;
        flow[v][rev] += routable;

        u = v;
      }
    }

    return maxFlow;
  }

  /**
   * Find the max flow between two terminal sets of vertices.
   */
  int maxFlow(std::vector<int> &ss, std::vector<int> &ts) {
    // check if the terminal sets are disjoint
    std::sort(ss.begin(), ss.end());
    std::sort(ts.begin(), ts.end());

    // std::vector<int> intersection;

    // std::set_intersection(ss.begin(), ss.end(), ts.begin(), ts.end(),
    //                       std::back_inserter(intersection));

    // if (!intersection.empty()) {
    //   LOG(ERROR) << "Sets of terminal vertices must be disjoint.";
    //   return -1;
    // }

    if (!fresh) {
      reset();
    }
    fresh = false;

    int maxFlow = 0;

    std::unordered_set<int> targets;
    for (auto t : ts) {
      targets.insert(t);
    }

    while (bfs.connectedResidual(ss, targets, flow)) {
      VLOG(4) << "Found path in residual graph.";

      auto [s, t] = bfs.getPair();

      auto routable = pathFlow(s, t);

      maxFlow += routable;

      VLOG(4) << "Can route " << routable << " flow between " << s << " and "
              << t << ".";

      int u = t;

      while (u != s) {
        auto i = bfs.getParent(u);

        auto v = graph->getEdge(u, i);
        auto rev = graph->getReverse(u, i);

        flow[u][i] -= routable;
        flow[v][rev] += routable;

        u = v;
      }
    }
    return maxFlow;
  }
};

#endif  // MAX_FLOW_H_
