#ifndef BFS_H_
#define BFS_H_

#include <algorithm>
#include <memory>
#include <queue>
#include <unordered_set>

/**
 * File containing a BFS algorithm to be used with graphs.
 *
 * TODO: properly parametrize the template, add BFS tree array?
 */

template <typename G>
class BFS {
 private:
  std::shared_ptr<G> graph;
  std::vector<bool> visited;
  std::vector<int> bfsLabels;

  std::vector<int> parents;

  // contains the index of the edge in the graph
  std::vector<int> tree;

  std::pair<int, int> found;

  bool ready;

  /**
   * Reset the data structures so another query can be made.
   */
  void reset() {
    std::fill(visited.begin(), visited.end(), false);
    std::fill(bfsLabels.begin(), bfsLabels.end(), -1);
    std::fill(tree.begin(), tree.end(), -1);
    ready = true;
  }

 public:
  BFS(std::shared_ptr<G> graph)
      : graph(graph),
        visited(graph->size(), false),
        bfsLabels(graph->size(), -1),
        parents(graph->size(), -1),
        tree(graph->size(), -1),
        ready(true) {}

  /**
   * Make a connectivity query for s and t. The algorithm does not terminate
   * early if t is found unless the terminateEarly condition is set to true.
   */
  bool connected(int s, int t) {
    assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");
    std::queue<int> queue;

    if (!ready) {
      reset();
    }

    ready = false;

    queue.push(s);
    visited[s] = true;
    bfsLabels[s] = 0;

    while (!queue.empty()) {
      auto u = queue.front();
      queue.pop();

      if (u == t) {
        return true;
      }

      for (auto i = 0; i < graph->getNumNeighbors(u); i++) {
        auto v = graph->getEdge(u, i);

        if (!visited[v]) {
          queue.push(v);
          visited[v] = true;
          bfsLabels[v] = bfsLabels[u] + 1;
          tree[v] = graph->getReverse(i);
        }
      }
    }

    return false;
  }

  /**
   * Run BFS on the component containing s.
   */
  void component(int s) {
    assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");
    std::queue<int> queue;

    if (!ready) {
      reset();
    }

    ready = false;

    queue.push(s);
    visited[s] = true;
    bfsLabels[s] = 0;

    while (!queue.empty()) {
      auto u = queue.front();
      queue.pop();

      for (auto i = 0; i < graph->getNumNeighbors(u); i++) {
        auto v = graph->getEdge(u, i);

        if (!visited[v]) {
          queue.push(v);
          visited[v] = true;
          bfsLabels[v] = bfsLabels[u] + 1;
          tree[v] = graph->getReverse(u, i);
        }
      }
    }
  }

  /**
   * Check if s and t are connected in the residual graph of a given flow. The
   * flow need not be valid.
   */
  bool connectedResidual(int s, int t,
                         std::vector<std::vector<int>> const &flow) {
    assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");

    std::queue<int> queue;

    if (!ready) {
      reset();
    }

    ready = false;

    queue.push(s);
    visited[s] = true;
    bfsLabels[s] = 0;

    while (!queue.empty()) {
      auto u = queue.front();

      // VLOG(0) << "Visited " << u;
      queue.pop();

      if (u == t) {
        return true;
      }

      for (int i = 0; i < graph->getNumNeighbors(u); i++) {
        auto v = graph->getEdge(u, i);
        auto capacity = graph->getEdgeWeight(u, i) - flow[u][i];

        if (!visited[v] && capacity > 0) {
          queue.push(v);
          visited[v] = true;
          bfsLabels[v] = bfsLabels[u] + 1;
          tree[v] = graph->getReverse(u, i);
        }
      }
    }

    return false;
  }

  /**
   * Check if s and t are connected in the residual graph of a given flow. The
   * flow need not be valid.
   */
  void componentResidual(int s, std::vector<std::vector<int>> const &flow) {
    assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");

    std::queue<int> queue;

    if (!ready) {
      reset();
    }

    ready = false;

    queue.push(s);
    visited[s] = true;
    bfsLabels[s] = 0;

    while (!queue.empty()) {
      auto u = queue.front();
      queue.pop();

      //   VLOG(0) << "Got " << u;

      for (int i = 0; i < graph->getNumNeighbors(u); i++) {
        // VLOG(0) << "i: " << i;
        // VLOG(0) << "numNeighbors: " << graph->getNumNeighbors(u);

        auto v = graph->getEdge(u, i);
        auto capacity = graph->getEdgeWeight(u, i) - flow[u][i];

        if (!visited[v] && capacity > 0) {
          queue.push(v);
          visited[v] = true;
          bfsLabels[v] = bfsLabels[u] + 1;
          tree[v] = graph->getReverse(u, i);
        }
      }
    }
  }

  //   void componentResidual(std::vector<int> const &ss,
  //                          std::vector<std::vector<int>> const &flow) {
  //     std::queue<int> queue;

  //     if (!ready) {
  //       reset();
  //     }

  //     ready = false;

  //     for (auto s : ss) {
  //       assert(0 <= s && s < graph->size() && "Graph must contain the
  //       vertex.");

  //       queue.push(s);
  //       visited[s] = true;
  //       bfsLabels[s] = 0;
  //     }

  //     while (!queue.empty()) {
  //       auto u = queue.front();
  //       queue.pop();

  //       for (int i = 0; i < graph->getNumNeighbors(u); i++) {
  //         auto v = graph->getEdge(u, i);
  //         auto capacity = graph->getEdgeWeight(u, i) - flow[u][i];

  //         if (!visited[v] && capacity > 0) {
  //           queue.push(v);
  //           visited[v] = true;
  //           bfsLabels[v] = bfsLabels[u] + 1;
  //           tree[v] = graph->getReverse(u, i);
  //         }
  //       }
  //     }
  //   }

  /**
   * Find a path in the residual graph between two sets of terminal vertices.
   */
  bool connectedResidual(std::vector<int> const &ss,
                         std::unordered_set<int> const &ts,
                         std::vector<std::vector<int>> const &flow) {
    std::queue<int> queue;

    if (!ready) {
      reset();
      std::fill(parents.begin(), parents.end(), -1);
    }

    ready = false;

    for (auto s : ss) {
      assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");

      queue.push(s);
      visited[s] = true;
      bfsLabels[s] = 0;
      parents[s] = s;
    }

    while (!queue.empty()) {
      auto u = queue.front();
      queue.pop();

      if (ts.count(u)) {
        found = {parents[u], u};
        return true;
      }

      for (int i = 0; i < graph->getNumNeighbors(u); i++) {
        auto v = graph->getEdge(u, i);
        auto capacity = graph->getEdgeWeight(u, i) - flow[u][i];

        if (!visited[v] && capacity > 0) {
          queue.push(v);
          visited[v] = true;
          bfsLabels[v] = bfsLabels[u] + 1;
          tree[v] = graph->getReverse(u, i);
          parents[v] = parents[u];
        }
      }
    }

    return false;
  }

  /**
   * Returns a boolean indicating whether the vertex was visited in the last
   * query.
   */
  bool getVisited(int u) const { return visited[u]; }

  /**
   * Gets the BFS label of vertex u. It is -1 if it was not visited in the last
   * query.
   */
  int getBFSLabel(int u) const { return bfsLabels[u]; }

  /**
   * Get the parent of the vertex u in the BFS tree.
   */
  int getParent(int u) const { return tree[u]; }

  std::pair<int, int> getPair() const { return found; }
};

#endif  // BFS_H_