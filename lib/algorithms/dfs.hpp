#ifndef DFS_H_
#define DFS_H_

#include <algorithm>
#include <memory>
#include <stack>

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
  std::vector<int> dfsLabels;
  std::vector<int> tree;
  bool ready;

  /**
   * Reset the data structures so another query can be made.
   */
  void reset() {
    std::fill(visited.begin(), visited.end(), false);
    std::fill(dfsLabels.begin(), dfsLabels.end(), -1);
    std::fill(tree.begin(), tree.end(), -1);
    ready = true;
  }

 public:
  BFS(std::shared_ptr<G> graph)
      : graph(graph),
        visited(graph->size(), false),
        dfsLabels(graph->size(), -1),
        tree(graph->size(), -1),
        ready(true) {}

  /**
   * Make a connectivity query for s and t. The algorithm does not terminate
   * early if t is found unless the terminateEarly condition is set to true.
   */
  bool connected(int s, int t, bool terminateEarly) {
    assert(0 <= s && s < graph->size() && "Graph must contain the vertex.");

    std::stack<int> stack;
    if (!ready) {
      reset();
    }

    ready = false;

    stack.push(s);
    visited[s] = true;
    bfsLabels[s] = 0;

    while (!stack.empty()) {
      auto u = stack.top();
      stack.pop();

      if (terminateEarly && u == t) {
        return true;
      }

      for (auto it = graph->cbeginEdges(u); it != graph->cendEdges(u); it++) {
        auto v = *it;

        if (!visited[v]) {
          stack.push(v);
          visited[v] = true;
          dfsLabels[v] = dfsLabels[u] + 1;
          tree[v] = u;
        }
      }
    }

    return visited[t];
  }

  /**
   * Returns a boolean indicating whether the vertex was visited in the last
   * query.
   */
  bool visited(int u) { return visited[u]; }

  /**
   * Checks whether the graph is connected.
   */
  bool isConnected() { return std::all_of(visited.begin(), visited.end()); }

  /**
   * Gets the BFS label of vertex u. It is -1 if it was not visited in the last
   * query.
   */
  bool getDFSLabel(int u) { return dfsLabels[u]; }

  /**
   * Get the parent of the vertex u in the DFS tree.
   */
  int getParent(int u) const { return tree[u]; }
};

#endif  // BFS_H_