#include <limits>
#include <memory>
#include <queue>

#include "bfs.hpp"

/**
 * Class containing an implementation fo the push-relabel max-flow algorithm,
 * using the same algorithmic logic as the algorithm found in KaHiP
 */
template <typename G>
class PushRelabel {
 private:
  static constexpr double UPDATE_FREQ = 0.51;
  static constexpr int RELABEL_WORK = 9;
  static constexpr int NTE_WORK = 4;

  std::shared_ptr<G> graph;

  std::vector<std::vector<int>> flow;
  std::vector<int> excess;
  std::vector<int> labels;

  std::vector<int> sizes;
  std::vector<bool> active;

  std::queue<int> queue;

  BFS<G> bfs;
  bool fresh;
  int work;

  void reset() {
    for (auto &f : flow) {
      std::fill(f.begin(), f.end(), 0);
    }
    std::fill(excess.begin(), excess.end(), 0);
    std::fill(labels.begin(), labels.end(), 0);
    std::fill(sizes.begin(), sizes.end(), 0);
    std::fill(active.begin(), active.end(), false);

    work = 0;
  }

  void enqueue(int u) {
    if (active[u]) {
      return;
    }
    if (excess[u] > 0) {
      active[u] = true;
      queue.push(u);
    }
  }

  void push(int u, int i) {
    // VLOG(0) << "u: " << u << " i: " << i;
    auto v = graph->getEdge(u, i);
    // VLOG(0) << "v: " << v;

    assert(excess[u] > 0 && "Node u must have excess flow.");

    auto rev = graph->getReverse(u, i);
    auto routable =
        std::min(excess[u], graph->getEdgeWeight(u, i) - flow[u][i]);

    // exit early if nothing can be pushed or the push is invalid as it would
    // push uphill
    if (routable == 0 || labels[u] <= labels[v]) {
      // VLOG(0) << "Can't push from " << u << " to neighbor " << v;
      return;
    }

    // VLOG(0) << "Pushing " << routable << " from " << u << " to " << v;

    flow[u][i] += routable;
    flow[v][rev] -= routable;

    excess[u] -= routable;
    excess[v] += routable;

    enqueue(v);
  }

  void relabel(int u) {
    assert(0 <= u && u < graph->size() &&
           "Node u must be a member of the graph.");
    assert(excess[u] > 0 && "Node u must have excess before relabelling.");

    work += RELABEL_WORK;

    int minLabel = 2 * graph->size();

    for (auto i = 0; i < graph->getNumNeighbors(u); i++) {
      if (graph->getEdgeWeight(u, i) - flow[u][i] > 0) {
        auto v = graph->getEdge(u, i);
        minLabel = std::min(minLabel, labels[v] + 1);
      }
      work++;
    }

    // VLOG(0) << "Relabeling " << u << " from " << labels[u];
    // VLOG(0) << "minLabel: " << minLabel;
    // VLOG(0) << "2 * n: " << 2 * graph->size();

    sizes[labels[u]]--;
    labels[u] = minLabel;
    sizes[labels[u]]++;

    enqueue(u);
  }

  void globalRelabel(int s, int t) {
    bfs.componentResidual(t, flow);

    for (int u = 0; u < graph->size(); u++) {
      if (u == s) {
        continue;
      }

      auto newDist = bfs.getBFSLabel(u);
      if (newDist >= 0) {
        sizes[labels[u]]--;
        labels[u] = newDist;
        sizes[labels[u]]++;
      }
    }
  }

  void gapHeuristic(int target) {
    for (int u = 0; u < graph->size(); u++) {
      if (labels[u] < target) {
        continue;
      }

      sizes[labels[u]]--;
      labels[u] = std::max(labels[u], graph->size());
      sizes[labels[u]]++;

      enqueue(u);
    }
  }

  void discharge(int u) {
    // try to push all excess out of the node
    for (int i = 0; i < graph->getNumNeighbors(u) && excess[u] > 0; i++) {
      push(u, i);
    }

    if (excess[u] > 0) {
      // gap heuristic
      if (sizes[labels[u]] == 1 && labels[u] < graph->size()) {
        gapHeuristic(labels[u]);
      } else {
        relabel(u);
      }
    }
  }

 public:
  PushRelabel(std::shared_ptr<G> graph)
      : graph(graph),
        flow(graph->size()),
        excess(graph->size(), 0),
        labels(graph->size(), 0),
        sizes(2 * graph->size(), 0),
        active(graph->size(), false),
        bfs(graph),
        fresh(true),
        work(0) {
    for (int u = 0; u < graph->size(); u++) {
      flow[u].resize(graph->getNumNeighbors(u), 0);
    }

    sizes[0] = graph->size() - 1;
    sizes[graph->size()] = 1;
  }

  int maxFlow(int s, int t) {
    assert(s < graph->size() && t < graph->size());

    if (!fresh) {
      reset();

      // after arrays have been cleaned need to set up sources
      sizes[0] = graph->size() - 1;
      sizes[graph->size()] = 1;
    }
    fresh = false;

    labels[s] = graph->size();
    active[s] = true;
    active[t] = true;

    for (int i = 0; i < graph->getNumNeighbors(s); i++) {
      excess[s] += graph->getEdgeWeight(s, i);

      push(s, i);
    }

    // VLOG(0) << "Finished initial push.";

    globalRelabel(s, t);

    // VLOG(0) << "Finished relabel.";

    const int workLeft = NTE_WORK * graph->size() + graph->getNumEdges();

    assert(excess[s] == 0 && "Excess at source must be zero.");

    while (!queue.empty()) {
      auto u = queue.front();
      //   VLOG(0) << "Got " << u << ". Discharging.";
      //   VLOG(0) << "Excess at " << u << ": " << excess[u];
      //   VLOG(0) << "Label at " << u << ": " << labels[u];
      queue.pop();

      active[u] = false;

      discharge(u);

      //   if (work > UPDATE_FREQ * workLeft) {
      //     VLOG(0) << "DOing global relable.";
      //     globalRelabel(s, t);
      //     work = 0;
      //   }
    }

    return excess[t];
  }
};
