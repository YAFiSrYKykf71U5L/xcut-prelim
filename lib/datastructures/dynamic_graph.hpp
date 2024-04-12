#ifndef DYNAMIC_GRAPH_H_
#define DYNAMIC_GRAPH_H_

#include <glog/logging.h>

#include <cassert>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #include "lib/datastructures/bucket_queue.hpp"

/**
 * A graph class that is functionally equivalent to the static array graphs but
 * uses vector<vector<V>> for its edge list, this lets us be more judicious with
 * adding edges.
 */
namespace DynamicGraph {

struct DFSResult {
  int numComponents;
  std::vector<int> labels;
  std::vector<int> sizes;

  DFSResult(int n, std::vector<int> labels, std::vector<int> sizes)
      : numComponents(n), labels(std::move(labels)), sizes(std::move(sizes)) {}
};

template <bool w, typename W>
struct WeightedMembers {};

template <typename W>
struct WeightedMembers<true, W> {
  /**
   * The vertex weight corresponding to the volume of the vertices in the
   * original graph.
   */
  std::vector<W> vertexWeights;
  std::vector<std::vector<W>> weights;
};

struct Edge {
  int from, to;

  Edge(int u, int v) : from(u), to(v) {}
};

template <typename wType>
struct WeightedEdge {
  int from, to;
  wType weight;

  WeightedEdge(int u, int v, wType w) : from(u), to(v), weight(w){};
};

/**
 * The type of W is only for weighted graphs.
 */
template <typename V, typename E, bool isWeighted, typename W = int>
class TemplateGraph : WeightedMembers<isWeighted, W> {
 public:
  constexpr static auto Weighted = isWeighted;

 private:
  /**
   * The weight type. Unweighted graphs use integer weights.
   */
  using WeightType = std::conditional_t<isWeighted, W, int>;

  /**
   * The type of a vector holding degree information for different nodes.
   */
  using DegreeVector =
      std::conditional_t<isWeighted, std::vector<W>, std::vector<int>>;

  using EdgeIt = typename std::vector<V>::iterator;
  using cEdgeIt = typename std::vector<V>::const_iterator;

  using wm = WeightedMembers<isWeighted, W>;

  /**
   * Number of vertices in the graph.
   */
  int size_;

  /**
   * Total volume of graph. In the case of unweighted it is the same as
   * numEdges, otherwise it is $\sum_{e \in E} w(e)$.
   */
  WeightType volume;

  /**
   * Number of pieces into which the graph has been cut. Initially set to 1.
   */
  int numCutComponents;

  /**
   * The volume / weight of each component.
   */
  DegreeVector volumes;
  DegreeVector componentWeights;

  /**
   * (Weighted) degrees of the vertices.
   */
  DegreeVector degrees;
  DegreeVector initialDegrees;

  /**
   * Begin of list of neighbors in the array.
   *
   * This must contain n+1 elements, where the last one in the list is equal to
   * m.
   */
  // std::vector<int> start;

  /**
   * Index of cut of expander decomposition that a given vertex belongs to.
   */
  std::vector<int> cut;

  int numClusters;

  /**
   * Cluster that the vertex belongs to. Relevant for the refinement.
   */
  std::vector<int> cluster;

  std::vector<int> clusterSizes;

  /**
   * The vectors holding degree information for the refinement step.
   */
  std::vector<int> externalDegrees;
  std::vector<int> internalDegrees;

  /**
   * Set holding all boundary vertices, i.e. those that have an outgoing edge
   * to another component. This isn't ideal but offers fast insert/delete.
   */
  std::unordered_set<int> boundary;

  /**
   * List of all edges. This is fixed after initialization, so the graph data
   * structure is not dynamic.
   */
  std::vector<std::vector<V>> edges;

  /**
   * Pointer to the reverse edge.
   */
  std::vector<std::vector<int>> reverse;

  bool flowSolving = false;

 public:
  /**
   * Construct an empty graph.
   */
  TemplateGraph() : size_(0), volume(0), numCutComponents(0) {}

  TemplateGraph(int n, const std::vector<E> &es)
      : size_(n),
        volume(2 * es.size()),
        numCutComponents(1),
        volumes(1, 2 * es.size()),
        componentWeights(1, 2 * es.size()),
        degrees(n, 0),
        initialDegrees(),
        // start(n, 0),
        cut(n, 0),
        cluster(n, 0),
        boundary(),
        edges(n),  // this and the rest are vec<vec<X>>
        reverse(n) {
    static_assert(!isWeighted);

    for (auto e : es) {
      degrees[e.from]++;
      degrees[e.to]++;
    }

    initialDegrees = degrees;

    // fill the edge array
    for (auto e : es) {
      reverse[e.to].push_back(edges[e.from].size());
      reverse[e.from].push_back(edges[e.to].size());

      edges[e.to].push_back(e.from);
      edges[e.from].push_back(e.to);

      // VLOG(0) << "Inserted edge: " << e.to << " " << e.from;
    }

    for (int u = 0; u < size(); u++) {
      for (int i = 0; i < edges[u].size(); i++) {
        auto v = edges[u][i];
        assert(u == edges[v][reverse[u][i]] && "Reverse edge must point back.");
      }
    }
  }

  /**
   * Alternative constructor
   */

  // TemplateGraph(int n, int m, std::vector<int> edges, std::vector<int>
  // degrees,
  //               std::vector<int> start)
  //     : size_(n),
  //       volume(2 * m),
  //       totalCutEdges(0),
  //       numCutComponents(1),
  //       volumes(1, 2 * m),
  //       componentWeights(1, 2 * m),
  //       degrees(std::move(degrees)),
  //       initialDegrees(),
  //       cut(n, 0),
  //       cluster(n, 0),
  //       boundary(),
  //       edges(std::move(edges)),
  //       reverse(2 * m),
  // {
  //   static_assert(!isWeighted);
  //   assert(degrees.size() == static_cast<size_t>(n) &&
  //          "Each node must have a degree.");
  //   assert(std::reduce(degrees.begin(), degrees.end(), 0) == 2 * m &&
  //          "Sum of degrees must equal volume.");

  //   // compute the reverse edge indices
  //   for (int u = 0; u < n; u++) {
  //     for (int j = start[u]; j < start[u + 1]; j++) {
  //       auto v = edges[j];

  //       // find the index of the reverse edge
  //       int i = start[v];
  //       while (edges[i] != u) {
  //         ++i;
  //       }

  //       reverse[j] = i;
  //     }
  //   }

  //   for (int u = 0; u < size(); u++) {
  //     for (int i = start[u]; i < start[u + 1]; i++) {
  //       assert(u == edges[reverse[i]] && "Reverse edge must point back.");
  //     }
  //   }
  // }

  /**
   * Weighted graph constructor.
   */
  TemplateGraph(int n, const std::vector<E> &es, const std::vector<W> &ws,
                const std::vector<W> &vws)
      : size_(n),
        volume(2 * std::accumulate(ws.begin(), ws.end(), 0)),
        numCutComponents(1),
        volumes(1, 2 * std::accumulate(ws.begin(), ws.end(), 0)),
        componentWeights(1, std::accumulate(vws.begin(), vws.end(), 0)),
        degrees(n, 0),
        initialDegrees(n, 0),
        cut(n, 0),
        cluster(n, 0),
        boundary(),
        edges(n),
        reverse(n) {
    static_assert(isWeighted);
    assert(ws.size() == es.size() && "Number of edges and weights must match.");

    wm::weights = std::vector<std::vector<W>>(n);
    wm::vertexWeights = vws;

    std::vector<int> numNeighbors(n, 0);

    for (size_t i = 0; i < es.size(); i++) {
      auto e = es[i];
      auto w = ws[i];

      numNeighbors[e.from]++;
      numNeighbors[e.to]++;

      degrees[e.from] += w;
      degrees[e.to] += w;
    }

    initialDegrees.assign(degrees.begin(), degrees.end());

    // fill the edge array
    for (size_t i = 0; i < es.size(); i++) {
      auto e = es[i];
      auto w = ws[i];

      wm::weights[e.to].push_back(w);
      wm::weights[e.from].push_back(w);

      reverse[e.to].push_back(edges[e.from].size());
      reverse[e.from].push_back(edges[e.to].size());

      edges[e.to].push_back(e.from);
      edges[e.from].push_back(e.to);
    }

    for (int u = 0; u < size(); u++) {
      for (size_t i = 0; i < edges[u].size(); i++) {
        auto v = edges[u][i];
        assert(u == edges[v][reverse[u][i]] && "Reverse edge must point back.");
      }
    }
  }

  /**
   * Edge begin and end operator. Iterates over all edges.
   */
  // EdgeIt begin() { return edges.begin(); }
  // EdgeIt end() { return edges.end(); }

  /**
   * Constant edge begin and end operator. Iterates over all edges.
   */
  // cEdgeIt cbegin() const { return edges.cbegin(); }
  // cEdgeIt cend() const { return edges.cend(); }

  /**
   * Iterate the neighbors of vertex u.
   */
  EdgeIt beginEdges(V u) { return edges[u].begin(); }
  EdgeIt endEdges(V u) { return edges[u].end(); }

  /**
   * Const-iterate the neighbors of vertex u.
   */
  cEdgeIt cbeginEdges(V u) const { return edges[u].cbegin(); }
  cEdgeIt cendEdges(V u) const { return edges[u].cend(); }

  /**
   * Iterate the weights of neighbors of vertex u.
   */
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::iterator> beginWeights(V u) {
    return wm::weights[u].begin();
  }

  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::iterator> endWeights(V u) {
    return wm::weights[u].end();
  }

  /**
   * Const-iterate the weights of neighbors of vertex u.
   */
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::const_iterator> cbeginWeights(
      V u) const {
    return wm::weights[u].cbegin();
  }
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::const_iterator> cendWeigths(
      V u) const {
    return wm::weights[u].cend();
  }

  template <bool w = isWeighted>
  std::enable_if_t<w, W> getVertexWeight(V u) {
    return wm::vertexWeights[u];
  }

  template <bool w = isWeighted>
  std::enable_if_t<!w, W> getVertexWeight(V u) {
    return getInitialDegree(u);
  }

  /**
   * Number of vertices in the graph.
   */
  int size() const { return edges.size(); }

  /**
   * Returns the global volume of the graph.
   */
  WeightType getVolume() const { return volume; }

  /**
   * Returns the degree of the vertex, taking into account that edges may have
   * been cut. Does not count cut edges.
   */
  WeightType getDegree(V u) const { return degrees[u]; }

  /**
   * Returns the degree of the vertex in the original graph, without taking the
   * cut structure into account.
   */
  WeightType getInitialDegree(V u) const { return initialDegrees[u]; }

  /**
   * Returns the (weighted) volume of the component containing u
   */
  WeightType getComponentVolume(V u) const { return volumes[cut[u]]; }

  /**
   * Returns the number of cut edges / cut edge weight incident to the
   * vertex. This should not be used for iteration checks.
   */
  WeightType getNumCutEdges(V u) const {
    return initialDegrees[u] - degrees[u];
  }

  /**
   * Returns the current component that u is assigned to.
   */
  int getCutIndex(V u) const { return cut[u]; }

  template <bool t = isWeighted>
  std::enable_if_t<t, WeightType> getEdgeWeight(V u, int i) const {
    return wm::weights[u][i];
  }

  template <bool t = isWeighted>
  std::enable_if_t<!t, WeightType> getEdgeWeight(V u, int i) const {
    // u is a virtual vertex
    if (u >= size_ || getEdge(u, i) >= size_) {
      return 1000;
    } else {
      return 1;
    }
  }

  int getNumNeighbors(V u) const { return edges[u].size(); }

  V getEdge(V u, int i) const { return edges[u][i]; }

  /**
   * Returns a const reference to the current cut array.
   */
  std::vector<int> const &getCutArray() const { return cut; }

  /**
   * Get the index of the reverse edge.
   */
  int getReverse(V u, int j) const { return reverse.at(u).at(j); };

  /**
   * Returns a const reference to the current cluster array.
   */
  std::vector<int> const &getClusterArray() const { return cluster; }

  int getNumClusters() const { return numClusters; }

  void setCluster(V u, int clusterLabel) {
    auto where = cluster[u];

    cluster[u] = clusterLabel;

    clusterSizes[where]--;
    clusterSizes[clusterLabel]++;
  }

  void newCluster(V u) {
    auto where = cluster[u];

    cluster[u] = numClusters;
    numClusters++;

    clusterSizes[where]--;
    clusterSizes.push_back(1);
  }

  void setClustering(int n, std::vector<int> const &parentArray) {
    assert(parentArray.size() == static_cast<size_t>(getNumComponents()) &&
           "Each component exactly one parent.");

    numClusters = n;
    clusterSizes = std::vector<int>(n, 0);

    for (int i = 0; i < size(); i++) {
      assert(parentArray[cut[i]] < numClusters &&
             "Num clusters must be correct.");

      cluster[i] = parentArray[cut[i]];

      clusterSizes[cluster[i]]++;
    }
  }

  /**
   * Returns a const reference to the current cluster array.
   */
  DegreeVector const &getInitialDegreeVector() const { return initialDegrees; }

  /**
   * Used to set the data structure's state to correspond to an externally
   * provided cut vector.
   */
  void setCutArray(int newComponents, std::vector<int> newCut) {
    assert(newCut.size() == static_cast<size_t>(size()) &&
           "Cut vector size must match size of graph.");

    numCutComponents = newComponents;
    cut = std::move(newCut);

    volumes = DegreeVector(numCutComponents, 0);
    componentWeights = DegreeVector(numCutComponents, 0);
    degrees = initialDegrees;

    for (int u = 0; u < size(); u++) {
      volumes[cut[u]] += initialDegrees[u];
      componentWeights[cut[u]] += getVertexWeight(u);

      for (int j = 0; j < edges[u].size(); j++) {
        if (!sameComponent(u, edges[u][j])) {
          degrees[u] -= getEdgeWeight(u, j);
        }
      }
    }
  }

  const std::vector<int> &getInitialDegreeArray() const {
    return initialDegrees;
  }

  // this could be computed on the fly and memorized
  DegreeVector getCutEdges() const {
    DegreeVector cutEdges(numCutComponents, 0);

    for (int u = 0; u < size(); u++) {
      for (int j = 0; j < edges[u].size(); j++) {
        if (!sameComponent(u, edges[u][j]))
          cutEdges[cut[u]] += getEdgeWeight(u, j);
      }
    }
    return cutEdges;
  }

  /**
   * Returns the number of pieces into which the graph has been cut.
   */
  int getNumComponents() const { return numCutComponents; }

  /**
   * Returns the number of edges that are cut by the current partition.
   */
  int getTotalNumCutEdges() const {
    return (volume - std::reduce(degrees.begin(), degrees.end(), 0)) / 2;
  }

  /**
   * Returns the number of total edges in the graph.
   */
  int getNumEdges() const {
    int res = 0;
    for (int i = 0; i < size(); i++) {
      res += edges[i].size();
    }
    return res / 2;
  }

  /**
   * Check whether two vertices are in the same component.
   */
  bool sameComponent(int u, int v) const { return cut[u] == cut[v]; }

  /**
   * Check whether two vertices are in the same cluster.
   */
  bool sameCluster(int u, int v) const { return cluster[u] == cluster[v]; }

  /**
   * Make the elements in the iterable range into a new cut component.
   *
   * Returns the index of the new component.
   *
   * TODO: make this more efficient, one iteration over the iterator should be
   * fine.
   */
  template <typename iterType>
  int updateCutVertices(iterType begin, iterType mid, iterType end) {
    WeightType newVolume = 0;
    WeightType newWeight = 0;
    int partitionIndex = cut[*begin];

    for (iterType it = mid; it != end; ++it) {
      cut[*it] = numCutComponents;
      newVolume += initialDegrees[*it];
      newWeight += getVertexWeight(*it);
    }

    assert(newVolume <= volumes[partitionIndex] &&
           "Volume of a sub-component can not be larger than parent component "
           "volume.");

    for (iterType it = begin; it != end; ++it) {
      int u = *it;

      WeightType newDegree = 0;

      for (size_t j = 0; j < edges[u].size(); j++) {
        if (sameComponent(u, edges[u][j])) {
          newDegree += getEdgeWeight(u, j);
        }
      }
    }

    volumes[partitionIndex] -= newVolume;
    volumes.push_back(newVolume);

    componentWeights[partitionIndex] -= newWeight;
    componentWeights.push_back(newWeight);

    assert(std::reduce(volumes.begin(), volumes.end(), 0) == volume &&
           "No volume can be lost.");

    return numCutComponents++;
  }

  /**
   * Make a new component from a connected component's iterator pair.
   */
  template <typename iterType>
  int makeNewComponentFromConnected(iterType begin, iterType end) {
    assert(dfs(begin, end).numComponents == 1 &&
           "Component must be connected.");

    WeightType newVolume = 0, newWeight = 0;
    int partitionIndex = cut[*begin];

    for (iterType it = begin; it != end; it++) {
      cut[*it] = numCutComponents;
      newVolume += initialDegrees[*it];
      newWeight += getVertexWeight(*it);
    }

    volumes[partitionIndex] -= newVolume;
    volumes.push_back(newVolume);

    componentWeights[partitionIndex] -= newWeight;
    componentWeights.push_back(newWeight);

    return numCutComponents++;
  }

  /**
   * Compute a new graph by contracting the components of the graph.
   *
   * Returns a weighted graph, where weights count the number of edges between
   * components.
   *
   * TODO: use more efficient data representation
   */
  std::shared_ptr<TemplateGraph<V, E, true, WeightType>> contractedGraph() {
    struct PairHash {
      std::size_t operator()(const std::pair<int, int> &p) const {
        return std::hash<size_t>{}(p.first) ^ std::hash<size_t>{}(p.second);
      }
    };

    std::unordered_map<std::pair<int, int>, WeightType, PairHash> edgeWeights;

    std::vector<Edge> newEdges;
    std::vector<WeightType> newWeights;

    if (numCutComponents < 2) {
      VLOG(8) << "Contracted graph to a singleton.";
      // if no cut has been made, return a singleton graph
      return std::make_shared<TemplateGraph<V, E, true, WeightType>>(
          numCutComponents, newEdges, newWeights, componentWeights);
    }

    VLOG(8) << "Contracting the graph, not a singleton.";

    for (int u = 0; u < size(); u++) {
      for (size_t i = 0; i < edges[u].size(); i++) {
        int v = edges[u][i];
        if (!sameComponent(u, v) &&
            cut[u] < cut[v]) {  // inactive means it crosses a cut,
                                // sorting so we don't double count
          auto p = std::make_pair(cut[u], cut[v]);

          if (edgeWeights.find(p) != edgeWeights.end())
            edgeWeights[p] += getEdgeWeight(u, i);
          else
            edgeWeights[p] = getEdgeWeight(u, i);
        }
      }
    }

    VLOG(8) << "Finished computing the edge weights.";

    for (const auto &[e, w] : edgeWeights) {
      assert(w > 0 && "Edge weights must be positive");

      newEdges.push_back(Edge(e.first, e.second));
      newWeights.push_back(w);
    }

    assert(newEdges.size() == newWeights.size() &&
           "Each edge must have a weight.");

    return std::make_shared<TemplateGraph<V, E, true, WeightType>>(
        numCutComponents, newEdges, newWeights, componentWeights);
  }

  /**
   * Depth first search on the graph to label all connected components.
   */
  template <typename iterType>
  DFSResult dfs(iterType begin, iterType end) {
    // TODO: make this space efficient
    std::vector<bool> visited(size(), false);
    std::vector<int> labels(size(), -1);
    std::vector<int> sizes;

    VLOG(8) << "Running DFS on component of size " << std::distance(begin, end);

    std::vector<int> queue;

    int label = 0;
    int size = 0;

    for (auto it = begin; it != end; it++) {
      auto u = *it;

      if (!visited[u]) {
        queue.push_back(u);

        while (!queue.empty()) {
          auto v = queue.back();
          queue.pop_back();

          if (!visited[v]) {
            visited[v] = true;
            labels[v] = label;
            size++;

            for (auto e = cbeginEdges(v); e != cendEdges(v); e++) {
              if (sameComponent(v, *e)) {
                queue.push_back(*e);
              }
            }
          }
        }
        assert(queue.empty() && "Queue must be empty.");

        label++;

        sizes.push_back(size);
        size = 0;
      }
    }

    VLOG(8) << "Found " << label << " connected components.";

    return DFSResult(label, std::move(labels), std::move(sizes));
  }

  /**
   * Run DFS on the whole graph. TODO: make more efficient, the iterator just
   * iterates on the range 1 to n, don't need the allocation.
   */
  DFSResult dfs() {
    std::vector<int> nodes(size());
    std::iota(nodes.begin(), nodes.end(), 0);
    return dfs(nodes.begin(), nodes.end());
  }

  // /**
  //  * Compute the boundary vertices of the graph.
  //  */
  // std::unordered_set<int> &getBoundary() {
  //   for (int u = 0; u < size(); u++) {
  //     if (degrees[u] != initialDegrees[u]) boundary.insert(u);
  //   }
  //   return boundary;
  // }

  /**
   * Compute the degree information for the refinement process and create a
   * vector of boundary vertices.
   */
  // std::vector<V> computeRefinementDegrees() {
  //   std::vector<V> boundaryVertices;

  //   std::vector<W> externalDegrees(size(), 0), internalDegrees(size(), 0);

  //   // VLOG(0) << "size: " << size();

  //   for (int u = 0; u < size(); u++) {
  //     for (int j = start[u]; j < start[u + 1]; j++) {
  //       if (cluster[u] != cluster[edges[j]]) {
  //         externalDegrees[u] += getEdgeWeight(j);
  //       } else {
  //         internalDegrees[u] += getEdgeWeight(j);
  //       }
  //     }

  //     // int gain = externalDegrees[u] - internalDegrees[u];
  //     // VLOG(0) << "gain of vertex " << u << ": " << gain;

  //     if (int gain = externalDegrees[u] - internalDegrees[u]; gain >= 0) {
  //       boundaryVertices.push_back(u);
  //     }
  //   }

  //   VLOG(0) << "Found " << boundaryVertices.size() << " boundary vertices.";

  //   return boundaryVertices;
  // }

  // /**
  //  * Compute the gain of a vertex.
  //  */
  // W computeGain(V u) {
  //   std::vector<W> gains(numClusters, 0);

  //   // VLOG(0) << "Computing gains for vertex " << u;
  //   // VLOG(0) << "numClusters: " << numClusters;

  //   for (int j = start[u]; j < start[u + 1]; j++) {
  //     V v = edges[j];
  //     int vIndex = cluster[v];

  //     gains[vIndex] += getEdgeWeight(j);
  //   }

  //   // VLOG(0) << "uwu";

  //   int maxGain = std::numeric_limits<int>::min();

  //   for (size_t i = 0; i < gains.size(); i++) {
  //     if (gains[i] > maxGain) {
  //       maxGain = gains[i];
  //     }
  //   }

  //   return maxGain - gains[cluster[u]];
  // }

  // int moveWhere(V u) {
  //   int to;

  //   std::vector<W> gains(numClusters, 0);

  //   for (int j = start[u]; j < start[u + 1]; j++) {
  //     V v = edges[j];
  //     int vIndex = cluster[v];

  //     gains[vIndex] += getEdgeWeight(j);
  //   }

  //   int maxGain = std::numeric_limits<int>::min();

  //   for (size_t i = 0; i < gains.size(); i++) {
  //     if (gains[i] > maxGain) {
  //       maxGain = gains[i];
  //       to = i;
  //     }
  //   }
  //   return to;
  // }

  // /**
  //  * Create a priority queue of vertices to be moved.
  //  */
  // BucketQueue<V, W> initializeQueue() {
  //   BucketQueue<V, W> refinementQueue;

  //   auto boundaryVector = computeRefinementDegrees();

  //   VLOG(0) << "Created boundary vector";

  //   for (auto v : boundaryVector) {
  //     auto gain = computeGain(v);
  //     refinementQueue.insert(v, gain);
  //   }

  //   return refinementQueue;
  // }

  // W gainTo(V u, int to) {
  //   W gain;

  //   for (int j = start[u]; j < start[u + 1]; j++) {
  //     auto v = edges[j];

  //     if (cluster[v] == to) {
  //       gain += getEdgeWeight(j);
  //     }
  //   }

  //   return gain;
  // }

  // /**
  //  * Perform a single move.
  //  */
  // bool makeSingleMove(V u, BucketQueue<V, W> &pQueue,
  //                     std::unordered_set<int> &moved) {
  //   int where = cluster[u];
  //   int to = moveWhere(u);

  //   // if the vertex is the only one in its partition don't move
  //   if (clusterSizes[where] == 1) {
  //     return false;
  //   }

  //   /**
  //    * TODO: think of a condition here that prevents empty nodes
  //    */

  //   // auto gain = gainTo(u, to);

  //   cluster[u] = to;
  //   moved.insert(u);

  //   clusterSizes[where]--;
  //   clusterSizes[to]++;

  //   // move vertex to this partition

  //   // update all neighbors, insert new boundary vertices into queue and
  //   // update gains for those already inside. maybe also delete some
  //   // non-boundary vertices.
  //   for (int j = start[u]; j < start[u + 1]; j++) {
  //     auto v = edges[j];
  //     auto gain = computeGain(v);

  //     if (moved.count(v)) continue;

  //     if (gain < 0) {
  //       pQueue.removeIfPresent(v);
  //     } else {
  //       pQueue.updateOrInsert(v, gain);
  //     }
  //   }

  //   return true;
  // }

  // /**
  //  * Compute the normalized cut of the current cluster set.
  //  */
  // double computeNCut() {
  //   DegreeVector clusterCuts(numClusters, 0), clusterVolumes(numClusters, 0);

  //   for (int u = 0; u < size(); u++) {
  //     clusterVolumes[cluster[u]] += getVertexWeight(u);

  //     for (int j = start[u]; j < start[u + 1]; j++) {
  //       auto v = edges[j];

  //       if (cluster[u] != cluster[v]) {
  //         clusterCuts[cluster[u]] += getEdgeWeight(j);
  //       }
  //     }
  //   }

  //   double val = 0.0;

  //   for (int i = 0; i < numClusters; i++) {
  //     val += (double)clusterCuts[i] / (double)clusterVolumes[i];
  //   }

  //   return val;
  // }

  // /**
  //  * Perform a single round of refinement.
  //  */
  // void refine() {
  //   VLOG(0) << "Beginning refinement.";
  //   // Initialize priority queue
  //   auto refinementQueue = initializeQueue();
  //   VLOG(0) << "Created queue.";

  //   bool stop = false;
  //   double target = computeNCut();

  //   std::unordered_set<int> moved;

  //   VLOG(0) << "Current target: " << target;

  //   int best_index = -1;
  //   int i = 0;

  //   std::vector<V> moves;
  //   std::vector<int> from;

  //   while (!stop) {
  //     if (refinementQueue.empty()) {
  //       break;
  //     }

  //     auto top = refinementQueue.popTop();
  //     auto topIndex = cluster[top];

  //     VLOG(0) << "Moving " << top;

  //     auto success = makeSingleMove(top, refinementQueue, moved);

  //     // if the move couldn't be made, don't need to do anythig else
  //     if (!success) {
  //       continue;
  //     }

  //     moves.push_back(top);
  //     from.push_back(topIndex);

  //     auto current = computeNCut();

  //     VLOG(0) << "Cost after move: " << current;

  //     if (current < target) {
  //       best_index = i;
  //       target = current;
  //     }

  //     if (i > best_index + 50) {
  //       stop = true;
  //     }

  //     i++;
  //   }

  //   VLOG(0) << "Concluded refinement after " << i
  //           << " steps. New objective value: " << target;

  //   VLOG(0) << "cluster sizes:";
  //   for (auto sz : clusterSizes) {
  //     VLOG(0) << sz;
  //   }

  //   // undo moves
  //   for (i--; i > best_index; i--) {
  //     auto v = moves.back();
  //     auto to = from.back();

  //     moves.pop_back();
  //     from.pop_back();

  //     cluster[v] = to;
  //   }

  //   VLOG(0) << "Target after undoing: " << computeNCut();
  //   // fine

  //   // TODO (?): clean up the data structure
  // }

  // /**
  //  * Refine multiple times
  //  */
  // void todo() {}

  /**
   * Add a source and sink to the
   */

  template <bool w = isWeighted>
  std::enable_if_t<!w, std::pair<V, V>> addSourceSink(std::vector<V> &ss,
                                                      std::vector<V> &ts) {
    auto s = size_;
    auto t = size_ + 1;

    reverse.emplace_back();

    for (int i = 0; i < ss.size(); i++) {
      auto u = ss[i];

      reverse[s].push_back(edges[u].size());

      edges[u].push_back(s);
      reverse[u].push_back(i);
    }

    reverse.emplace_back();

    for (int i = 0; i < ts.size(); i++) {
      auto u = ts[i];

      reverse[t].push_back(edges[u].size());

      edges[u].push_back(t);
      reverse[u].push_back(i);
    }

    edges.push_back(std::move(ss));
    edges.push_back(std::move(ts));

    assert(size_ + 2 == size());
    assert(edges.size() == size());
    assert(reverse.size() == size());

    return std::make_pair(s, t);
  }

  template <bool w = isWeighted>
  std::enable_if_t<w, std::pair<V, V>> addSourceSink(std::vector<V> &ss,
                                                     std::vector<V> &ts) {
    auto s = size_;
    auto t = size_ + 1;

    reverse.emplace_back();
    wm::weights.emplace_back();

    for (int i = 0; i < ss.size(); i++) {
      auto u = ss[i];

      reverse[s].push_back(edges[u].size());
      wm::weights[s].push_back(initialDegrees[u]);

      edges[u].push_back(s);
      reverse[u].push_back(i);
      wm::weights[u].push_back(initialDegrees[u]);
    }

    reverse.emplace_back();
    wm::weights.emplace_back();

    for (int i = 0; i < ts.size(); i++) {
      auto u = ts[i];

      reverse[t].push_back(edges[u].size());
      wm::weights[t].push_back(initialDegrees[u]);

      edges[u].push_back(t);
      reverse[u].push_back(i);
      wm::weights[u].push_back(initialDegrees[u]);
    }

    edges.push_back(std::move(ss));
    edges.push_back(std::move(ts));

    // VLOG(0) << size_;
    // VLOG(0) << size();
    assert(size_ + 2 == size());
    assert(edges.size() == size());
    assert(reverse.size() == size());

    return std::make_pair(s, t);
  }

  template <bool w = isWeighted>
  std::enable_if_t<!w> removeSourceSink() {
    auto &ts = edges.back();

    for (auto u : ts) {
      edges[u].pop_back();
      reverse[u].pop_back();
    }

    edges.pop_back();
    reverse.pop_back();

    auto &ss = edges.back();

    for (auto u : ss) {
      edges[u].pop_back();
      reverse[u].pop_back();
    }

    edges.pop_back();
    reverse.pop_back();
  }

  template <bool w = isWeighted>
  std::enable_if_t<w> removeSourceSink() {
    auto &ts = edges.back();

    for (auto u : ts) {
      edges[u].pop_back();
      reverse[u].pop_back();
      wm::weights[u].pop_back();
    }

    edges.pop_back();
    reverse.pop_back();
    wm::weights.pop_back();

    auto &ss = edges.back();

    for (auto u : ss) {
      edges[u].pop_back();
      reverse[u].pop_back();
      wm::weights[u].pop_back();
    }

    edges.pop_back();
    reverse.pop_back();
    wm::weights.pop_back();
  }
};

using Graph = TemplateGraph<int, Edge, false>;
using WeightedGraph = TemplateGraph<int, Edge, true>;

}  // namespace DynamicGraph

#endif  // DYNAMIC_GRAPH_H_
