#ifndef ARRAY_GRAPH_H_
#define ARRAY_GRAPH_H_

#include <glog/logging.h>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "lib/datastructures/bucket_queue.hpp"

namespace ArrayGraph {

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
  std::vector<W> weights;
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
  WeightType vertexWeight;

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
   * Begin of list of neighbors in the array. It indicates the starting point of
   * the neighbors in the edge arary. The neighbors end at start[u+1].
   *
   * This must contain n+1 elements, where the last one in the list is equal to
   * m.
   */
  std::vector<int> start;

  /**
   * Array holding the end of the *active* edges. All edges connecting between u
   * and its neighbors in the same component are in the range [start[u],
   * end[u]).
   */
  std::vector<int> endActive;

  /**
   * Index of cut of expander decomposition that a given vertex belongs to.
   */
  std::vector<int> cut;

  int dfsCalls;
  std::vector<int> dfsVisited;
  std::vector<int> dfsStack;

  int numClusters;

  /**
   * Cluster that the vertex belongs to. Relevant for the refinement.
   */
  std::vector<int> cluster;

  std::vector<int> clusterSizes;

  std::unordered_set<int> tabooVertices;

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
  std::vector<V> edges;

  /**
   * Pointer to the reverse edge.
   */
  // std::vector<int> reverse;

  // these are for internal use, we don't need a client to know about the edge
  // ids
  V getEdge(int i) const { return edges[i]; }

  template <bool t = isWeighted>
  std::enable_if_t<t, WeightType> getEdgeWeight(int i) const {
    return wm::weights[i];
  }

  template <bool t = isWeighted>
  std::enable_if_t<!t, WeightType> getEdgeWeight(int i) const {
    return 1;
  }

  /**
   * Get the index of the reverse edge.
   */
  // int getReverse(int j) const { return reverse[j]; };

 public:
  /**
   * Construct an empty graph.
   */
  TemplateGraph() : size_(0), volume(0), numCutComponents(0) {}

  TemplateGraph(int n, const std::vector<E> &es)
      : size_(n),
        volume(2 * es.size()),
        vertexWeight(n),
        numCutComponents(1),
        volumes(1, 2 * es.size()),
        componentWeights(1, 2 * es.size()),
        degrees(n, 0),
        initialDegrees(),
        start(n, 0),
        endActive(n, 0),
        cut(n, 0),
        dfsVisited(n, 0),
        dfsCalls(1),
        cluster(n, 0),
        boundary(),
        edges(2 * es.size())
  // reverse(2 * es.size())  {
  {
    static_assert(!isWeighted);

    for (auto e : es) {
      degrees[e.from]++;
      degrees[e.to]++;
    }

    initialDegrees = degrees;

    int sum = 0;
    for (int i = 0; i < size_; i++) {
      start[i] = sum;
      sum += degrees[i];
      endActive[i] = sum;
    }
    start.push_back(
        sum);  // add a final element that designates the end of the vector

    assert(sum == volume && "Number of edges must match");

    std::vector<int> counts(start);

    // fill the edge array
    for (auto e : es) {
      // reverse[counts[e.to]] = counts[e.from];
      // reverse[counts[e.from]] = counts[e.to];

      edges[counts[e.to]++] = e.from;
      edges[counts[e.from]++] = e.to;
    }

    assert(counts.back() == sum && "Edges must be inserted correctly.");

    // for (int u = 0; u < size(); u++) {
    //   for (int i = start[u]; i < start[u + 1]; i++) {
    //     assert(u == edges[reverse[i]] && "Reverse edge must point back.");
    //   }
    // }

    // // since we currently only work with connected graphs
    // assert(std::all_of(degrees.begin(), degrees.end(),
    //                    [](int a) { return a > 0; }) &&
    //        "No degree zero vertices.");
  }

  /**
   * Alternative constructor
   */

  TemplateGraph(int n, int m, std::vector<int> edges, std::vector<int> degrees,
                std::vector<int> start)
      : size_(n),
        volume(2 * m),
        vertexWeight(n),
        numCutComponents(1),
        volumes(1, 2 * m),
        componentWeights(1, 2 * m),
        degrees(degrees),
        initialDegrees(std::move(degrees)),
        start(std::move(start)),
        endActive(n, 0),
        cut(n, 0),
        dfsVisited(n, 0),
        dfsCalls(1),
        cluster(n, 0),
        boundary(),
        edges(std::move(edges))
  // reverse(2 * m) {
  {
    static_assert(!isWeighted);
    VLOG(4) << "deg size: " << this->degrees.size();
    VLOG(2) << "n: " << n;

    assert(this->degrees.size() == static_cast<size_t>(n) &&
           "Each node must have a degree.");
    assert(std::reduce(this->degrees.begin(), this->degrees.end(), 0) ==
               2 * m &&
           "Sum of degrees must equal volume.");

    // compute the reverse edge indices
    for (int u = 0; u < n; u++) {
      // for (int j = this->start[u]; j < this->start[u + 1]; j++) {
      //   auto v = this->edges[j];

      //   // find the index of the reverse edge
      //   int i = this->start[v];
      //   while (this->edges[i] != u) {
      //     ++i;
      //   }

      //   reverse[j] = i;
      // }

      endActive[u] = this->start[u + 1];
    }

    // for (int u = 0; u < size(); u++) {
    //   for (int i = this->start[u]; i < this->start[u + 1]; i++) {
    //     assert(u == this->edges[reverse[i]] && "Reverse edge must point
    //     back.");
    //   }
    // }
  }

  /**
   * Weighted graph constructor.
   */
  TemplateGraph(int n, const std::vector<E> &es, const std::vector<W> &ws,
                const std::vector<W> &vws)
      : size_(n),
        volume(2 * std::accumulate(ws.begin(), ws.end(), 0)),
        vertexWeight(std::accumulate(vws.begin(), vws.end(), 0)),
        numCutComponents(1),
        volumes(1, 2 * std::accumulate(ws.begin(), ws.end(), 0)),
        componentWeights(1, std::accumulate(vws.begin(), vws.end(), 0)),
        degrees(n, 0),
        initialDegrees(n, 0),
        start(n, 0),
        endActive(n, 0),
        cut(n, 0),
        dfsVisited(n, 0),
        dfsCalls(1),
        cluster(n, 0),
        boundary(),
        edges(2 * es.size())
  // reverse(2 * es.size()) {
  {
    static_assert(isWeighted);
    assert(ws.size() == es.size() && "Number of edges and weights must match.");

    wm::weights = std::vector<W>(2 * es.size());
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

    size_t sum = 0;
    for (int i = 0; i < size_; i++) {
      start[i] = sum;
      sum += numNeighbors[i];
      endActive[i] = sum;
    }
    start.push_back(
        sum);  // add a final element that designates the end of the vector

    assert(sum == 2 * es.size() && "Number of edges must match");

    std::vector<int> counts(start);

    // fill the edge array
    for (size_t i = 0; i < es.size(); i++) {
      auto e = es[i];
      auto w = ws[i];

      wm::weights[counts[e.to]] = w;
      wm::weights[counts[e.from]] = w;

      // reverse[counts[e.to]] = counts[e.from];
      // reverse[counts[e.from]] = counts[e.to];

      assert(e.from != e.to && "Self loops are not allowed.");

      edges[counts[e.to]++] = e.from;
      edges[counts[e.from]++] = e.to;
    }

    assert(static_cast<size_t>(counts.back()) == sum &&
           "Edges must be inserted correctly.");

    // for (int u = 0; u < size(); u++) {
    //   for (int i = start[u]; i < start[u + 1]; i++) {
    //     assert(u == edges[reverse[i]] && "Reverse edge must point back.");
    //   }
    // }
  }

  /**
   * Alternative constructor
   */

  TemplateGraph(int n, int m, std::vector<int> edges, std::vector<int> weights,
                std::vector<int> degrees, std::vector<int> start,
                std::vector<int> vertexWeights)
      : size_(n),
        volume(std::accumulate(degrees.begin(), degrees.end(), 0)),
        vertexWeight(
            std::accumulate(vertexWeights.begin(), vertexWeights.end(), 0)),
        numCutComponents(1),
        volumes(1, std::accumulate(weights.begin(), weights.end(), 0)),
        componentWeights(
            1, std::accumulate(vertexWeights.begin(), vertexWeights.end(), 0)),
        degrees(degrees),
        initialDegrees(std::move(degrees)),
        start(std::move(start)),
        endActive(n, 0),
        cut(n, 0),
        dfsVisited(n, 0),
        dfsCalls(1),
        cluster(n, 0),
        boundary(),
        edges(std::move(edges))
  // reverse(2 * m) { {
  {
    static_assert(isWeighted);

    wm::weights = std::move(weights);
    wm::vertexWeights = std::move(vertexWeights);

    assert(this->degrees.size() == static_cast<size_t>(n) &&
           "Each node must have a degree.");

    // compute the reverse edge indices
    for (int u = 0; u < n; u++) {
      // VLOG(0) << "u: " << u;
      // for (int j = this->start[u]; j < this->start[u + 1]; j++) {
      //   auto v = this->edges[j];

      //   // find the index of the reverse edge
      //   int i = this->start[v];
      //   while (this->edges[i] != u) {
      //     ++i;
      //   }

      //   reverse[j] = i;
      // }

      endActive[u] = this->start[u + 1];
    }

    // for (int u = 0; u < size(); u++) {
    //   for (int i = this->start[u]; i < this->start[u + 1]; i++) {
    //     assert(u == this->edges[reverse[i]] && "Reverse edge must point
    //     back.");
    //   }
    // }
  }

  /**
   * Edge begin and end operator. Iterates over all edges.
   */
  EdgeIt begin() { return edges.begin(); }
  EdgeIt end() { return edges.end(); }

  /**
   * Constant edge begin and end operator. Iterates over all edges.
   */
  cEdgeIt cbegin() const { return edges.cbegin(); }
  cEdgeIt cend() const { return edges.cend(); }

  /**
   * Iterate the neighbors of vertex u.
   */
  EdgeIt beginEdges(V u) { return edges.begin() + start[u]; }
  EdgeIt endActiveEdges(V u) { return edges.begin() + endActive[u]; }
  EdgeIt endEdges(V u) { return edges.begin() + start[u + 1]; }

  /**
   * Const-iterate the neighbors of vertex u.
   */
  cEdgeIt cbeginEdges(V u) const { return edges.cbegin() + start[u]; }
  cEdgeIt cendActiveEdges(V u) const { return edges.cbegin() + endActive[u]; }
  cEdgeIt cendEdges(V u) const { return edges.cbegin() + start[u + 1]; }

  /**
   * Iterate the weights of neighbors of vertex u.
   */
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::iterator> beginWeights(V u) {
    return wm::weights.begin() + start[u];
  }

  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::iterator> endWeights(V u) {
    return wm::weights.begin() + start[u + 1];
  }

  /**
   * Const-iterate the weights of neighbors of vertex u.
   */
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::const_iterator> cbeginWeights(
      V u) const {
    return wm::weights.cbegin() + start[u];
  }
  template <bool w = isWeighted>
  std::enable_if_t<w, typename std::vector<W>::const_iterator> cendWeigths(
      V u) const {
    return wm::weights.cbegin() + start[u + 1];
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
  int size() const { return size_; }

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

  /**
   * Get the number of neighboring vertices.
   */
  int getNumNeighbors(V u) const { return start[u + 1] - start[u]; }

  /**
   * Get the number of active neighbors.
   */
  int getNumActiveNeighbors(V u) const { return start[u] - endActive[u]; }

  V getEdge(V u, int i) const { return edges[start[u] + i]; }

  template <bool t = isWeighted>
  std::enable_if_t<t, WeightType> getEdgeWeight(V u, int i) const {
    return wm::weights[start[u] + i];
  }

  template <bool t = isWeighted>
  std::enable_if_t<!t, WeightType> getEdgeWeight(V u, int i) const {
    return 1;
  }

  // /**
  //  * Get the index of the reverse edge.
  //  */
  // int getReverse(V u, int j) const {
  //   auto v = edges[start[u] + j];

  //   return reverse[start[u] + j] - start[v];
  // };

  /**
   * Returns a const reference to the current cut array.
   */
  std::vector<int> const &getCutArray() const { return cut; }

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

    /**
     * TODO: need to fix the endactive part here
     */
    LOG(ERROR) << "setCutArray is not currently implemented correctly.";

    for (int u = 0; u < size(); u++) {
      volumes[cut[u]] += initialDegrees[u];
      componentWeights[cut[u]] += getVertexWeight(u);

      for (int j = start[u]; j < start[u + 1]; j++) {
        if (!sameComponent(u, edges[j])) {
          degrees[u] -= getEdgeWeight(j);
        }
      }
    }
  }

  /**
   * Reset the cut array so a new expander decomposition can be computed.
   */
  void resetCutArray() {
    numCutComponents = 1;

    cut = std::vector<int>(size(), 0);

    volumes = DegreeVector(1, volume);
    componentWeights = DegreeVector(1, vertexWeight);
    degrees = initialDegrees;
    for (int i = 0; i < endActive.size(); i++) {
      endActive[i] = start[i + 1];
    }
  }

  const std::vector<int> &getInitialDegreeArray() const {
    return initialDegrees;
  }

  WeightType getCutEdges(V u) const { return initialDegrees[u] - degrees[u]; }

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
  int getNumEdges() const { return edges.size() / 2; }

  /**
   * Check whether two vertices are in the same component.
   */
  bool sameComponent(int u, int v) const { return cut[u] == cut[v]; }

  /**
   * Check whether two vertices are in the same cluster.
   */
  bool sameCluster(int u, int v) const { return cluster[u] == cluster[v]; }

  // start DFS from u and give all encountered vertices the label label
  // also update the weight and volume of the new component
  // return the number of vertices in the component
  int dfsLoop(int u, int componentLabel, int newLabel) {
    int numVertices = 0;

    int newWeight = 0;
    int newVolume = 0;

    dfsStack.push_back(u);

    while (!dfsStack.empty()) {
      auto v = dfsStack.back();
      dfsStack.pop_back();

      if (dfsVisited[v] < dfsCalls) {
        dfsVisited[v] = dfsCalls;

        cut[v] = newLabel;  // BUG HERE

        newVolume += initialDegrees[v];
        newWeight += getVertexWeight(v);
        numVertices++;

        for (auto e = cbeginEdges(v); e != cendEdges(v); e++) {
          if (cut[*e] == componentLabel) {
            dfsStack.push_back(*e);
          }
        }
      }
    }

    // finished the dfs starting from u
    // can update the vertices now
    componentWeights[newLabel] = newWeight;
    volumes[newLabel] = newVolume;

    assert(dfsStack.empty() && "DFS stack must be empty");

    return numVertices;
  }

  template <typename iterType>
  void fixDegrees(iterType begin, iterType end) {
    for (iterType it = begin; it != end; it++) {
      int u = *it;

      WeightType newDegree = 0;

      int l = start[u];
      int r = endActive[u] - 1;

      while (r >= l) {
        auto v = edges[l];

        if (!sameComponent(u, v)) {
          // need to move the edge to the end
          std::swap(edges[l], edges[r]);

          if constexpr (Weighted) {
            std::swap(wm::weights[l], wm::weights[r]);
          }

          r--;
        } else {
          // edge can stay
          newDegree += getEdgeWeight(l);
          l++;
        }
      }
      endActive[u] = l;

      degrees[u] = newDegree;
    }
  }

  template <typename iterType>
  std::pair<std::vector<int>, std::vector<int>> makeNewComponentsFromCut(
      iterType begin, iterType mid, iterType end) {
    // make a new component
    VLOG(2) << "Making new component be component " << numCutComponents;

    int leftIndex = cut[*begin];
    int rightIndex = numCutComponents;

    for (iterType it = mid; it != end; ++it) {
      cut[*it] = numCutComponents;
    }

    volumes.push_back(0);
    componentWeights.push_back(0);

    numCutComponents++;

    // now fix the component information, doing dfs
    std::vector<int> firstSizes, secondSizes;

    auto size = dfsLoop(*begin, leftIndex, leftIndex);

    VLOG(2) << "Ran DFS on left hand side. Connected component is of size "
            << size;

    firstSizes.push_back(size);

    for (iterType it = begin; it != mid; ++it) {
      // not visited yet in this iteration
      if (dfsVisited[*it] < dfsCalls) {
        // this idiom needs its own function
        volumes.push_back(0);
        componentWeights.push_back(0);

        size = dfsLoop(*it, leftIndex, numCutComponents++);
        firstSizes.push_back(size);
      }
    }

    size = dfsLoop(*mid, rightIndex, rightIndex);

    VLOG(2) << "Ran DFS on right hand side. Connected component is of size "
            << size;

    secondSizes.push_back(size);

    for (iterType it = mid; it != end; ++it) {
      if (dfsVisited[*it] < dfsCalls) {
        volumes.push_back(0);
        componentWeights.push_back(0);

        size = dfsLoop(*it, rightIndex, numCutComponents++);
        secondSizes.push_back(size);
      }
    }

    assert(std::reduce(volumes.begin(), volumes.end(), 0) == volume &&
           "No volume can be lost.");

    fixDegrees(begin, end);

    VLOG(2) << "Found " << firstSizes.size()
            << " components in the first part.";
    VLOG(2) << "Found " << secondSizes.size()
            << " components in the second part.";
    VLOG(2) << "Num Components: " << numCutComponents;

    dfsCalls++;

    return std::make_pair(firstSizes, secondSizes);
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
   * TODO: rewrite this function
   */
  TemplateGraph<V, E, true, WeightType> *contractedGraph() {
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
      return new TemplateGraph<V, E, true, WeightType>(
          numCutComponents, newEdges, newWeights, componentWeights);
    }

    VLOG(8) << "Contracting the graph, not a singleton.";

    for (int u = 0; u < size(); u++) {
      for (int i = start.at(u); i < start.at(u + 1); i++) {
        int v = edges.at(i);
        if (!sameComponent(u, v) &&
            cut[u] < cut[v]) {  // inactive means it crosses a cut,
                                // sorting so we don't double count
          auto p = std::make_pair(cut[u], cut[v]);

          if (edgeWeights.find(p) != edgeWeights.end())
            edgeWeights[p] += getEdgeWeight(i);
          else
            edgeWeights[p] = getEdgeWeight(i);
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

    return new TemplateGraph<V, E, true, WeightType>(
        numCutComponents, newEdges, newWeights, componentWeights);
  }

  /**
   * Faster contraction method
   */
  TemplateGraph<V, E, true, WeightType> *contractedGraphArray() {
    // get the nodes sorted by component
    auto perm = std::vector<int>(size(), 0);
    std::iota(perm.begin(), perm.end(), 0);

    std::sort(perm.begin(), perm.end(),
              [this](int u, int v) { return cut[u] < cut[v]; });

    auto edgeTargets = std::vector<int>(numCutComponents, -1);
    auto targetPosition = std::vector<int>(numCutComponents, 0);

    auto newEdgeArray = std::vector<int>(edges.size());
    auto newWeightArray = std::vector<WeightType>(edges.size());
    auto newDegreeArray = std::vector<int>(numCutComponents, 0);
    auto newStartArray = std::vector<int>(numCutComponents, 0);
    auto newVertexWeightArray = std::vector<WeightType>(numCutComponents, 0);

    int currentComponent = 0, edgeCounter = 0;

    for (auto v : perm) {
      // new component
      if (currentComponent < cut[v]) {
        // components must be sequential and non-empty
        assert(cut[v] == currentComponent + 1);

        currentComponent++;
        newStartArray[currentComponent] = edgeCounter;
      }

      newVertexWeightArray[currentComponent] += getVertexWeight(v);

      // iterate over possible edges
      for (auto j = endActive[v]; j != start[v + 1]; j++) {
        auto u = edges[j];

        auto target = cut[u];

        // can only be edges outside of component
        assert(cut[u] != cut[v]);

        newDegreeArray[currentComponent] += getEdgeWeight(j);

        // check if we have seen this edge target already
        if (edgeTargets[target] < currentComponent) {
          // if not, insert it in the array
          newEdgeArray[edgeCounter] = target;
          newWeightArray[edgeCounter] = getEdgeWeight(j);

          edgeTargets[target] = currentComponent;
          targetPosition[target] = edgeCounter;

          edgeCounter++;
        } else {
          newWeightArray[targetPosition[target]] += getEdgeWeight(j);
        }
      }
    }

    newStartArray.push_back(edgeCounter);

    assert(newStartArray.size() == numCutComponents + 1);

    newEdgeArray.resize(edgeCounter);
    newWeightArray.resize(edgeCounter);

    assert(newEdgeArray.size() % 2 == 0 && "Must be even.");

    return new TemplateGraph<V, E, true, WeightType>(
        numCutComponents, newEdgeArray.size() / 2, std::move(newEdgeArray),
        std::move(newWeightArray), std::move(newDegreeArray),
        std::move(newStartArray), std::move(newVertexWeightArray));
  }

  // TemplateGraph<V, E, true, WeightType> *contractedGraphHash() {}

  /**
   * Depth first search on the graph to label all connected components.
   */
  template <typename iterType>
  DFSResult dfs(iterType begin, iterType end) {
    /**
     * TODO: replace this function with one based on dfsLoop
     */

    // getting ready to delete this
    // assert(numCutComponents == 1 && "Only using this in setup.");

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
   *
   *
   *
   *
   *  Everything below here is refinement stuff
   *
   *
   *
   *
   */
  DegreeVector clusterCuts, clusterVolumes;

  void updateAfterMove(int u, int from, int to) {
    assert(tabooVertices.count(u) == 0);

    clusterSizes[from]--;
    clusterSizes[to]++;

    // now fix volumes
    clusterVolumes[from] -= getVertexWeight(u);
    clusterVolumes[to] += getVertexWeight(u);

    // now fix edges
    for (int j = start[u]; j < start[u + 1]; j++) {
      auto v = edges[j];

      if (cluster[v] == from) {
        // edges that were internal before, now cross between clusters
        clusterCuts[from] += getEdgeWeight(j);
        clusterCuts[to] += getEdgeWeight(j);
      } else if (cluster[v] == to) {
        // edges that become internal
        clusterCuts[from] -= getEdgeWeight(j);
        clusterCuts[to] -= getEdgeWeight(j);
      } else {
        // edges that stay external
        clusterCuts[from] -= getEdgeWeight(j);
        clusterCuts[to] += getEdgeWeight(j);
      }
    }
  }

  /**
   * Returns a const reference to the current cluster array.
   */
  std::vector<int> const &getClusterArray() const { return cluster; }

  int getNumClusters() const { return numClusters; }

  void setCluster(V u, int clusterLabel) {
    auto where = cluster[u];

    cluster[u] = clusterLabel;

    updateAfterMove(u, where, clusterLabel);

    assert(clusterSizes[where] > 0 && "Cluster must not become empty.");
  }

  void newCluster(V u) {
    auto where = cluster[u];
    auto clusterLabel = numClusters;

    cluster[u] = numClusters;
    clusterSizes.push_back(0);
    clusterVolumes.push_back(0);
    clusterCuts.push_back(0);

    numClusters++;

    updateAfterMove(u, where, clusterLabel);

    assert(clusterSizes[where] > 0 && "Cluster must not become empty.");
  }

  /**
   * Makes a vertex taboo. This means it won't count towards cluster size, which
   * should prevent it from being the only node in a cluster.
   */
  void makeTaboo(V u) {
    tabooVertices.insert(u);

    auto where = cluster[u];

    clusterSizes[where]--;

    assert(clusterSizes[where] > 0 &&
           "Taboo node can not be alone in cluster.");
  }

  void setClustering(int n, std::vector<int> const &parentArray) {
    assert(parentArray.size() == static_cast<size_t>(getNumComponents()) &&
           "Each component exactly one parent.");

    numClusters = n;
    clusterSizes = std::vector<int>(n, 0);
    clusterCuts = std::vector<W>(n, 0);
    clusterVolumes = std::vector<W>(n, 0);

    for (int i = 0; i < size(); i++) {
      assert(parentArray[cut[i]] < numClusters &&
             "Num clusters must be correct.");

      cluster[i] = parentArray[cut[i]];

      clusterSizes[cluster[i]]++;
    }

    // now fix the clusters
    clusterCuts = std::vector<W>(n, 0);
    clusterVolumes = std::vector<W>(n, 0);

    for (int u = 0; u < size(); u++) {
      auto where = cluster[u];

      clusterVolumes[where] += getVertexWeight(u);

      for (int j = start[u]; j < start[u + 1]; j++) {
        auto v = edges[j];
        if (cluster[u] != cluster[v]) {
          clusterCuts[where] += getEdgeWeight(j);
        }
      }
    }

    for (auto s : clusterSizes) {
      assert(s > 0 && "Cluster size may not be zero after projecting.");
    }

    // printDebugRefine();
  }

  /**
   * Compute the degree information for the refinement process and create a
   * vector of boundary vertices.
   */
  std::vector<V> computeRefinementDegrees() {
    std::vector<V> boundaryVertices;

    std::vector<W> externalDegrees(size(), 0), internalDegrees(size(), 0);

    for (int u = 0; u < size(); u++) {
      for (int j = start[u]; j < start[u + 1]; j++) {
        auto v = edges[j];

        if (cluster[u] != cluster[v]) {
          // intra cluster edge
          externalDegrees[u] += getEdgeWeight(j);
        } else {
          // inter cluster edge
          internalDegrees[u] += getEdgeWeight(j);
        }
      }

      // insert those that have positive gain
      if (int gain = externalDegrees[u] - internalDegrees[u]; gain >= 0) {
        boundaryVertices.push_back(u);
      } else if (true) {
        // additionally see what happens if we insert those that are on boundary
        // and in large clusters
        if (externalDegrees[u] > 0 &&
            clusterVolumes[cluster[u]] >= vertexWeight / numClusters) {
          boundaryVertices.push_back(u);
        }
      }
    }

    return boundaryVertices;
  }

  /**
   * Compute the gain of a vertex.
   */
  W computeGain(V u) {
    std::vector<W> gains(numClusters, 0);

    for (int j = start[u]; j < start[u + 1]; j++) {
      V v = edges[j];
      int vIndex = cluster[v];

      gains[vIndex] += getEdgeWeight(j);
    }

    int maxGain = std::numeric_limits<int>::min();

    for (size_t i = 0; i < gains.size(); i++) {
      if (gains[i] > maxGain) {
        maxGain = gains[i];
      }
    }

    return maxGain - gains[cluster[u]];
  }

  int moveWhere(V u) {
    assert(numClusters > 1);

    int to = 0;

    std::vector<W> gains(numClusters, 0);

    for (int j = start[u]; j < start[u + 1]; j++) {
      V v = edges[j];
      int vIndex = cluster[v];

      gains[vIndex] += getEdgeWeight(j);
    }

    int maxGain = std::numeric_limits<int>::min();

    for (size_t i = 0; i < gains.size(); i++) {
      if (gains[i] > maxGain && cluster[u] != i) {
        maxGain = gains[i];
        to = i;
      }
    }
    return to;
  }

  /**
   * Create a priority queue of vertices to be moved.
   */
  BucketQueue<V, W> initializeQueue() {
    BucketQueue<V, W> refinementQueue;

    auto boundaryVector = computeRefinementDegrees();

    for (auto v : boundaryVector) {
      // don't insert taboo vertices
      if (tabooVertices.count(v) > 0) {
        continue;
      }

      auto gain = computeGain(v);

      VLOG(8) << "Inserting " << v << " with gain " << gain;

      refinementQueue.insert(v, gain);
    }

    return refinementQueue;
  }

  W gainTo(V u, int to) {
    W gain;

    for (int j = start[u]; j < start[u + 1]; j++) {
      auto v = edges[j];

      if (cluster[v] == to) {
        gain += getEdgeWeight(j);
      }
    }

    return gain;
  }

  /**
   * Perform a single move.
   */
  bool makeSingleMove(V u, BucketQueue<V, W> &pQueue,
                      std::unordered_set<int> &moved) {
    assert(tabooVertices.count(u) == 0 && "Must not move a taboo vertex.");

    int where = cluster[u];
    int to = moveWhere(u);

    // if the vertex is the only one in its partition don't move
    if (clusterSizes[where] == 1) {
      return false;
    }

    cluster[u] = to;
    moved.insert(u);

    updateAfterMove(u, where, to);

    assert(clusterSizes[where] > 0 && "Cluster can not become empty.");

    // move vertex to this partition

    // update all neighbors, insert new boundary vertices into queue and
    // update gains for those already inside. maybe also delete some
    // non-boundary vertices.
    for (int j = start[u]; j < start[u + 1]; j++) {
      auto v = edges[j];
      auto gain = computeGain(v);

      if (moved.count(v)) continue;

      if (gain < 0) {
        pQueue.removeIfPresent(v);
      } else {
        pQueue.updateOrInsert(v, gain);
      }
    }

    return true;
  }

  /**
   * Compute the normalized cut of the current cluster set.
   */
  double computeNCut() {
    // DegreeVector clusterCuts(numClusters, 0), clusterVolumes(numClusters, 0);

    // for (int u = 0; u < size(); u++) {
    //   clusterVolumes[cluster[u]] += getVertexWeight(u);

    //   for (int j = start[u]; j < start[u + 1]; j++) {
    //     auto v = edges[j];

    //     if (cluster[u] != cluster[v]) {
    //       clusterCuts[cluster[u]] += getEdgeWeight(j);
    //     }
    //   }
    // }

    double val = 0.0;

    for (int i = 0; i < numClusters; i++) {
      val += (double)clusterCuts[i] / (double)clusterVolumes[i];
    }

    return val;
  }

  /**
   * Perform a single round of refinement.
   * Returns the improvement.
   */
  std::pair<double, double> refineOneRound() {
    VLOG(1) << "Beginning refinement.";
    //  Initialize priority queue
    auto refinementQueue = initializeQueue();

    VLOG(1) << "Boundary size: " << refinementQueue.size();

    bool stop = false;
    double target = computeNCut();

    assert(target == target && "Target must not be NaN.");

    double targetBeforeRefine = target;

    std::unordered_set<int> moved;

    VLOG(1) << "Current target: " << target;

    int best_index = -1;
    int i = 0;

    std::vector<V> moves;
    std::vector<int> from;

    while (!stop) {
      if (refinementQueue.empty()) {
        break;
      }

      auto top = refinementQueue.popTop();
      auto topIndex = cluster[top];

      VLOG(8) << "Moving " << top;

      // if the move is taboo, we skip
      if (tabooVertices.count(top) > 0) {
        continue;
      }

      auto success = makeSingleMove(top, refinementQueue, moved);

      // if the move couldn't be made, don't need to do anythig else
      if (!success) {
        VLOG(8) << "Move unsuccessful.";
        continue;
      }

      moves.push_back(top);
      from.push_back(topIndex);

      auto current = computeNCut();

      // VLOG(0) << "Cost after move: " << current;

      if (current < target) {
        best_index = i;
        target = current;
      }

      if (i > best_index + 100) {
        stop = true;
      }

      i++;
    }

    // undo moves
    for (i--; i > best_index; i--) {
      auto v = moves.back();
      auto where = cluster[v];
      auto to = from.back();

      moves.pop_back();
      from.pop_back();

      cluster[v] = to;

      updateAfterMove(v, where, to);
    }

    return {target, targetBeforeRefine - target};
  }

  /**
   * Refine multiple times
   */
  void refine() {
    int maxRounds = 20;  // std::max(maxRounds, )

    double overallImprovement = 0;
    bool changed = false;

    auto target = computeNCut();

    if (target != target) {
      LOG(ERROR) << "Target value is NaN before refinement.";
      exit(-43);
    }

    // printDebugRefine();

    for (int i = 0; i < maxRounds || changed; i++) {
      auto [target, improvement] = refineOneRound();

      // if (target != target) {
      //   LOG(ERROR) << "Target value is NaN";
      //   exit(-42);
      // }

      VLOG(2) << "Round " << i << ". Improvement: " << improvement;

      overallImprovement += improvement;

      if (improvement == 0.0) {
        break;
      }

      if (improvement < 0.001 * target) {
        changed = false;
      } else {
        // TODO: investigate why this causes infinite (?) loop
        changed = true;
      }
    }

    VLOG(1) << "Total improvement: " << overallImprovement;

    for (auto s : clusterSizes) {
      assert(s > 0 &&
             "Clusters can not be empty (or taboo only) after refinement.");
    }
  }

  void printDebugRefine() {
    for (int i = 0; i < numClusters; i++) {
      VLOG(1) << "Edges: " << clusterCuts[i] << " Volume: " << clusterVolumes[i]
              << " Size: " << clusterSizes[i];
    }
  }
};

using Graph = TemplateGraph<int, Edge, false>;
using WeightedGraph = TemplateGraph<int, Edge, true>;

}  // namespace ArrayGraph

#endif  // ARRAY_GRAPH_H_
