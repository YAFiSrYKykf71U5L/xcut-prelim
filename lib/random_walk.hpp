#ifndef RANDOM_WALK_H_
#define RANDOM_WALK_H_

#include <glog/logging.h>
#include <glog/stl_logging.h>

#include <iterator>
#include <limits>
#include <random>
#include <unordered_set>

#include "datastructures/array_graph.hpp"
#include "datastructures/graph_hierarchy.hpp"
#include "datastructures/sparsifier.hpp"

namespace RandomWalk {

struct Bound {
  using iterType = std::vector<int>::iterator;

  iterType begin;
  iterType end;

  Bound(iterType b, iterType e) : begin(b), end(e) {}
};

template <typename G, typename W>
class RandomWalk {
 private:
  G *graph;

  std::vector<W> nodeDistribution;
  std::vector<W> oldNodeDistribution;

  /**
   * p describes how likely the random walk is to stay at the current node.
   * q is 1.0 - p.
   */
  W p, q;

  /**
   * Pointer to the randomness source. Storing multiple random number generators
   * is costly. Ideally not more than 1 per thread.
   */
  std::mt19937 *randomGen;

  /**
   * A standard normal distribution.
   */
  std::normal_distribution<W> gaussian;

  /**
   * Sample new values for the nodes and subtract the average, so that a
   * perfectly converged random walk will have 0 on all nodes.
   */
  void resetDistribution() {
    W avg = 0.0;

    for (size_t i = 0; i < nodeDistribution.size(); i++) {
      nodeDistribution[i] = gaussian(*randomGen);
      avg += nodeDistribution[i];
    }

    avg = avg / nodeDistribution.size();

    for (size_t i = 0; i < nodeDistribution.size(); i++)
      nodeDistribution[i] -= avg;
  }

  template <typename iterType>
  void resetDistribution(iterType begin, iterType end) {
    W avg = 0.0;

    for (auto it = begin; it != end; it++) {
      nodeDistribution[*it] = gaussian(*randomGen);
      avg += nodeDistribution[*it];
    }

    avg = avg / std::distance(begin, end);

    for (auto it = begin; it != end; it++) nodeDistribution[*it] -= avg;
  }

  /**
   * Compute the value of a node based on the values stored in
   * oldNodeDistribution. This can be used when computing an iteration of the
   * random walk to update the new nodeDistribution.
   */
  W computeDistributionValue(int u) const {
    W val = 0.0;

    /* Using template deduction we save a multiplication here for unweighted
     * graph types. if constexpr is evaluated at compile time and prunes
     * unused branches.
     */
    if (graph->getDegree(u) == 0) return oldNodeDistribution[u];

    if constexpr (!G::Weighted) {
      for (auto e = graph->cbeginEdges(u); e != graph->cendActiveEdges(u);
           e++) {
        val += oldNodeDistribution[*e];
      }
    } else {
      auto e = graph->cbeginEdges(u);
      auto w = graph->cbeginWeights(u);
      for (; e != graph->cendActiveEdges(u); e++, w++) {
        val += *w * oldNodeDistribution[*e];
      }
    }

    return p * oldNodeDistribution[u] + q * val / graph->getDegree(u);
  }

 public:
  RandomWalk(G *g, std::mt19937 *rand)
      : graph(g),
        nodeDistribution(g->size()),
        oldNodeDistribution(g->size()),
        randomGen(rand),
        gaussian() {
    p = 0.1;
    q = 1.0 - p;

    resetDistribution();
  }

  void setLaziness(W value) {
    if (p < 0.0 || p > 1.0) {
      VLOG(1) << "p needs to be a probability, i.e. in the range [0, 1].";
      return;
    }

    p = value;
    q = 1.0 - value;
  }

  void iterate() {
    std::swap(nodeDistribution, oldNodeDistribution);

    for (int u = 0; u < graph->size(); u++) {
      nodeDistribution[u] = computeDistributionValue(u);
    }
  }

  /**
   * Update the distribution in nodeDistribution for the given range of values.
   */
  template <typename iterType>
  void iterate(iterType begin, iterType end) {
    for (iterType it = begin; it != end; it++) {
      int u = *it;

      nodeDistribution[u] = computeDistributionValue(u);
    }
  }

  /**
   * This is used to iterate on only "active" components, i.e. once we no longer
   * need to iterate on a component, for example because it is an expander, we
   * might want to skip this component.
   */
  void iterate(std::vector<Bound> const &bounds) {
    std::swap(nodeDistribution, oldNodeDistribution);

    for (auto &bound : bounds) iterate(bound.begin, bound.end);
  }

  int iterate(std::vector<Bound> const &bounds, int n) {
    if (n < 0) {
      VLOG(1) << "The number of iterations must be positive, got " << n << ".";
      return n;
    }

    while (n--) {
      iterate(bounds);
    }

    return n;
  }

  /**
   * Perform n iterations of the random walk.
   */
  int iterate(int n) {
    if (n < 0) {
      LOG(ERROR) << "The number of iterations must be positive, got " << n
                 << ".";
      return n;
    }

    while (n--) {
      iterate();
    }

    return n;
  }

  int iterateDiffusion(int n) {
    if (n < 0) {
      LOG(ERROR) << "The number of iterations must be positive, got " << n
                 << ".";
      return n;
    }

    while (n--) {
      iterate();
    }

    return n;
  }

  /**
   * Used to access the current node distribution.
   */
  W &operator[](const int u) { return nodeDistribution[u]; }
  const W &operator[](const int u) const { return nodeDistribution[u]; }

  /**
   * Assigns a new graph and resets the distribution.
   *
   * This lets us reuse the already assigned memory for the node distributions
   * (hopefully).
   */
  int assignGraph(G *g) {
    graph = g;

    nodeDistribution.resize(graph->size());
    oldNodeDistribution.resize(graph->size());

    resetDistribution();

    return graph->size();
  }

  /**
   * Ensure that the average of the distribution is 0 on the given iterator
   * range.
   */
  template <typename iterType>
  W recenterDistribution(iterType begin, iterType end) {
    assert(std::distance(begin, end) > 0 && "Distance must be positive.");

    // resetDistribution(begin, end);
    W avg = 0.0;

    for (auto it = begin; it != end; it++) {
      avg += nodeDistribution[*it];
    }

    avg /= static_cast<W>(std::distance(begin, end));

    for (auto it = begin; it != end; it++) {
      nodeDistribution[*it] -= avg;
      // nodeDistribution[*it] *= 2;
    }

    return avg;
  }

  /**
   * Potential of the random walk, measured in the \ell_\infty norm.
   */
  template <typename iterType>
  W potential(iterType begin, iterType end) const {
    W minimum = std::numeric_limits<W>::max(),
      maximum = std::numeric_limits<W>::lowest();

    for (auto it = begin; it != end; it++) {
      minimum = std::min(minimum, nodeDistribution[*it]);
      maximum = std::max(maximum, nodeDistribution[*it]);
    }

    return maximum - minimum;
  }

  /**
   * Potential of the random walk, measured in the squared \ell_2 norm.
   */
  template <typename iterType>
  W L2Potential(iterType begin, iterType end) const {
    W potential = 0.0;
    W avg = 0.0;

    for (auto it = begin; it != end; it++) {
      avg += nodeDistribution[*it];
    }

    avg /= static_cast<W>(std::distance(begin, end));

    for (auto it = begin; it != end; it++) {
      auto u = *it;

      auto divergence = nodeDistribution[u] - avg;

      potential += divergence * divergence;
    }

    return potential;
  }

  /**
   * Potential of the random walk, measured in the squared \ell_2 norm.
   */
  template <typename iterType>
  W L2PotentialWeighted(iterType begin, iterType end) const {
    W potential = 0.0;
    W avg = 0.0;

    for (auto it = begin; it != end; it++) {
      avg += nodeDistribution[*it];
    }

    avg /= static_cast<W>(std::distance(begin, end));

    for (auto it = begin; it != end; it++) {
      auto u = *it;

      auto divergence = nodeDistribution[u] - avg;

      potential += graph->getInitialDegree(u) * divergence * divergence;
    }

    return potential;
  }

  void reset() { resetDistribution(); }
};

using wType = double;

enum CutType {
  Cut,
  Expander,
  Inconclusive,
};

template <typename G, typename W>
class SparseCut {
 private:
  G *graph;

  std::vector<int> perm;

 public:
  SparseCut(G *g) : graph(g), perm(g->size()) {
    std::iota(perm.begin(), perm.end(), 0);
  }
};

template <typename G, typename W>
class ExpanderDecomposition {
 private:
  /**
   * A pointer to the graph which we are decomposing.
   */
  G *graph;

  /**
   * A vector storing the sorted vertices. Each component is sorted by the
   * values of the random walk and delimited by a bound from the componentBounds
   * or expanderBounds vectors.
   */
  std::vector<int> perm;

  /**
   * A vector used in the levelCut routine
   */
  int cutCalls;
  std::vector<int> cut;

  /**
   * This stores for each component where it ends in the graph.
   */
  std::vector<Bound> componentBounds;
  std::vector<Bound> expanderBounds;

  RandomWalk<G, W> randomWalk;

  /**
   * Parameters of the random walk solver. Maybe automatic tuning would be
   * useful?
   *
   * itersBetweenCheck provides a tradeoff between quality and running time,
   * lower values make computing the expander hierarchy more time-consuming but
   * also lead to better quality results.
   */
  static constexpr int itersBetweenCheck = 1;
  static constexpr int maxIters = 10;

  W potentialThreshold;
  W ratio;
  int boost;

  /**
   * Compute the best level cut in a range of vertices.
   *
   * Returns the sparsity of the cut and the index of the cut's position in the
   * range.
   */
  template <typename iterType>
  std::pair<wType, int> levelCut(iterType begin, iterType end) {
    // const int boundaryLinkedness = 100;

    const int componentSize = std::distance(begin, end);

    assert(componentSize > 1 &&
           "Component must contain at least two vertices.");

    int componentWeight = graph->getComponentVolume(*begin);

    int cutSize = 0, bestIndex = -1, bestCutSize = 0;  // bestCutWeight = 0;
    W bestSparsity = std::numeric_limits<W>::infinity();

    // sort the elements in the indicated range as necessary
    std::sort(begin, end,
              [this](int u, int v) { return randomWalk[u] < randomWalk[v]; });

    // sorting does not invalidate iterators
    int i = 1;
    int cutWeight = 0;
    for (iterType it = begin; it != end - 1; it++) {
      auto u = *it;
      cut[u] = cutCalls;

      if (boost == 1) {
        cutWeight += graph->getInitialDegree(u);
      } else {
        // we want to weight cut edges more
        cutWeight += graph->getDegree(u) + boost * graph->getNumCutEdges(u);
      }

      if constexpr (!G::Weighted) {
        for (auto e = graph->cbeginEdges(u); e != graph->cendActiveEdges(u);
             e++) {
          if (cut[*e] == cutCalls)
            cutSize--;
          else
            cutSize++;
        }
      } else {
        auto e = graph->cbeginEdges(u);
        auto w = graph->cbeginWeights(u);

        for (; e != graph->cendActiveEdges(u); e++, w++) {
          if (cut[*e] == cutCalls)
            cutSize -= *w;
          else
            cutSize += *w;
        }
      }

      assert(cutSize > 0 &&
             "Cutting a connected component must disconnect edges.");

      W sparsity =
          (W)cutSize / (W)std::min(cutWeight, componentWeight - cutWeight);

      assert(sparsity >= 0.0 && "Sparsity can never be negative");

      if (sparsity < bestSparsity) {
        bestSparsity = sparsity;
        bestCutSize = cutSize;
        bestIndex = i;
      }
      i++;
    }
    assert(bestSparsity > 0.0 && "Sparsity can not be 0 in connected graph.");

    VLOG(4) << "Best sparsity is " << bestSparsity << " found in iteration "
            << bestIndex;
    VLOG(4) << "The cut disconnects " << bestCutSize << " edges.";

    cutCalls++;

    return std::make_pair(bestSparsity, bestIndex);
  }

  /**
   * Compute the level cut of a set of vertices, but pick the most balanced cut
   * that fulfills the threshold criteria given.
   *
   * Returns the sparsity of the cut and the index of the cut's position in the
   * range.
   */
  template <typename iterType>
  std::pair<wType, int> balancedThresholdCut(iterType begin, iterType end,
                                             double phi) {
    // const int boundaryLinkedness = 100;

    const int componentSize = std::distance(begin, end);
    const int middle = componentSize / 2;
    const int quarter = componentSize / 4;
    const int threequarter = middle + quarter;

    assert(componentSize > 1 &&
           "Component must contain at least two vertices.");

    int componentWeight = graph->getComponentVolume(*begin);

    int cutSize = 0, bestIndex = -1, bestCutSize = 0;  // bestCutWeight = 0;
    W bestSparsity = std::numeric_limits<W>::infinity();

    int bestBalancedIndex = -1, bestBalancedCutSize = 0;
    W bestBalancedSparsity = std::numeric_limits<W>::infinity();

    // sort the elements in the indicated range as necessary
    std::sort(begin, end,
              [this](int u, int v) { return randomWalk[u] < randomWalk[v]; });

    // sorting does not invalidate iterators
    int i = 1;
    int cutWeight = 0;
    for (iterType it = begin; it != end - 1; it++) {
      auto u = *it;
      cut[u] = cutCalls;

      if (boost == 1) {
        cutWeight += graph->getInitialDegree(u);
      } else {
        // we want to weight cut edges more
        cutWeight += graph->getDegree(u) + boost * graph->getNumCutEdges(u);
      }

      if constexpr (!G::Weighted) {
        for (auto e = graph->cbeginEdges(u); e != graph->cendActiveEdges(u);
             e++) {
          if (cut[*e] == cutCalls)
            cutSize--;
          else
            cutSize++;
        }
      } else {
        auto e = graph->cbeginEdges(u);
        auto w = graph->cbeginWeights(u);
        for (; e != graph->cendActiveEdges(u); e++, w++) {
          if (cut[*e] == cutCalls)
            cutSize -= *w;
          else
            cutSize += *w;
        }
      }

      assert(cutSize > 0 &&
             "Cutting a connected component must disconnect edges.");

      W sparsity =
          (W)cutSize / (W)std::min(cutWeight, componentWeight - cutWeight);

      assert(sparsity >= 0.0 && "Sparsity can never be negative");

      if (sparsity < bestBalancedSparsity && i >= quarter &&
          i <= threequarter) {
        bestBalancedSparsity = sparsity;
        bestBalancedCutSize = cutSize;
        bestBalancedIndex = i;
      }

      if (sparsity < bestSparsity) {
        bestSparsity = sparsity;
        bestCutSize = cutSize;
        bestIndex = i;
      }

      i++;
    }
    if (bestBalancedSparsity < phi && bestIndex < quarter ||
        bestIndex > threequarter) {
      bestSparsity = bestBalancedSparsity;
      bestIndex = bestBalancedIndex;
      bestCutSize = bestBalancedCutSize;
    }

    assert(bestSparsity > 0.0 && "Sparsity can not be 0 in connected graph.");

    VLOG(4) << "Best sparsity is " << bestSparsity << " found in iteration "
            << bestIndex;
    VLOG(4) << "The cut disconnects " << bestCutSize << " edges.";

    cutCalls++;

    return std::make_pair(bestSparsity, bestIndex);
  }

 public:
  ExpanderDecomposition(G *g, std::mt19937 *rand, W potential, int boost)
      : graph(g),
        perm(g->size()),
        cutCalls(1),
        cut(g->size()),
        componentBounds(),
        expanderBounds(),
        randomWalk(g, rand),
        potentialThreshold(potential),
        boost(boost) {
    // if graph has more than one cut component we need to reset it
    if (graph->getNumComponents() > 1) {
      VLOG(4) << "Resetting the graph's cut structure.";
      graph->resetCutArray();
    }

    std::iota(perm.begin(), perm.end(), 0);

    // create a bound for every connected component
    auto dfsResult = graph->dfs(perm.begin(), perm.end());

    VLOG(1) << "Graph has " << dfsResult.numComponents
            << " connected components.";

    if (dfsResult.numComponents > 1) {
      auto comp = [&dfsResult](int u, int v) {
        return dfsResult.labels[u] < dfsResult.labels[v];
      };

      std::sort(perm.begin(), perm.end(), comp);

      // make new bounds
      auto const &sizes = dfsResult.sizes;

      int runningSum = 0;
      for (int i = 0; i < dfsResult.numComponents; i++) {
        if (i > 0) {
          graph->makeNewComponentFromConnected(
              perm.begin() + runningSum, perm.begin() + runningSum + sizes[i]);
        }

        componentBounds.emplace_back(perm.begin() + runningSum,
                                     perm.begin() + runningSum + sizes[i]);

        runningSum += sizes[i];
      }

    } else {
      VLOG(1) << "Initializing the single bound.";
      componentBounds.emplace_back(perm.begin(), perm.end());
    }

    assert(componentBounds.size() ==
               static_cast<size_t>(graph->getNumComponents()) &&
           "Each connected component is separate cut component.");

    assert(componentBounds.size() ==
               static_cast<size_t>(dfsResult.numComponents) &&
           expanderBounds.size() == 0 && "Bounds sizes are incorrect.");
  }

  /**
   * Compute an expander decomposition using the approach of Saranurak and
   * Wang, i.e. find a sparse cut, remove it and recurse.
   *
   * Instead of the cut-matching game we use spectral methods.
   */
  double computeExpanderDecomp(W targetSparsity) {
    bool done = false;
    int numFinal = 0, numIters = 0;

    VLOG(1) << "Computing an expander decomposition of the graph with target "
               "sparsity "
            << targetSparsity;

    while (!done) {
      // iterate the random walk

      // active only iteration
      randomWalk.iterate(componentBounds, itersBetweenCheck);

      VLOG(1) << "Iteration: " << ++numIters;

      VLOG(1) << "Total number of active components: "
              << componentBounds.size();
      VLOG(1) << "Finalized vertices: " << numFinal;

      // compute cuts on each component
      size_t i = 0;
      while (i < componentBounds.size()) {
        auto &bound = componentBounds[i];

        auto componentSize = std::distance(bound.begin, bound.end);

        assert(componentSize > 0 && "Bound must not be empty.");

        VLOG(2) << "Processing component " << i << " of size " << componentSize;

        // deal with size 1 components
        if (componentSize == 1) {
          VLOG(2) << "Component is a singleton.";
          // move component to expander array
          expanderBounds.push_back(std::move(bound));

          std::swap(componentBounds[i], componentBounds.back());
          componentBounds.pop_back();

          numFinal++;

          continue;
        }

        // const auto [sparsity, index] =
        //     balancedThresholdCut(bound.begin, bound.end, targetSparsity);
        const auto [sparsity, index] = levelCut(bound.begin, bound.end);
        VLOG(2) << "Found cut of sparsity " << sparsity << ".";

        if (componentSize == 2 && sparsity >= targetSparsity) {
          VLOG(2) << "Size 2 component does not contain sparse cut.";

          expanderBounds.push_back(std::move(bound));

          std::swap(componentBounds[i], componentBounds.back());
          componentBounds.pop_back();

          numFinal += 2;

          continue;
        }

        if (sparsity < targetSparsity) {
          VLOG(2) << "Splitting the graph.";

          auto begin = bound.begin;
          auto mid = bound.begin + index;
          auto end = bound.end;

          VLOG(2) << "New components are of size "
                  << std::distance(bound.begin, mid) << " and "
                  << std::distance(mid, bound.end);

          // make the cut in the graph
          auto [firstSizes, secondSizes] =
              graph->makeNewComponentsFromCut(bound.begin, mid, end);

          // sorting function for the permutation
          auto comp = [this](int u, int v) {
            return graph->getCutIndex(u) < graph->getCutIndex(v);
          };

          // and update the bounds
          if (firstSizes.size() > 1) {
            std::sort(begin, mid, comp);

            bound.end = begin + firstSizes[0];

            int runningSum = firstSizes[0];

            for (int i = 1; i < firstSizes.size(); i++) {
              componentBounds.emplace_back(begin + runningSum,
                                           begin + runningSum + firstSizes[i]);
              runningSum += firstSizes[i];
            }

            assert(begin + runningSum == mid);
          } else {
            bound.end = mid;
          }

          if (secondSizes.size() > 1) {
            std::sort(mid, end, comp);

            int runningSum = 0;
            for (int i = 0; i < secondSizes.size(); i++) {
              componentBounds.emplace_back(mid + runningSum,
                                           mid + runningSum + secondSizes[i]);

              runningSum += secondSizes[i];
            }

            assert(mid + runningSum == end);
          } else {
            componentBounds.emplace_back(mid, end);
          }
        } else {
          // if no sparse cut is found, check the potential

          auto potential =
              randomWalk.L2PotentialWeighted(bound.begin, bound.end);
          // auto potential = randomWalk.L2Potential(bound.begin, bound.end);
          // auto potential = randomWalk.potential(bound.begin, bound.end);
          // compute cut

          // check result
          // if expander
          VLOG(2) << "Potential of component " << i << ": " << potential;

          assert(potential >= 0.0 && "Potential must be positive.");

          if (potential < potentialThreshold) {
            VLOG(2) << "Component is expanding. Finalizing it.";
            numFinal += std::distance(bound.begin, bound.end);

            expanderBounds.push_back(std::move(bound));

            std::swap(componentBounds[i], componentBounds.back());
            componentBounds.pop_back();

            continue;
          }
        }

        i++;
      }

      if (componentBounds.empty()) {
        done = true;
      }

      assert(componentBounds.size() + expanderBounds.size() ==
                 static_cast<size_t>(graph->getNumComponents()) &&
             "Each bound must have a component");
    }

    assert(expanderBounds.size() ==
               static_cast<size_t>(graph->getNumComponents()) &&
           "The number of bounds and components must match.");

    VLOG(0) << "Finished computing an expander decomposition. There are "
            << graph->getNumComponents() << " components and "
            << graph->getTotalNumCutEdges() << " cut edges.";

    ratio =
        2 * (double)graph->getTotalNumCutEdges() / (double)graph->getVolume();
    double nodeRatio = (double)graph->getNumComponents() / graph->size();

    VLOG(0) << "The ratio of cut edges is " << ratio << ".";
    VLOG(0) << "The ratio of node contraction is " << nodeRatio << ".";

    return nodeRatio;
  }

  void setPotential(W potential) { potentialThreshold = potential; }
};

template <typename W, typename G, typename GW>
class ExpanderHierarchy {
 private:
  G *graph;
  GW *weightedGraph;
  GraphHierarchy<G, GW> graphStack;

  std::mt19937 *randomSource;

  CutSparsifier::CutSparsifier<int> sparsifier;

  W potentialThreshold;
  int boost;

 public:
  ExpanderHierarchy(G *g, std::mt19937 *rand, W potential, int boost)
      : graph(g),
        graphStack(g),
        randomSource(rand),
        sparsifier(g),
        potentialThreshold(potential),
        boost(boost) {}

  ExpanderHierarchy(G *g, std::mt19937 *rand)
      : ExpanderHierarchy(g, rand, 0.000001, 1.0) {}

  ~ExpanderHierarchy() { graphStack.del(); }

  /**
   * Compute an expander decomposition with an exponential backoff
   */
  template <typename H>
  W tunedExpanderDecomp(H *g, W &targetSparsity, double backoff,
                        double ratioTarget) {
    W ratio;
    int iters = 0;

    VLOG(0) << "Start autotuned expander decomp.";

    do {
      auto decomposition = ExpanderDecomposition<H, double>(
          g, randomSource, potentialThreshold, boost);

      // if too many disconnected components, this can be an issue here

      ratio = decomposition.computeExpanderDecomp(targetSparsity);

      if (ratio > ratioTarget) {
        targetSparsity *= backoff;
        g->resetCutArray();

        VLOG(0) << "Did not hit ratio target. Backing off to "
                << targetSparsity;
      }
      if (iters++ > 20) {
        LOG(ERROR) << "Too many iterations to determine a good sparsity.";
        exit(1);
      }
    } while (ratio > ratioTarget);

    VLOG(0) << "Finish autotuned expander decomp in " << iters
            << " iterations.";

    return ratio;
  }

  /**
   * Compute an expander hierarchy over the whole graph using a given target
   * sparsity.
   */
  void computeExpanderHierarchy(W targetSparsity) {
    if (targetSparsity == 0) {
      VLOG(0) << "Computing an expander hierarchy, automatically tuning the "
                 "sparsity.";
    } else {
      VLOG(0) << "Computing an expander hierarchy on the graph with target "
                 "sparsity "
              << targetSparsity;
    }

    VLOG(0) << "Computing the decomposition of the initial graph..";

    if (graph->getVolume() == 0) {
      VLOG(0) << "Graph is empty";
      return;
    }

    double reductionRatio;
    W target = 0.3;

    if (targetSparsity == 0.0) {
      reductionRatio = tunedExpanderDecomp(graph, target, 0.8, 0.95);
    } else {
      auto decomposition = ExpanderDecomposition<G, double>(
          graph, randomSource, potentialThreshold, boost);
      reductionRatio = decomposition.computeExpanderDecomp(targetSparsity);
    }

    VLOG(0) << "Finished initial deomposition.";

    if (reductionRatio > 0.99) {
      LOG(ERROR) << "The hierarchy will (probably) not terminate for this "
                    "parameter choice.";
      exit(42);
    }

    // get which component the vertices belong to
    std::vector<int> parents = graph->getCutArray();

    // weights of edges in sparsifier are the degrees of vertices
    std::vector<int> weights = graph->getInitialDegreeArray();

    sparsifier.addLayer(graph->getNumComponents(), parents, weights);

    if (graph->getNumComponents() >= 0.995 * graph->size()) {
      VLOG(0) << "The hierarchy will (probably) not terminate for this "
                 "parameter choice.";
      exit(42);
    }

    weightedGraph = graph->contractedGraphArray();

    while (true) {
      if (weightedGraph->size() == 1) {
        VLOG(0) << "The graph is a singleton. The hierarchy has been computed.";
        graphStack.push(weightedGraph);
        break;
      }
      if (weightedGraph->getVolume() == 0) {
        VLOG(0) << "The graph is empty. We are done.";
        graphStack.push(weightedGraph);
        break;
      }

      if (targetSparsity == 0.0) {
        reductionRatio = tunedExpanderDecomp(weightedGraph, target, 0.8, 0.95);
      } else {
        auto decomposition = ExpanderDecomposition<GW, double>(
            weightedGraph, randomSource, potentialThreshold, boost);
        reductionRatio = decomposition.computeExpanderDecomp(targetSparsity);
      }

      VLOG(0) << "Finished computing the decomposition, contracting "
                 "the graph...";

      // update the sparsifier
      std::vector<int> parents = weightedGraph->getCutArray();
      std::vector<int> weights = weightedGraph->getInitialDegreeArray();
      sparsifier.addLayer(weightedGraph->getNumComponents(), parents, weights);

      // if (weightedGraph->getNumComponents() == weightedGraph->size()) {
      //   VLOG(0)
      //       << "The hierarchy will not terminate for this parameter choice.";
      //   exit(42);
      // }

      if (reductionRatio > 0.99) {
        LOG(ERROR) << "The hierarchy will (probably) not terminate for this "
                      "parameter choice.";
        exit(42);
      }

      // compute the contracted graph
      graphStack.push(weightedGraph);

      weightedGraph = weightedGraph->contractedGraphArray();

      VLOG(0) << "Contracted graph has " << weightedGraph->size()
              << " vertices and " << weightedGraph->getNumEdges() << " edges.";
    }
    VLOG(0) << "Finalizing the sparsifier.";
    sparsifier.finalize();
  }

  CutSparsifier::CutSparsifier<int> const &getSparsifier() const {
    return sparsifier;
  }

  GraphHierarchy<G, GW> getHierarchy() const { return graphStack; }

  void setPotential(W potential) { potentialThreshold = potential; }
};
}  // namespace RandomWalk

#endif  // RANDOM_WALK_H_
