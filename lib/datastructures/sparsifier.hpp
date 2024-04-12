#ifndef SPARSIFIER_H_
#define SPARSIFIER_H_

#include <glog/logging.h>
#include <glog/stl_logging.h>

#include <functional>
#include <string>
#include <vector>

#include "array_graph.hpp"
#include "dynamic_graph.hpp"

namespace CutSparsifier {

struct CutResult {
  double cost;
  std::vector<std::pair<int, int>> edges;
  std::vector<std::pair<int, int>> tabooNodes;
};

struct DPCell {
  double cost;
  int cut;
  int volume;
  std::vector<std::pair<int, int>> best;

  DPCell() : cost(std::numeric_limits<double>::infinity()), cut(0), volume(0) {}
  DPCell(double cost, int cut, int volume,
         std::vector<std::pair<int, int>> best)
      : cost(cost), cut(cut), volume(volume), best(std::move(best)) {}

  static std::string bestToString(std::vector<std::pair<int, int>> const &v) {
    if (v.empty()) {
      return "[]";
    }

    std::string res = "[";
    for (auto &[a, b] : v) {
      res += "(" + std::to_string(a) + ", " + std::to_string(b) + "), ";
    }

    res.pop_back();
    res.back() = ']';

    return res;
  }

  std::string to_string() const {
    return "{ cost: " + std::to_string(cost) + ", cut: " + std::to_string(cut) +
           ", volume: " + std::to_string(volume) +
           ", best: " + bestToString(best) + " }";
  }
};

struct DPTable {
  /**
   * The dynamic programming table.
   */
  std::vector<std::vector<DPCell>> row;
  std::vector<std::vector<DPCell>> oldRow;

  int level;
  int numParts;

  /**
   * Initialize the dynamic programming table with the entries for the leaf
   * nodes with no cuts made.
   */
  DPTable(int n, int k, std::vector<int> const &weights)
      : row(n, std::vector<DPCell>(k)),
        oldRow(n, std::vector<DPCell>(k)),
        level(0),
        numParts(k) {
    VLOG(0) << "Initializing DP table.";

    for (size_t i = 0; i < row.size(); i++) {
      row[i][0].cost = 0.0;
      row[i][0].volume = weights[i];

      cutParent(i, 0, weights[i]);
    }

    level++;

    // printRow(n, 2);
  }

  /**
   * After having combined child solutions, check if cutting the parent edge of
   * cell (n, k) gives a better solution for cell (n, k+1) than the ones
   * obtained from combining child solutions.
   *
   * we is the edge weight of the parent edge in the cut sparsifier.
   */
  void cutParent(int n, int k, int we) {
    auto const &cell = row[n][k];

    // if the solution has zero volume, cutting would yield an empty component
    if (cell.volume == 0) return;

    double cutCost = cell.cost + (double)(cell.cut + we) / (cell.volume);
    // TODO check what happens when we add the other term here

    // if cutting is a better strategy, replace the cell
    if (cutCost < row[n][k + 1].cost) {
      auto &targetCell = row[n][k + 1];

      targetCell.cost = cutCost;
      targetCell.cut = we;
      targetCell.volume = 0;

      targetCell.best = cell.best;
      targetCell.best.emplace_back(level, n);
    }
  }

  /**
   * Attempt to cut the parents of each cell, and generate new solutions if that
   * works.
   */
  void cutParents(int n, int we) {
    for (int k = 0; k < numParts - 1; k++) {
      cutParent(n, k, we);
    }
  }

  /**
   * Combine the solution of two cells
   *
   * TODO: we probably don't need to do as much array copying
   */
  // DPCell combineCells(int n1, int k1, int n2, int k2) {
  //   auto const &cell1 = row[n1][k1], &cell2 = row[n2][k2];

  //   std::vector<std::pair<int, int>> best = cell1.best;
  //   best.insert(best.end(), cell2.best.begin(), cell2.best.end());

  //   return DPCell(cell1.cost + cell2.cost, cell1.cut + cell2.cut,
  //                 cell1.volume + cell2.volume, best);
  // }

  /**
   * Combine two DPCell objects into one. If the left object has k cuts and the
   * right has n cuts, the new object will be one where k + n cuts have been
   * made.
   */
  DPCell combineCells(DPCell const &l, DPCell const &r) {
    std::vector<std::pair<int, int>> best = l.best;
    best.insert(best.end(), r.best.begin(), r.best.end());

    return DPCell(l.cost + r.cost, l.cut + r.cut, l.volume + r.volume, best);
  }

  DPCell combineCellsFinal(DPCell const &l, DPCell const &r) {
    std::vector<std::pair<int, int>> best = l.best;
    best.insert(best.end(), r.best.begin(), r.best.end());

    if (l.volume + r.volume == 0) {
      return DPCell(std::numeric_limits<double>::infinity(), 0, 0, {});
    } else {
      return DPCell(l.cost + r.cost +
                        (double)(l.cut + r.cut) / (double)(l.volume + r.volume),
                    l.cut + r.cut, l.volume + r.volume, best);
    }
  }

  std::vector<std::vector<DPCell>> combineNeighbors(
      std::vector<std::vector<DPCell>> const &row) {
    VLOG(6) << "Combining high degree nodes.";
    std::vector<std::vector<DPCell>> res;
    for (size_t i = 0; i + 1 < row.size(); i += 2) {
      std::vector<DPCell> col(numParts);
      for (int k = 0; k < numParts; k++) {
        for (int l = 0; l <= k; l++) {
          auto newCell = combineCells(row[i][l], row[i + 1][k - l]);

          if (newCell.cost < col[k].cost ||
              (newCell.cost == col[k].cost && newCell.volume > col[k].volume)) {
            col[k] = newCell;
          }
        }
      }
      res.push_back(std::move(col));
    }

    // for vectors with uneven amounts of elements need to add the last
    // untouched element
    if (row.size() % 2 != 0) {
      res.push_back(row.back());
    }

    return res;
  }

  std::vector<std::vector<DPCell>> combineNeighborsFinal(
      std::vector<std::vector<DPCell>> const &row) {
    VLOG(6) << "Combining high degree nodes.";
    std::vector<std::vector<DPCell>> res;
    for (size_t i = 0; i + 1 < row.size(); i += 2) {
      std::vector<DPCell> col(numParts);
      for (int k = 0; k < numParts; k++) {
        for (int l = 0; l <= k; l++) {
          auto newCell = combineCellsFinal(row[i][l], row[i + 1][k - l]);

          if (newCell.cost < col[k].cost ||
              (newCell.cost == col[k].cost && newCell.volume > col[k].volume)) {
            col[k] = newCell;
          }
        }
      }
      res.push_back(std::move(col));
    }

    // for vectors with uneven amounts of elements need to add the last
    // untouched element
    if (row.size() % 2 != 0) {
      res.push_back(row.back());
    }

    return res;
  }

  /**
   * Calculate the value of a cell with many children by 'virtually' binarizing
   * the tree.
   */
  void computeNode(int n, std::vector<int> const &children) {
    assert(level == 0 ||
           (children.size() > 0 && "Node must have at least one child."));
    // if there is only one child, we just copy it over, nothing to do anyways
    if (children.size() == 1) {
      VLOG(6) << "Node " << n << " has 1 child, copying on up.";
      auto child = children[0];

      row[n] = oldRow[child];
      return;
    }

    VLOG(6) << "Creating initial child cells.";

    std::vector<std::vector<DPCell>> cells;
    for (auto child : children) {
      cells.push_back(oldRow[child]);
    }

    // printRow(cells, cells.size(), numParts);

    VLOG(6) << "Computed the initial child cells of node " << n
            << ". Array has size " << cells.size();

    while (cells.size() > 1) {
      cells = combineNeighbors(cells);
      VLOG(6) << "Finished iteration. New size " << cells.size();
      // printRow(cells, cells.size(), numParts);
    }

    row[n] = cells[0];
  }

  void computeFinalCell(std::vector<int> const &children) {
    assert(children.size() > 1 &&
           "Top level node should have more than one child.");

    std::vector<std::vector<DPCell>> cells;
    for (auto child : children) {
      cells.push_back(oldRow[child]);
    }

    // printRow(cells, cells.size(), numParts);

    VLOG(6) << "Computed the initial child cells of final node"
            << ". Array has size " << cells.size();

    while (cells.size() > 1) {
      cells = combineNeighborsFinal(cells);
      VLOG(6) << "Finished iteration. New size " << cells.size();
      // printRow(cells, cells.size(), numParts);
    }

    row[0] = cells[0];
  }

  void computeFinal(std::vector<std::vector<int>> const &children) {
    assert(children.size() == 1);

    VLOG(6) << "Computing final cell's solution.";

    std::swap(row, oldRow);

    computeFinalCell(children[0]);

    level++;
  }

  /**
   * Compute a row of the DP Table.
   */
  void computeRow(std::vector<int> const &weights,
                  std::vector<std::vector<int>> const &children) {
    assert(weights.size() == children.size() &&
           "Size of weights and children arrays must match.");

    VLOG(6) << "Computing row " << level;

    std::swap(row, oldRow);

    for (size_t i = 0; i < weights.size(); i++) {
      computeNode(i, children[i]);

      // if we're not in the top node, cut the parent edges
      if (weights.size() > 1) {
        cutParents(i, weights[i]);
      }
    }
    VLOG(6) << "Cut the parents.";

    // printRow(weights.size(), numParts);

    level++;
  }

  DPCell getCell(int n, int k) { return row[n][k]; }

  void printRow(int rowSize, int maxNumParts) {
    for (int i = 0; i < rowSize; i++) {
      for (int k = 0; k < maxNumParts; k++) {
        VLOG(0) << row[i][k].to_string();
      }
    }
  }

  void printRow(std::vector<std::vector<DPCell>> const &row, int rowSize,
                int maxNumParts) {
    for (int i = 0; i < rowSize; i++) {
      for (int k = 0; k < maxNumParts; k++) {
        VLOG(0) << row[i][k].to_string();
      }
    }
  }

  /**
   * This should do it for the dynamic program, but we might have to hedge
   * against parts becoming too small.
   * TODO: Play around with the target functions.
   */
};

template <typename W>
class CutSparsifier {
 private:
  /**
   * We store the tree from the bottom up, with pointers to the parents.
   *
   * The weight vector stores the at the same position the weight of the edge.
   */
  // std::shared_ptr<ArrayGraph::Graph> graph;

  std::vector<std::vector<int>> parents;
  std::vector<std::vector<W>> weights;

  /**
   * Vector holding the children of node (i, j)
   */
  std::vector<std::vector<std::vector<int>>> children;

  std::vector<std::vector<int>> subtreeSizes;
  std::vector<std::vector<int>> subtreeVolumes;

  int totalSize, totalVolume;
  bool done;

  /**
   * Returns a vector such that the vector at index (i, j) contains the
   * children of the tree node (i, j).
   */
  void computeChildren() {
    VLOG(1) << "Parents size: " << parents.size();

    // bottom layer has no children
    for (size_t i = 0; i < parents.size(); i++) {
      std::vector<std::vector<int>> currentChildren(parents[i].size());

      if (i > 0) {
        // iterate over the children and add them to the right parent vectors
        for (size_t j = 0; j < parents[i - 1].size(); j++) {
          currentChildren[parents[i - 1][j]].push_back(j);
        }
      }

      children.push_back(std::move(currentChildren));
    }

    assert(children.size() == parents.size());

    // final layer from the (virtual) root node
    if (parents.back().size() != 1) {
      std::vector<std::vector<int>> last(1);
      for (size_t i = 0; i < parents.back().size(); i++) {
        last[0].push_back(i);
      }
      children.push_back(std::move(last));
    }

    /*
     * This is only for debugging
     */
    for (size_t i = 0; i < children.size() - 1; i++) {
      assert(children[i].size() == parents[i].size());
    }

    for (size_t i = 0; i < parents.size(); i++) {
      for (size_t j = 0; j < parents[i].size(); j++) {
        for (auto c : children[i][j]) {
          assert(static_cast<size_t>(parents[i - 1][c]) == j);
        }
      }
    }

    for (size_t i = 0; i < parents.size() - 1; i++) {
      for (size_t j = 0; j < parents[i].size(); j++) {
        auto const &cs = children[i + 1][parents[i][j]];
        assert(std::find(cs.begin(), cs.end(), j) != cs.end());
      }
    }

    for (size_t i = 1; i < parents.size(); i++) {
      for (size_t j = 0; j < parents[i].size(); j++) {
        int sum = 0;

        for (auto c : children[i][j]) {
          sum += subtreeVolumes[i - 1][c];
        }

        assert(sum == subtreeVolumes[i][j]);
      }
    }
  }

  /**
   * Create a vector to hold labels used by a cut algorithm.
   */
  template <typename L>
  std::vector<std::vector<L>> makeLabelVector(L l) {
    std::vector<std::vector<L>> labels;

    for (size_t i = 0; i < parents.size(); i++) {
      labels.emplace_back(parents[i].size(), l);
    }

    return labels;
  }

  /**
   * Function to compute the cut indicator vector for graph bipartitions.
   */
  std::vector<std::vector<int>> makeBipartitionVector(int i, int j) {
    auto cut = makeLabelVector<int>(0);

    paintChildren(i, j, cut, 1);

    return cut;
  }

  /**
   * Given a vector containing node labels, paint all children of node (i, j)
   * with label l.
   *
   * TODO: make iterative
   */
  void paintChildren(int i, int j, std::vector<std::vector<int>> &labels,
                     int l) {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    labels[i][j] = l;

    for (auto c : children[i][j]) {
      paintChildren(i - 1, c, labels, l);
    }
  }

  /**
   * Given a vector containing node labels, paint all children of node (i, j)
   * with label l that currently contain label k.
   *
   * TODO: make iterative
   */
  void paintChildrenIfMatch(int i, int j, std::vector<std::vector<int>> &labels,
                            int l, int k) {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    if (labels[i][j] != k) return;

    labels[i][j] = l;

    for (auto c : children[i][j]) {
      paintChildrenIfMatch(i - 1, c, labels, l, k);
    }
  }

  template <typename V, typename H>
  void paintChildrenAvoid(int i, int j, std::vector<std::vector<int>> &labels,
                          int l, std::unordered_set<V, H> const &avoid) {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    if (avoid.count(std::make_pair(i, j)) == 1) return;

    labels[i][j] = l;

    for (auto c : children[i][j]) {
      paintChildrenAvoid(i - 1, c, labels, l, avoid);
    }
  }
  std::vector<int> makeCutVectorFromEdges(
      std::vector<std::pair<int, int>> &edges) {
    auto cut = makeLabelVector(0);

    auto cmp = [](std::pair<int, int> &a, std::pair<int, int> &b) {
      return a.first > b.first;
    };

    // sort descendingly by level
    std::sort(edges.begin(), edges.end(), cmp);

    int color = 1;

    // VLOG(0) << "Number of layers: " << parents.size();

    for (auto &p : parents) {
      // VLOG(0) << "Layer size: " << p.size();
    }

    for (auto [i, j] : edges) {
      // VLOG(0) << "Painting (" << i << ", " << j << ")";
      paintChildren(i, j, cut, color);
      color++;
    }

    return cut[0];
  }

  /**
   * Change all parents of index (i, j) using the modifierFunction.
   */
  void modifyParents(int i, int j, std::vector<std::vector<int>> &property,
                     std::function<int(int)> const &modifierFunction) {
    auto parentI = i + 1, parentJ = parents[i][j];

    if (parentJ == -1) return;

    assert(modifierFunction(property[parentI][parentJ]) >= 0);

    property[parentI][parentJ] = modifierFunction(property[parentI][parentJ]);

    modifyParents(parentI, parentJ, property, modifierFunction);
  }

  void modifyParentsIfMatch(int i, int j,
                            std::vector<std::vector<int>> &property,
                            std::function<int(int)> const &modifierFunction,
                            std::vector<std::vector<int>> &labels, int k) {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    auto parentI = i + 1, parentJ = parents[i][j];

    // order matters, otherwise it crashes
    if (parentJ == -1 || labels[parentI][parentJ] != k) return;

    assert(modifierFunction(property[parentI][parentJ]) >= 0);

    property[parentI][parentJ] = modifierFunction(property[parentI][parentJ]);

    modifyParentsIfMatch(parentI, parentJ, property, modifierFunction, labels,
                         k);
  }

  /**
   * Generic minimum search function. It returns the minimum target function
   * and the index of the node where this minimum is achieved.
   */
  std::pair<double, std::pair<int, int>> minimumSearch(
      const std::function<double(int, int)> &targetFunction) {
    double best = std::numeric_limits<double>::infinity();
    int bestI = -1, bestJ = -1;

    for (size_t i = 0; i < parents.size(); i++) {
      // the last layer might give badly formed expressions otherwise (NaN as
      // we get 0 / 0)
      if (parents[i].size() == 1) continue;

      for (size_t j = 0; j < parents[i].size(); j++) {
        double current = targetFunction(i, j);

        if (current < best) {
          best = current;
          bestI = i;
          bestJ = j;
        }
      }
    }

    return {best, {bestI, bestJ}};
  }

  /**
   * Get a set of taboo vertices, i.e. those that have all children cut away.
   */
  std::vector<std::pair<int, int>> tabooVertices(
      std::vector<std::pair<int, int>> &edges) {
    std::vector<std::pair<int, int>> tabooSet;

    // sort them
    std::sort(edges.begin(), edges.end(),
              [](std::pair<int, int> &a, std::pair<int, int> &b) {
                return a.first < b.first ||
                       (a.first == b.first && a.second < b.second);
              });

    size_t i = 0;
    size_t j = 0;

    while (i < edges.size()) {
      std::unordered_map<int, int> vertices;
      while (i < edges.size() && edges[i].first == j) {
        // VLOG(0) << "Edge: " << edges[i];

        auto v = parents[j][edges[i].second];  // parent index is (j + 1, v)

        // VLOG(0) << "Parent: (" << j + 1 << ", " << v << ")";

        if (vertices.find(v) != vertices.end()) {
          vertices[v]++;
        } else {
          vertices[v] = 1;
        }

        i++;
      }

      for (auto [v, numChildren] : vertices) {
        // VLOG(0) << v << " " << numChildren;
        if (numChildren == children[j + 1][v].size()) {
          tabooSet.emplace_back(j + 1, v);
        }
      }

      j++;
    }

    if (tabooSet.size() > 0) {
      // VLOG(0) << "Recursive call.";
      auto rec = tabooVertices(tabooSet);

      tabooSet.insert(tabooSet.end(), rec.begin(), rec.end());
    }

    // VLOG(0) << "Found " << tabooSet.size() << " taboo vertices. They are:";

    return tabooSet;
  }

 public:
  template <typename G>
  CutSparsifier(G *g)
      : totalSize(g->size()), totalVolume(g->getVolume()), done(false) {
    // bottom of subtree sizes tree is all ones
    subtreeSizes.emplace_back(g->size(), 1);

    // compute bottom layer of subtree volumes
    std::vector<int> degrees(g->size());

    for (int i = 0; i < g->size(); i++) {
      degrees[i] = g->getInitialDegree(i);
    }

    subtreeVolumes.push_back(std::move(degrees));
  }

  /**
   * This constructor can be used when reading stored sparsifiers from files.
   */
  CutSparsifier(std::vector<std::vector<int>> &ps,
                std::vector<std::vector<W>> &ws)
      : parents(std::move(ps)), weights(std::move(ws)) {
    // need to calculate subtree sizes, subtree volumes
    totalSize = parents.front().size();
    totalVolume =
        std::reduce(weights.front().begin(), weights.front().end(), 0);

    subtreeSizes.emplace_back(totalSize, 1);
    subtreeVolumes.push_back(parents.front());

    for (int i = 1; i < parents.size(); i++) {
      std::vector<int> sizes(parents[i].size(), 0),
          volumes(parents[i].size(), 0);
      for (int j = 1; j < parents[i - 1].size(); j++) {
        sizes[parents[i - 1][j]] += subtreeSizes.back()[j];
        volumes[parents[i - 1][j]] += subtreeVolumes.back()[j];
      }

      subtreeSizes.push_back(std::move(sizes));
      subtreeVolumes.push_back(std::move(volumes));
    }

    computeChildren();
    done = true;
  }

  /**
   * Using this function we can construct the sparsifier in a bottom up
   * manner.
   */
  void addLayer(int numParents, std::vector<int> &ps, std::vector<W> &ws) {
    assert(ps.size() == subtreeSizes.back().size() &&
           "Number of parent edges must match number of nodes in the graph.");
    assert(ps.size() == ws.size() &&
           "Number of edges and number of weights must match.");

    std::vector<int> sizes(numParents, 0), volumes(numParents, 0);

    for (size_t i = 0; i < ps.size(); i++) {
      sizes[ps[i]] += subtreeSizes.back()[i];
      volumes[ps[i]] += subtreeVolumes.back()[i];
    }

    for (int i = 0; i < numParents; i++) {
      assert(std::reduce(volumes.begin(), volumes.end(), 0) == totalVolume);
    }

    subtreeSizes.push_back(std::move(sizes));
    subtreeVolumes.push_back(std::move(volumes));

    parents.push_back(std::move(ps));
    weights.push_back(std::move(ws));
  }

  /**
   * Finalizes the tree. Sets the parent of each vertex that exists on the
   * current layer to -1 to indicate that it is a top level layer.
   */
  void finalize() {
    auto numComponents = subtreeSizes.back().size();

    // no clue why i am comparing to parents.back.size this really needs a
    // rewrite what a mess :////
    if (parents.back().size() > 1) {
      LOG(INFO) << "Graph is disconnected, adding " << numComponents
                << " cost 0 edges.";
      parents.emplace_back(numComponents, -1);
      weights.emplace_back(numComponents, 0);
    }

    computeChildren();
    done = true;
  }

  CutResult sparseCut() {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    auto [best, index] = minimumSearch([this](int i, int j) {
      return (double)weights[i][j] /
             (double)std::min(subtreeSizes[i][j],
                              totalSize - subtreeSizes[i][j]);
    });

    auto bestI = index.first, bestJ = index.second;

    VLOG(1) << "Best weight: " << weights[bestI][bestJ];
    VLOG(1) << "Best size: " << subtreeSizes[bestI][bestJ];

    return {best, {{bestI, bestJ}}, {}};
  }

  CutResult lowConductanceCut() {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    auto [best, index] = minimumSearch([this](int i, int j) {
      return (double)weights[i][j] /
             (double)std::min(subtreeVolumes[i][j],
                              totalVolume - subtreeVolumes[i][j]);
    });

    auto bestI = index.first, bestJ = index.second;

    VLOG(1) << "Best weight: " << weights[bestI][bestJ];
    VLOG(1) << "Best volume: " << subtreeVolumes[bestI][bestJ];

    return {best, {{bestI, bestJ}}, {}};
  }

  CutResult normalizedCut() {
    assert(done &&
           "Can only use this function if the sparsifier has been finalized.");

    // the greedy function works perfectly fine finding the optimal edge when k
    // = 2, no need to make another here.

    return greedyNormalizedCut(2);
  }

  /**
   * Try and approximate the objective ncut_k(A_1, \dots, A_k) = cut(A_1,
   * A_1^c) / vol(A_1) + \dots + cut(A_k, A_k^c) / vol(A_k) By greedily taking
   * the cheapest edges.
   *
   * Also it's algorithmically just not that great.
   */
  CutResult greedyNormalizedCut(int k) {
    auto cut = makeLabelVector<int>(0);

    // holds the component values
    std::vector<int> cutSizes{0}, componentVolumes{totalVolume};

    // holds the information we need to compute things
    std::vector<std::vector<int>> belowVolumes = subtreeVolumes,
                                  belowCuts = makeLabelVector<int>(0);

    // lookup table to store the cut edges
    struct PairHash {
      std::size_t operator()(const std::pair<int, int> &p) const {
        return std::hash<size_t>{}(p.first) ^ std::hash<size_t>{}(p.second);
      }
    };

    std::unordered_set<std::pair<int, int>, PairHash> cutEdgeIndices;
    std::vector<std::pair<int, int>> cutEdgeVector;

    // this holds the terms that compute to the ncut objective
    // should always contain sum of cutSizes / componentVolumes
    // of components other than the one at the index
    std::vector<double> componentPrecomputation(1, 0.0);

    double target = 0.0;

    for (int n = 0; n < k - 1; n++) {
      assert(std::reduce(componentVolumes.begin(), componentVolumes.end(), 0) ==
                 totalVolume &&
             "Component volumes must form a partition.");

      auto [val, index] = minimumSearch(
          [&cutEdgeIndices, &cut, &cutSizes, &componentVolumes, &belowVolumes,
           &belowCuts, &componentPrecomputation, this](int i, int j) {
            if (cutEdgeIndices.count(std::make_pair(i, j)) == 1) {
              return std::numeric_limits<double>::infinity();
            }

            size_t componentIndex = cut[i][j];

            return (double)(weights[i][j] + belowCuts[i][j]) /
                       (double)belowVolumes[i][j] +
                   (double)(weights[i][j] + cutSizes[componentIndex] -
                            belowCuts[i][j]) /
                       (double)(componentVolumes[componentIndex] -
                                belowVolumes[i][j]) +
                   componentPrecomputation[componentIndex];
          });

      auto oldIndex = cut[index.first][index.second];

      // values of the old component
      auto oldVolume = componentVolumes[oldIndex];
      auto oldCut = cutSizes[oldIndex];
      auto oldValue = (double)oldCut / (double)oldVolume;

      auto newVolume = belowVolumes[index.first][index.second];
      auto newWeight = weights[index.first][index.second];
      auto newBelowCut = belowCuts[index.first][index.second];

      VLOG(1) << "new Volume: " << newVolume;

      // create new component
      paintChildrenIfMatch(index.first, index.second, cut, n + 1, oldIndex);
      // paintChildrenAvoid(index.first, index.second, cut, n + 1,
      // cutEdgeIndices);

      componentVolumes[oldIndex] -= newVolume;
      componentVolumes.push_back(newVolume);

      assert(std::all_of(componentVolumes.begin(), componentVolumes.end(),
                         [](int a) { return a >= 0; }) &&
             "All volumes must be non-negative.");

      cutSizes[oldIndex] += newWeight - newBelowCut;
      cutSizes.push_back(newBelowCut + newWeight);

      auto x = (double)cutSizes[oldIndex] / (double)componentVolumes[oldIndex];
      auto y = (double)cutSizes.back() / (double)componentVolumes.back();

      // compute component values
      double newValue = 0.0;
      for (int i = 0; i <= n; i++) {
        if (i != oldIndex) {
          componentPrecomputation[i] -= oldValue;
          componentPrecomputation[i] += x + y;
        } else {
          componentPrecomputation[i] += y;
        }
        newValue += (double)cutSizes[i] / (double)componentVolumes[i];
      }
      componentPrecomputation.push_back(newValue);

      assert(std::all_of(cutSizes.begin(), cutSizes.end(),
                         [](int a) { return a >= 0; }) &&
             "All cut values must be non-negative.");

      // change parents
      modifyParentsIfMatch(
          index.first, index.second, belowVolumes,
          [newVolume](int val) { return val - newVolume; }, cut, oldIndex);

      modifyParentsIfMatch(
          index.first, index.second, belowCuts,
          [newWeight, newBelowCut](int val) {
            return val + newWeight - newBelowCut;
          },
          cut, oldIndex);

      // assert(target <= val && "Target value must be monotonic in k.");
      target = val;

      cutEdgeIndices.insert(index);
      cutEdgeVector.push_back(index);
    }

    for (int i = 0; i < k; i++) {
      VLOG(1) << "Edges: " << cutSizes[i] << " Volume: " << componentVolumes[i];
    }

    return {target, cutEdgeVector, tabooVertices(cutEdgeVector)};
  }

  /**
   * Solving Normalized cut using a bottom up dynamic program.
   */
  // std::pair<double, std::vector<int>>
  CutResult normalizedCutDynamic(int k) {
    VLOG(2) << "Creating a DP table of size " << k;

    DPTable table(totalSize, k, subtreeVolumes[0]);

    VLOG(0) << "Computing the DP table";

    // compute all cells except for the final one
    for (size_t i = 1; i < children.size() - 1; i++) {
      VLOG(2) << "Level " << i << " of " << children.size() - 1
              << "(Size: " << weights.size() << ")";
      table.computeRow(weights[i], children[i]);
    }

    // compute the final row
    // table.computeRow(weights[children.size() - 1],
    //                  children[children.size() - 1]);
    table.computeFinal(children.back());

    VLOG(1) << "Computed a normalized cut, k = " << k;

    auto cell = table.getCell(0, k - 1);

    VLOG(1) << "Result: ";
    VLOG(1) << cell.to_string();

    return {cell.cost, cell.best, tabooVertices(cell.best)};
  }

  /**
   * Compute the number of components in the graph.
   */
  int numComponents() { return parents.back().size(); }

  /**
   * Check if the graph which we sparsify is globally connected.
   */
  bool connected() { return numComponents() == 1; }

  /**
   * Answer pairwise connectivity queries for vertices u and v, by going up in
   * the hierarchy and checking whether they have a common parent.
   */
  bool connected(int u, int v) {
    // exit early if the whole graph is connected.
    if (connected()) return true;

    for (int i = 0; i < parents.size(); i++) {
      u = parents[u];
      v = parents[v];

      if (u == v) {
        return true;
      }
    }
    return false;
  }

  /**
   * Construct a weighted graph from the sparsifier.
   *
   * TODO: Figure out which vertex weights to use?
   */
  DynamicGraph::WeightedGraph *toGraph() {
    assert(done && "Can only construct a graph from a finished sparsifier.");

    int n = 0, layerSize = 0;
    std::vector<DynamicGraph::Edge> es;
    std::vector<int> ws;

    for (size_t i = 0; i < parents.size() - 1; i++) {
      assert(n == layerSize);

      layerSize += parents[i].size();

      for (size_t j = 0; j < parents[i].size(); j++) {
        es.emplace_back(n, parents[i][j] + layerSize);
        ws.emplace_back(weights[i][j]);

        n++;
      }
    }

    // to count the top vertex;
    n++;

    return new DynamicGraph::WeightedGraph(n, es, ws, std::vector<int>(n, 1));
  }

  void printDebug() {
    for (size_t i = 0; i < parents.size(); i++) {
      VLOG(0) << parents[i];
      VLOG(0) << weights[i];
    }
  }

  std::pair<int, W> getParentEdge(int i, int j) const {
    return std::make_pair(parents[i][j], weights[i][j]);
  }

  int getLayerSize(int i) const { return parents[i].size(); }

  int size() const { return parents.size(); }
};
}  // namespace CutSparsifier

#endif  // SPARSIFIER_H_
