#ifndef CUT_EVALUATORS_H_
#define CUT_EVALUATORS_H_

#include <memory>
#include <vector>

#include "lib/datastructures/array_graph.hpp"

/**
 * Utility file for computing the objective value of certain cut functions.
 */

/**
 * Compute the size of a cut going between two components.
 */
int cutSize(ArrayGraph::Graph *graph, std::vector<int> const &cut);

double sparseCut(ArrayGraph::Graph *graph, std::vector<int> const &cut);

double lowConductanceCut(ArrayGraph::Graph *graph, std::vector<int> const &cut);

double normalizedCut(ArrayGraph::Graph *graph, std::vector<int> const &cut,
                     int numParts);

#endif