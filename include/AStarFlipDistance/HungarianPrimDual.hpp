

#ifndef A_STAR_FOR_FLIPDISTANCE_HUNGARIANPRIMDUAL_HPP
#define A_STAR_FOR_FLIPDISTANCE_HUNGARIANPRIMDUAL_HPP


#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <climits>
#include "Global.hpp"

/*
 * This code is a port of the code by Kevin L.Stern with additional dynamic features
 * https://github.com/KevinStern/software-and-algorithms
 */

 /* Copyright (c) 2012 Kevin L. Stern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * An implementation of the Hungarian algorithm for solving the assignment
 * problem. An instance of the assignment problem consists of a number of
 * workers along with a number of jobs and a cost matrix which gives the cost of
 * assigning the i'th worker to the j'th job at position (i, j). The goal is to
 * find an assignment of workers to jobs so that no job is assigned more than
 * one worker and so that no worker is assigned to more than one job in such a
 * manner so as to minimize the total cost of completing the jobs.
 * <p>
 *
 * An assignment for a cost matrix that has more workers than jobs will
 * necessarily include unassigned workers, indicated by an assignment value of
 * -1; in no other circumstance will there be unassigned workers. Similarly, an
 * assignment for a cost matrix that has more jobs than workers will necessarily
 * include unassigned jobs; in no other circumstance will there be unassigned
 * jobs. For completeness, an assignment for a square cost matrix will give
 * exactly one unique worker to each job.
 * <p>
 *
 * This version of the Hungarian algorithm runs in time O(n^3), where n is the
 * maximum among the number of workers and the number of jobs.
 *
 * @author Kevin L. Stern
 */




class HungarianPrimDual {
private:
    std::vector<std::vector<triangle_int>> costMatrix;
    triangle_int rows, cols, dim;
    std::vector<triangle_int> labelByWorker, labelByJob;
    std::vector<triangle_int> minSlackWorkerByJob;
    std::vector<triangle_int> minSlackValueByJob;
    std::vector<triangle_int> matchJobByWorker, matchWorkerByJob;
    std::vector<triangle_int> parentWorkerByCommittedJob;
    std::vector<bool> committedWorkers;

    std::vector<triangle_int> in_which_line_is_the_edge;
    std::vector<triangle_int> parent_labelByWorker, parent_labelByJob;
    std::vector<triangle_int> parent_matchJobByWorker, parent_matchWorkerByJob;

    std::vector<triangle_int> edges_used_in_cost_matrix;



public:
    explicit HungarianPrimDual(std::vector<std::vector<triangle_int>>& _costMatrix);

    HungarianPrimDual()=default;

    void computeInitialFeasibleSolution();

    std::vector<triangle_int> execute();


    triangle_int parent_execute(const std::vector<int>& edges,triangle_int number_of_overall_edges);

    triangle_int child_execute(triangle_int edge_given_by_parent,const std::vector<triangle_int>& new_costs_in_row);

    void executePhase();

    triangle_int fetchUnmatchedWorker();

    /**
     * Find a valid matching by greedily selecting among zero-cost matchings. This
     * is a heuristic to jump-start the augmentation algorithm.
     */
    void greedyMatch();

    void initializePhase(triangle_int w);

    /**
     * Helper method to record a matching between worker w and job j.
     */
    void match(triangle_int w, triangle_int j);

    /**
     * Reduce the cost matrix by subtracting the smallest element of each row from
     * all elements of the row as well as the smallest element of each column from
     * all elements of the column. Note that an optimal assignment for a reduced
     * cost matrix is optimal for the original cost matrix.
     */
    void reduce();

    /**
     * Update labels with the specified slack by adding the slack value for
     * committed workers and by subtracting the slack value for committed jobs. In
     * addition, update the minimum slack values appropriately.
     */
    void updateLabeling(triangle_int slack);
};

#endif //A_STAR_FOR_FLIPDISTANCE_HUNGARIANPRIMDUAL_HPP
