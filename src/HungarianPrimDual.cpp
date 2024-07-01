
#include "AStarFlipDistance/HungarianPrimDual.hpp"

HungarianPrimDual::HungarianPrimDual(std::vector<std::vector<triangle_int>> &_costMatrix) {
    this->dim = std::max(_costMatrix.size(), _costMatrix[0].size());
    this->rows = _costMatrix.size();
    this->cols = _costMatrix[0].size();
    costMatrix=_costMatrix;


    labelByWorker.resize(this->dim, 0.0);
    labelByJob.resize(this->dim, 0.0);
    minSlackWorkerByJob.resize(this->dim, 0);
    minSlackValueByJob.resize(this->dim, 0.0);
    committedWorkers.resize(this->dim, false);
    parentWorkerByCommittedJob.resize(this->dim, 0);
    matchJobByWorker.resize(this->dim, -1);
    matchWorkerByJob.resize(this->dim, -1);
}

void HungarianPrimDual::computeInitialFeasibleSolution() {
    for (triangle_int j = 0; j < dim; j++) {
        labelByJob[j] = 10000;
    }
    for (triangle_int w = 0; w < dim; w++) {
        for (triangle_int j = 0; j < dim; j++) {
            if (costMatrix[w][j] < labelByJob[j]) {
                labelByJob[j] = costMatrix[w][j];
            }
        }
    }
}

void HungarianPrimDual::executePhase() {
    while (true) {
        triangle_int minSlackWorker = -1, minSlackJob = -1;
        triangle_int minSlackValue = 10000;
        for (triangle_int j = 0; j < dim; j++) {
            if (parentWorkerByCommittedJob[j] == -1) {
                if (minSlackValueByJob[j] < minSlackValue) {
                    minSlackValue = minSlackValueByJob[j];
                    minSlackWorker = minSlackWorkerByJob[j];
                    minSlackJob = j;
                }
            }
        }
        if (minSlackValue > 0) {
            updateLabeling(minSlackValue);
        }
        parentWorkerByCommittedJob[minSlackJob] = minSlackWorker;
        if (matchWorkerByJob[minSlackJob] == -1) {

            triangle_int committedJob = minSlackJob;
            triangle_int parentWorker = parentWorkerByCommittedJob[committedJob];
            while (true) {
                triangle_int temp = matchJobByWorker[parentWorker];
                match(parentWorker, committedJob);
                committedJob = temp;
                if (committedJob == -1) {
                    break;
                }
                parentWorker = parentWorkerByCommittedJob[committedJob];
            }
            return;
        } else {

            triangle_int worker = matchWorkerByJob[minSlackJob];
            committedWorkers[worker] = true;
            for (triangle_int j = 0; j < dim; j++) {
                if (parentWorkerByCommittedJob[j] == -1) {
                    triangle_int slack = costMatrix[worker][j] - labelByWorker[worker] - labelByJob[j];
                    if (minSlackValueByJob[j] > slack) {
                        minSlackValueByJob[j] = slack;
                        minSlackWorkerByJob[j] = worker;
                    }
                }
            }
        }
    }
}

std::vector<triangle_int> HungarianPrimDual::execute() {

    computeInitialFeasibleSolution();
    greedyMatch();

    {
        triangle_int w = fetchUnmatchedWorker();
        while (w < dim) {
            initializePhase(w);
            executePhase();
            w = fetchUnmatchedWorker();
        }
    }

    std::vector<triangle_int> result(matchJobByWorker.begin(), matchJobByWorker.begin() + rows);
    for (triangle_int w = 0; w < result.size(); w++) {
        if (result[w] >= cols) {
            result[w] = -1;
        }
    }

    return result;
}


triangle_int HungarianPrimDual::fetchUnmatchedWorker() {
    triangle_int w;
    for (w = 0; w < dim; w++) {
        if (matchJobByWorker[w] == -1) {
            break;
        }
    }
    return w;
}

void HungarianPrimDual::greedyMatch() {
    for (triangle_int w = 0; w < dim; w++) {
        for (triangle_int j = 0; j < dim; j++) {
            if (matchJobByWorker[w] == -1 && matchWorkerByJob[j] == -1 &&
                costMatrix[w][j] - labelByWorker[w] - labelByJob[j] == 0) {
                match(w, j);
            }
        }
    }
}

void HungarianPrimDual::initializePhase(triangle_int w) {
    std::fill(committedWorkers.begin(), committedWorkers.end(), false);
    std::fill(parentWorkerByCommittedJob.begin(), parentWorkerByCommittedJob.end(), -1);
    committedWorkers[w] = true;
    for (triangle_int j = 0; j < dim; j++) {
        minSlackValueByJob[j] = costMatrix[w][j] - labelByWorker[w] - labelByJob[j];
        minSlackWorkerByJob[j] = w;
    }
}

void HungarianPrimDual::match(triangle_int w, triangle_int j) {
    matchJobByWorker[w] = j;
    matchWorkerByJob[j] = w;
}

void HungarianPrimDual::reduce() {
    for (triangle_int w = 0; w < dim; w++) {
        triangle_int min = 10000;
        for (triangle_int j = 0; j < dim; j++) {
            if (costMatrix[w][j] < min) {
                min = costMatrix[w][j];
            }
        }
        for (triangle_int j = 0; j < dim; j++) {
            costMatrix[w][j] -= min;
        }
    }
    std::vector<triangle_int> min(dim, 10000);
    for (triangle_int w = 0; w < dim; w++) {
        for (triangle_int j = 0; j < dim; j++) {
            if (costMatrix[w][j] < min[j]) {
                min[j] = costMatrix[w][j];
            }
        }
    }
    for (triangle_int w = 0; w < dim; w++) {
        for (triangle_int j = 0; j < dim; j++) {
            costMatrix[w][j] -= min[j];
        }
    }
}

void HungarianPrimDual::updateLabeling(triangle_int slack) {
    for (triangle_int w = 0; w < dim; w++) {
        if (committedWorkers[w]) {
            labelByWorker[w] += slack;
        }
    }
    for (triangle_int j = 0; j < dim; j++) {
        if (parentWorkerByCommittedJob[j] != -1) {
            labelByJob[j] -= slack;
        } else {
            minSlackValueByJob[j] -= slack;
        }
    }
}


triangle_int HungarianPrimDual::parent_execute(const std::vector<int>& edges, triangle_int number_of_overall_edges) {
    computeInitialFeasibleSolution();
    greedyMatch();

    {
        triangle_int w = fetchUnmatchedWorker();
        while (w < dim) {
            initializePhase(w);
            executePhase();
            w = fetchUnmatchedWorker();
        }
    }

    std::vector<triangle_int> result(matchJobByWorker.begin(), matchJobByWorker.begin() + rows);
    for (triangle_int w = 0; w < result.size(); w++) {
        if (result[w] >= cols) {
            result[w] = -1;
        }
    }

    triangle_int value=0;
    for(triangle_int i=0;i<result.size();i++){
        value+=costMatrix[i][result[i]];
    }

    parent_labelByJob=std::vector<triangle_int>(labelByJob);
    parent_labelByWorker=std::vector<triangle_int>(labelByWorker);
    parent_matchJobByWorker=std::vector<triangle_int>(matchJobByWorker);
    parent_matchWorkerByJob=std::vector<triangle_int>(matchWorkerByJob);
    edges_used_in_cost_matrix.resize(number_of_overall_edges,-1);
    for(triangle_int i=0;i<edges.size();i++){
        edges_used_in_cost_matrix[edges[i]]=i;
    }
    return value;
}

triangle_int HungarianPrimDual::child_execute(triangle_int edge_given_by_parent, const std::vector<triangle_int> &new_costs_in_row) {

    labelByWorker=std::vector<triangle_int>(parent_labelByWorker);
    labelByJob=std::vector<triangle_int>(parent_labelByJob);
    matchJobByWorker=std::vector<triangle_int>(parent_matchJobByWorker);
    matchWorkerByJob=std::vector<triangle_int>(parent_matchWorkerByJob);

    triangle_int row_given_by_parent=edges_used_in_cost_matrix[edge_given_by_parent];


    triangle_int worker=row_given_by_parent;
    triangle_int job=matchJobByWorker[row_given_by_parent];
    matchWorkerByJob[matchJobByWorker[row_given_by_parent]]=-1;
    matchJobByWorker[row_given_by_parent]=-1;


    std::vector<triangle_int> old_cost=costMatrix[row_given_by_parent];
    costMatrix[row_given_by_parent]=new_costs_in_row;


    triangle_int mini=10000;
    for(triangle_int i=0;i<costMatrix[row_given_by_parent].size();i++){
        if(costMatrix[row_given_by_parent][i]-labelByJob[i]<mini){
            mini=costMatrix[row_given_by_parent][i]-labelByJob[i];
        }
    }
    labelByWorker[row_given_by_parent]=mini;






    {
        triangle_int alpha= costMatrix[worker][job]-labelByWorker[worker]-labelByJob[job];
        //std::cout<<alpha<<std::endl;
        if (alpha==0){
            match(worker,job);
        }
    }


    {
        triangle_int w = fetchUnmatchedWorker();
        while (w < dim) {
            initializePhase(w);
            executePhase();
            w = fetchUnmatchedWorker();
        }
    }



    triangle_int value=0;
    for(triangle_int i=0;i<matchJobByWorker.size();i++){
        value+=costMatrix[i][matchJobByWorker[i]];
    }

    costMatrix[row_given_by_parent]=old_cost;
    return value;
}

