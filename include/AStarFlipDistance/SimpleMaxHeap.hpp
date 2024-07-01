

#ifndef A_STAR_FOR_FLIPDISTANCE_SIMPLEMAXHEAP_HPP
#define A_STAR_FOR_FLIPDISTANCE_SIMPLEMAXHEAP_HPP

#endif //A_STAR_FOR_FLIPDISTANCE_SIMPLEMAXHEAP_HPP

#include <iostream>
#include <vector>
#include <algorithm>


class MaxHeap {
private:
    std::vector<std::pair<int,int>> heap;
    std::vector<int> node_positions;

    void heapifyUp(int idx);

    void heapifyDown(int idx);

public:

    MaxHeap()=default;

    MaxHeap(int maxNodes);

    void insert(int node, int key);

    void changeKey(int node, int newKey);

    void increaseKey(int node, int newKey);

    void decreaseKey(int node, int newKey);

    std::pair<int,int> extractMax();

    bool empty();

    bool contains(int node);
};

