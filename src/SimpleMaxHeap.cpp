
#include "AStarFlipDistance/SimpleMaxHeap.hpp"

void MaxHeap::heapifyUp(int idx) {
    while (idx > 0) {
        int parentIdx = (idx - 1) / 2;
        if (heap[parentIdx].first < heap[idx].first) {
            std::swap(heap[parentIdx], heap[idx]);
            std::swap(node_positions[heap[parentIdx].second], node_positions[heap[idx].second]);
            idx = parentIdx;
        } else {
            break;
        }
    }
}

void MaxHeap::heapifyDown(int idx) {
    int leftChildIdx = 2 * idx + 1;
    int rightChildIdx = 2 * idx + 2;
    int largest = idx;

    if (leftChildIdx < heap.size() && heap[leftChildIdx].first > heap[largest].first) {
        largest = leftChildIdx;
    }
    if (rightChildIdx < heap.size() && heap[rightChildIdx].first > heap[largest].first) {
        largest = rightChildIdx;
    }

    if (largest != idx) {
        std::swap(heap[idx], heap[largest]);
        std::swap(node_positions[heap[idx].second], node_positions[heap[largest].second]);
        heapifyDown(largest);
    }
}

void MaxHeap::insert(int node, int key) {
    if(node_positions[node]!=-1){
        return;
    }
    heap.emplace_back(key,node);
    node_positions[node] = heap.size() - 1;
    heapifyUp(heap.size() - 1);
}

void MaxHeap::changeKey(int node, int newKey) {
    auto tmp= heap[node_positions[node]];
    if(newKey == tmp.first){
        return;
    }
    if (node_positions[node] != -1) {
        if (newKey >= tmp.first) {
            increaseKey(node, newKey);
        } else if (newKey <= tmp.first) {
            decreaseKey(node, newKey);
        }
    }
    if(node_positions[node] == -1){
        insert(node,newKey);
    }
}

void MaxHeap::increaseKey(int node, int newKey) {
    auto tmp= heap[node_positions[node]];
    if (node_positions[node] != -1 && newKey > tmp.first) {
        int idx = node_positions[node];
        heap[idx] = std::pair(newKey,node);
        heapifyUp(idx);
    }
}

void MaxHeap::decreaseKey(int node, int newKey) {
    auto tmp= heap[node_positions[node]];
    if (node_positions[node] != -1 && newKey < tmp.first) {
        int idx = node_positions[node];
        heap[idx] = std::pair(newKey,node);
        heapifyDown(idx);
    }
}

std::pair<int, int> MaxHeap::extractMax() {
    auto max = heap[0];
    node_positions[max.second] = -1;
    heap[0] = heap.back();
    heap.pop_back();
    node_positions[heap[0].second] = 0;
    heapifyDown(0);
    return max;
}

bool MaxHeap::empty() {
    return heap.empty();
}

bool MaxHeap::contains(int node) {
    return node_positions[node] != -1;
}

MaxHeap::MaxHeap(int maxNodes) : node_positions(maxNodes, -1) {}
