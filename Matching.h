#pragma once
#ifndef MATCHING_H
#define MATCHING_H
#include "Cluster.h"
#include <map>

class Matching {

public:
    int dist;
    std::vector<Cluster> list1;
    std::vector<Cluster> list2;
    // contains pairs of indexes into list1 and list2 to represent matching betweem clusters
    std::vector<std::pair<int, int>> matching;
    std::map<int, Eigen::MatrixXi> best_rotations;
    Matching(std::vector<Cluster> list1, std::vector<Cluster> list2, std::vector<std::pair<int, int>> matching);
    void addPair(std::pair<Cluster,Cluster> pair);
    int calculateDistance();
    void assignColors();
    bool isSolution();
    bool operator < (const Matching& other);
};

#endif