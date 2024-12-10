#include "Matching.h"
#include "Util.h"

Matching::Matching(std::vector<Cluster> list1, std::vector<Cluster> list2, std::vector<std::pair<int, int>> matching) {
    this->list1 = list1;
    this->list2 = list2;
    this->matching = matching;
    calculateDistance();

}

void Matching::addPair(std::pair<Cluster, Cluster> pair) {
    list1.push_back(pair.first);
    list2.push_back(pair.second);
    matching.push_back(std::pair<int, int>(list1.size() - 1, list2.size() - 1));
    calculateDistance();
}

// returns the minimum distance between each matched pair, and stores the best rotation matrix for each pairing
int Matching::calculateDistance() {
    best_rotations = std::map<int, Eigen::MatrixXi>();
    int dist = 0;
    for (auto match: matching) {
        Eigen::MatrixXi best_rotation_matrix = Eigen::MatrixXi(2, 2);
        dist += calculateMinDistance(list1[match.first], list2[match.second], &best_rotation_matrix);
        best_rotations[match.first] = best_rotation_matrix;
    }

    this->dist = dist;

    return dist;
}


// assigns the same color for each cluster pair in this matching
void Matching::assignColors() {
    int colors[] = { 255, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255 , 0 };
    int colorIndex = 0;
    for (int i = 0; i < matching.size(); i++) {
        std::pair<int, int> cluster_pair = matching[i];
        list1[cluster_pair.first].setColor(colors[colorIndex], colors[colorIndex + 1], colors[colorIndex + 2]);
        list2[cluster_pair.second].setColor(colors[colorIndex], colors[colorIndex + 1], colors[colorIndex + 2]);
        colorIndex += 3;
    }
}

bool Matching::isSolution() {
    return calculateDistance() == 0;
}

bool Matching::operator < (const Matching& other) {
    return (this->dist < other.dist);
}