#pragma once
#ifndef CLUSTER_H
#define CLUSTER_H
#include <Eigen/Dense>
#include <vector>

class Cluster {

public:
    std::vector<std::pair<int, int>> pixels;
    Eigen::MatrixXi pixelsAsMatrix;
    int r; int g; int b;


    // these are all used by ligibl for visualization of the cluster
    Eigen::MatrixXd V;
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXd V3;
    Eigen::MatrixXd V4;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
    int viewer_mesh_id;

    Cluster();
    Cluster(std::vector<std::pair<int, int>> pixels);
    Cluster(Eigen::MatrixXi pixelsAsMatrix);
    void initialize(std::vector<std::pair<int, int>> pixels);
    void initialize(Eigen::MatrixXi pixelsAsMatrix);
    void setColor(int r, int g, int b);
    void addPixels(std::vector<std::pair<int, int>>& pixels);
    void removePixels(std::vector<std::pair<int, int>>& pixels);
    std::pair<int, int> getBoundingBoxCenter();
    void generateRenderMeshes();

private:

};



#endif