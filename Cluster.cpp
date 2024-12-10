#include "Cluster.h"
#include <random>
#include <queue>
#include <Eigen/Geometry> 

Cluster::Cluster() {

}


Cluster::Cluster(std::vector<std::pair<int, int>> pixels) {
    initialize(pixels);
}

Cluster::Cluster(Eigen::MatrixXi pixelsAsMatrix) {
    initialize(pixelsAsMatrix);
}

void Cluster::initialize(std::vector<std::pair<int, int>> pixels) {
    this->pixels = pixels;
    this->pixelsAsMatrix = Eigen::MatrixXi(pixels.size(), 2);
    int row = 0;
    for (auto & pixel : pixels) {
        this->pixelsAsMatrix.row(row) << pixel.first, pixel.second;
        row++;
    }
}

void Cluster::initialize(Eigen::MatrixXi pixelsAsMatrix) {
    this->pixelsAsMatrix = pixelsAsMatrix;
    for (int i = 0; i < pixelsAsMatrix.rows(); i++) {
        this->pixels.push_back(std::pair<int, int>(pixelsAsMatrix.coeff(i, 0), pixelsAsMatrix.coeff(i, 1)));
    }
}

void Cluster::setColor(int r, int g, int b) {
    this->r = r;
    this->g = g;
    this->b = b;
}

void Cluster::addPixels(std::vector<std::pair<int, int>>& pixels) {
    for (auto px : pixels) {
        this->pixels.push_back(px);
    }
    initialize(this->pixels);
}

void Cluster::removePixels(std::vector<std::pair<int, int>> & pixels) {
    for (auto pixel : pixels) {
        auto it = std::find(this->pixels.begin(), this->pixels.end(), pixel);

        if (it != this->pixels.end()) {
            this->pixels.erase(it);
        }
    }
    initialize(this->pixels);


}

std::pair<int, int> Cluster::getBoundingBoxCenter() {
    int maxY = 0;
    int minY = 10;
    int maxX = 0;
    int minX = 10;

    for (auto pixel : pixels) {
        if (pixel.first > maxX) {
            maxX = pixel.first;
        }

        if (pixel.first < minX) {
            minX = pixel.first;
        }

        if (pixel.second > maxY) {
            maxY = pixel.second;
        }

        if (pixel.second < minY) {
            minY = pixel.second;
        }
    }

    int height = maxY - minY + 1;
    int width = maxX - minX + 1;

    return std::pair<int, int>(minX + (width / 2), minY + (height / 2));
}

void Cluster::generateRenderMeshes() {
    // render each pixel as two triangles that make up a square
    // so this generates 4 points and 2 faces for each pixel

    V = Eigen::MatrixXd(pixels.size() * 4, 2);
    V1 = Eigen::MatrixXd(pixels.size(), 2);
    V2 = Eigen::MatrixXd(pixels.size(), 2);
    V3 = Eigen::MatrixXd(pixels.size(), 2);
    V4 = Eigen::MatrixXd(pixels.size(), 2);
    F = Eigen::MatrixXi(pixels.size() * 2, 3);

    int Vrow = 0;
    int Virow = 0;
    int Frow = 0;
    for (auto pixel : pixels) {
        V.row(Vrow) << pixel.first - 0.5, pixel.second - 0.5;
        V1.row(Virow) << pixel.first - 0.5, pixel.second - 0.5;

        V.row(Vrow + 1) << pixel.first - 0.5, pixel.second + 0.5;
        V2.row(Virow) << pixel.first - 0.5, pixel.second + 0.5;

        V.row(Vrow + 2) << pixel.first + 0.5, pixel.second + 0.5;
        V3.row(Virow) << pixel.first + 0.5, pixel.second + 0.5;

        V.row(Vrow + 3) << pixel.first + 0.5, pixel.second - 0.5;
        V4.row(Virow) << pixel.first + 0.5, pixel.second - 0.5;

        F.row(Frow) << Vrow, Vrow + 1, Vrow + 2;
        F.row(Frow + 1) << Vrow, Vrow + 3, Vrow + 2;

        Vrow += 4;
        Virow++;
        Frow += 2;
    }

    // color matrix
    C = Eigen::MatrixXd(pixels.size() * 2, 3);
    for (int i = 0; i < pixels.size() * 2; i++) {
        C.row(i) << r, g, b;
    }
}

