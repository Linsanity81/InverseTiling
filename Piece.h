#pragma once
#ifndef PIECE_H
#define PIECE_H
#include <Eigen/Dense>
#include <vector>
#include <set>

using namespace std;

struct Pixel {

    Pixel(int x, int y) {
        this->x = x;
        this->y = y;
    }

    int x;
    int y;
    int label;

    // float accessibility;
    // float betweenessCentrality;
    // float closenessCentrality;

    bool operator==(const Pixel &other) const {
        return (x == other.x && y == other.y);
    }

    bool operator<(const Pixel& other) const {
        return (x < other.x || (!(other.x < x) && y < other.y));
    }

    void print();
};

struct Color {
    int r;
    int g;
    int b;
};

class Piece {

public:
    bool isBlocked = false;

    vector<Pixel> pixels;
    Eigen::MatrixXi pixelsAsMatrix;
    Color color;
    vector<vector<int>> pieceShapeAsMatrix;

    int pieceID;
    int pieceClusterID;

    // these are all used by ligibl for visualization of the cluster
    Eigen::MatrixXd V;
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXd V3;
    Eigen::MatrixXd V4;
    Eigen::MatrixXi F;
    Eigen::MatrixXd C;
    int viewer_mesh_id;

    // these are for rendering boundary
    Eigen::MatrixXd V_boundary;
    Eigen::MatrixXi F_boundary;
    int viewer_boundary_mesh_id;

    // these are all used for saving obj files
    Eigen::MatrixXd V_obj;
    Eigen::MatrixXi F_obj;

    Piece();
    Piece(const Pixel& pixel);
    Piece(const Eigen::MatrixXi& pixelsAsMatrix, const Color& color);
    Piece(const vector<Pixel> &pixels,const Color &color);
    Piece(const vector<Pixel> &pixels,const Color &color, const int & pieceID);
    Piece(const vector<Pixel> &pixels);

    void initializeMatrix();
    void initializePieceShapeAsMatrix();
    void generateRenderMeshes();
    void generateSavingMeshes();
    bool hasPixel(const Pixel& pixel) const;
    void setColor(int r, int g, int b);
    void setColor(const Color &c);
    void addPixel(const Pixel& px);
    void deletePixel(const Pixel & px);
    Piece rotate(const Eigen::MatrixXi& rotation_matrix, int px_index) const;
    vector<Piece> get_all_possible_flips_and_rotations() const;
    bool isSamePiece(const Piece& other);
    string printPixels();
    void printPieceMatrix();

private:

};



#endif