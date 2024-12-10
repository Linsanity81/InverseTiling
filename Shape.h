#pragma once
#ifndef SHAPE_H
#define SHAPE_H

#include <set>
#include "Piece.h"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

class Shape {

public:
    int volRow;
    int volCol;

    bool isDisconnected;
    bool valid;
    int score;
    int backtrack_count;
    float maxAvgSizeRatio;
    float minAvgSizeRatio;
    int maxTileSize;
    int minTileSize;

    vector<float> accessibilityList;
    vector<float> blockabilityList;
    vector<float> closenessCentralityList;
    vector<float> betweenessCentralityList;
    float remainingGrowingEvaluation;

    set<Pixel> all_pixels;
    vector<Piece> pieces;
    vector<Pixel> assignedPixels;

    vector<Piece> current_templates;
    vector<Piece> blocked_templates;

    vector< vector<Piece> > pieceListsByCluster;

    Eigen::MatrixXd V_grid;
    Eigen::MatrixXi F_grid;

    string growthLog;

    Shape();
    Shape(const set<Pixel> &all_pixels);
    Shape(const set<Pixel> &all_pixels, const int row, const int col);
    Shape(const vector<Piece>& pieces);
    Shape(const vector<Piece>& pieces, const int row, const int col);

    void clear();

    void addPiece(const Piece& piece);
    bool hasPixel(const Pixel& pixel) const;
    bool pixelInPieces(const Pixel & pixel);
    void assignColors();
    int getAllPiecesSize();
    vector<Pixel> getUnfilledPixels();
    int getNumUnfilledPx();
    void addLogGrowPiece(int piece_idx, Color piece_color);
    void addLogCurrentShapeState();
    void printTemplates();
    string printAllPixels();

    void assignPieceID();
    void assignPieceClusterID();
    std::vector<std::vector<int>> get2DPuzzleMatrix();
    std::vector<std::vector<int>> get2DPuzzleTemplateMatrix();
    string getOutputFolderPath(string prefix);
    string getOutputFileFullName(string prefix);

    int getSmallestPieceSize();
    int getLargestPieceSize();
    bool isSatisfyMinSizeRequirement();
    int getPixelID(const Pixel & pixel);
    int getPieceID(const Piece & piece);
    int getSmallSizePiecesNum(int currMinTileSize);
    int getLargeSizePiecesNum(int currMaxTileSize);
    double getInstanceNumDeviation();
    double getAverageEnlargeability();

    void updateTemplates(bool isLooseSizeRequirement);
    void updateCentrality();

    void saveEachPieceOBJs(string prefix);
    void saveEachClassTemplate(string prefix);

    void generateGridMesh();
};



#endif