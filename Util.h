#pragma once
#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Shape.h"
#include "Cluster.h"

// Structure to represent a color
using Color_eigen = Eigen::RowVector3d;

// helper functions
void print2DMatrix(vector<vector<int>> & matrix);
Cluster translateToOrigin(Cluster& cluster);
Cluster translateToCenterOfOtherCluster(Cluster& cluster, Cluster& otherCluster);
Cluster copyPaste(Cluster& cluster, Cluster& otherCluster, Eigen::MatrixXi& best_rotation_matrix);
vector<Piece> getNeighbouringAndCurrentPieces(Piece& piece, Shape & shape);
vector<Piece> getNeighbouringAndCurrentPieces_random(Piece& piece, Shape & shape);
vector<Pixel> getNeighbouringPixels(const Piece& piece, const Shape& shape);
vector<Pixel> getNeighbouringPixels(const Pixel& pixel, const Shape& shape);
vector<Pixel> getNeighbouringPixels(const Pixel& pixel, const Shape& shape, int & distinctPieceNum);
vector<Pixel> getPossibleNeighbours(const Pixel& pixel);
vector<Pixel> getPossibleCircleNeigbours(const Pixel & pixel);
vector<Pixel> getValidNeigbours(Pixel & pixel, Shape & shape);
vector<Pixel> getUnoccupiedNeigbouringPixels(Pixel & pixel, Shape & shape);
vector<Pixel> getValidCircleNeigbours(Pixel & pixel, Shape & shape);
int getPieceIndex(const Pixel &pixel, const Shape &shape);
bool isPixelInPieces(const Pixel& pixel, const vector<Piece>& pieces);
bool isPixelInPieces(const Pixel& pixel, const vector<Piece>& pieces, int & pieceID);
int getAccessibilityValue(const Pixel& pixel,const Shape& shape);
int getBlockabilityValue(const Pixel& pixel,const Shape& shape);
void getTileInstanceList(const Pixel& pixel,const Shape& shape, set<int> & pieceIDList);
float getAvgAccessibilityValue(Shape & shape);
Pixel getPixelWithLeastAccessibilityValue(const vector<Pixel>& pixels, const Shape& shape);
Pixel getRandomPixelBasedOnAccessibilityValue(const vector<Pixel>& pixels, const Shape& shape);
Pixel getRandomPixel(const vector<Pixel>& pixels);
int getRandomNumber(int begin, int end);
bool isSameShape_old(const Piece& piece1, const Piece& piece2);
vector<vector<int>> rotate90(const vector<vector<int>>& original);
bool isSameShape(Piece & piece1, Piece & piece2);
bool isSubset(const vector<vector<int>>& subset, const vector<vector<int>>& superset);
bool isSubsetShape(Piece & piece1, Piece & piece2);
double distBetweenPx(const Pixel& px1, const Pixel& px2);
int distBetweenPx(Pixel & px1, Pixel & px2, Shape & shape);
bool isBlocked(const Piece& piece, const Shape& shape);
bool isPieceinTemplates(Piece& piece, vector<Piece>& templates);
bool isPieceIsSubsetShapeinTemplates(Piece & piece, vector<Piece> & templates);
int differenceInPixels(Cluster& cluster, Cluster& otherCluster);
int calculateMinDistance(Cluster& cluster, Cluster& otherCluster, Eigen::MatrixXi* bestRotationMatrix);
int getCurrentTemplateInstanceNum(Piece & current_template_piece, Shape & shape);
Piece getTemplateWithLargestGroup(vector<Piece> & template_pieces, Shape & shape);
vector<Piece> getTemplateInstances(Piece & current_template_piece, Shape & shape);
vector<Piece> getTemplateInstances(Piece & current_template_piece, Shape & shape, vector<int> & pieceIDList);
vector<Piece> getBlockedTemplateInstances(Piece & current_template_piece, Shape & shape);
vector<Piece> getUnblockedTemplateInstances(Piece & current_template_piece, Shape & shape);
bool isPieceInPieceList(Piece & piece, vector<Piece> & pieces);
bool isPixelInPixelList(Pixel & pixel, vector<Pixel> & pixels);
bool isPixelInPixelList(Pixel & pixel, vector<pair<Pixel, int>> & pixels);
Piece getPieceFromPixel(Pixel & pixel, vector<Piece> & pieces);
vector<Pixel> getAlonePixels(Shape & shape, vector<int> & pieceIDs);
int getPixelAccessibility(Pixel & pixel, Shape & shape, int depth);
int getPixelBlockability(Pixel & pixel, Shape & shape, int depth);
void updateAccessibilityAndBlockability(Shape & shape);
void updateAccessibility(Shape & shape);
void updateBlockability(Shape & shape);

// shape i/o operations
void saveShape2File(Shape & shape, string folderPath, string filePath, double runningTime);
Shape readShapeFromFile(string filePath);
Shape readSeedFromFile(string filePath);

// random staff
vector<int> GetRandomObjIndexList(vector<float> possibList, float alpha, int objNum);
int GetRandomObjIndex(vector<float> possibList, float alpha);
vector<float> PossibExpMapping(vector<float> possibList, float alpha);
float GetRandomNumber(int seedIndex);

// normalization
void normalizeVector(std::vector<float>& vec);

// randomly generate N colors
void generateUnsaturatedColor(Color_eigen& color, double minSaturation, double maxSaturation, double minValue, double maxValue, const std::vector<Color_eigen>& existingColors) ;


#endif