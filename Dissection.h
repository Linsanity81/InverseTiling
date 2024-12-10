#pragma once
#ifndef DISSECTION_H
#define DISSECTION_H

#include "Shape.h"

using namespace std;

Shape generatePiecesLocally_distribution(Shape & shape, int pieceNumThreshold);
Shape generatePiecesLocally(Shape & shape, int pieceNumThreshold);
Shape optimizePiecesSizeLocally(Shape & shape, int minTileSize, int maxTileSize);
Shape getSubShape(vector<Piece> & eliminationPieces, Shape & shape, int & pieceNum, int & K, vector<int> & eliminationPieceIDList );
Shape generatePieces(Shape & shape, int N, int K, int G, int T, bool autoSave, string pureFileName, 
                     int seedNum, float minAvgSizeRatio, float maxAvgSizeRatio, int minTileSize, int maxTileSize,
                     int minSeedDist, bool isUniformDistribution, int shapeCandisNum,
                     double & runningTime, int & bestScore, vector<Piece> & existingTemplates, bool & isFound, bool isLooseSizeRequirement);
Shape generatePiecesWithGivenSeeds(Shape & shape, int N, int K, int G, int T, bool autoSave, string pureFileName, 
                     int seedNum, float minAvgSizeRatio, float maxAvgSizeRatio, int minTileSize, int maxTileSize,
                     int minSeedDist, bool isUniformDistribution, int shapeCandisNum,
                     double & runningTime, int & bestScore, vector<Piece> & existingTemplates, bool & isFound, bool isLooseSizeRequirement);
Shape generateNSeedPieces(int N, Shape & shape, int minSeedDist = 2, bool isUniformDistribution = true);
Shape growAllPiecesByOnePixel(const Shape& shape, bool use_access_val, vector<Shape>& shape_state);
Shape growPiecesByOnePixel(Shape & shape, int backtrackNum, bool & isDisconnected, vector<Piece> & existingTemplates, int minTileSize, int maxTileSize, bool isLooseSizeRequirement, vector<Shape> & growingStates);
vector<Pixel> getPossiblePixelsToGrow(Piece & piece, Shape & shape);
Piece growPiece(const Shape& shape, int piece_index, bool use_access_val);
Piece growPiece(const Shape& shape, Piece & piece);
Piece growFirstSelectedPiece(Shape& shape, vector<Piece> & current_templates, Piece & piece, float alpha);
vector<Pixel> getPossiblePixelsGrow2Match(Shape & shape, Piece & piece, Piece & piece_template);
Piece growPiece2Match(Shape & shape, Piece & piece, Piece & piece_template, bool & isMatched, float alpha);
Piece growPieceWithSameShape(const Shape& shape, int piece_index, vector<Pixel> neighbouring_pixels, const Piece& existingPiece, bool use_access_val);
void fillDisconnectedPixels(Shape & shape);
void updateTemplates(Shape & shape);

// selection functions
Piece selectGuidingPieceClassTemplate(vector<Piece> & piece_templates, Shape & shape, float alpha, bool & isDisconnected);
Piece selectFirstGrowingPieceInGuidingPieceClass(vector<Piece> & guiding_piece_instances, Shape & shape, float alpha);
Piece selectEliminationPieceTemplate(vector<Piece> & piece_templates, Shape & shape, int pieceNumThreshold, bool & noPossiblePieceTemplate);

#endif