#include "Util.h"
#include <random>
#include <cmath>

void print2DMatrix(vector<vector<int>> & matrix)
{
    // Get the number of rows and columns in the matrix
    int numRows = matrix.size();
    int numCols = matrix[0].size(); // Assuming all rows have the same number of columns

    // Print the 2D matrix
    std::cout << "Piece:" << std::endl;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    cout << endl;
}


Cluster translateToOrigin(Cluster & cluster) {
    std::pair<int, int> center = cluster.getBoundingBoxCenter();
    Eigen::Vector2i translation;
    translation << -center.first, -center.second;

    Eigen::MatrixXi translated_pixels = cluster.pixelsAsMatrix.rowwise() + translation.transpose();
    Cluster newCluster = Cluster(translated_pixels);
    newCluster.setColor(cluster.r, cluster.g, cluster.b);

    return newCluster;
}

Cluster translateToCenterOfOtherCluster(Cluster & cluster, Cluster& otherCluster) {
    std::pair<int, int> center = cluster.getBoundingBoxCenter();
    std::pair<int, int> otherCenter = otherCluster.getBoundingBoxCenter();

    Eigen::Vector2i translation;
    translation << otherCenter.first - center.first, otherCenter.second - center.second;


    Eigen::MatrixXi translated_pixels = cluster.pixelsAsMatrix.rowwise() + translation.transpose();
    return Cluster(translated_pixels);
}

//given the rotation matrix resulting in the lowest distance, this function returns a new cluster that is a copy-paste to otherCluster. 
Cluster copyPaste(Cluster & cluster, Cluster& otherCluster, Eigen::MatrixXi& best_rotation_matrix) {
    Cluster copyPastedCluster = Cluster(cluster.pixels);
    copyPastedCluster = translateToOrigin(copyPastedCluster);
    //rotate
    copyPastedCluster = Cluster((best_rotation_matrix * copyPastedCluster.pixelsAsMatrix.transpose()).transpose());
    //translate
    copyPastedCluster = translateToCenterOfOtherCluster(copyPastedCluster, otherCluster);

    //copy same color
    copyPastedCluster.setColor(cluster.r, cluster.g, cluster.b);

    return copyPastedCluster;
}

vector<Piece> getNeighbouringAndCurrentPieces(Piece& piece, Shape & shape)
{
    vector<Piece> pieceList;
    set<int> pieceIDList;

    pieceIDList.insert(shape.getPieceID(piece));

    for (int i = 0; i < piece.pixels.size(); ++i)
    {
        vector<Pixel> currNeigbourPixels = getPossibleNeighbours(piece.pixels[i]);

        for (int j = 0; j < currNeigbourPixels.size(); ++j)
        {
            int pieceID;

            if (shape.hasPixel(piece.pixels[i]) && isPixelInPieces(currNeigbourPixels[j], shape.pieces, pieceID))
            {
                pieceIDList.insert(pieceID);
            }
        }
    }

    for (auto & pieceID : pieceIDList)
    {
        pieceList.push_back(shape.pieces[pieceID]);
    }

    return pieceList;
}

vector<Piece> getNeighbouringAndCurrentPieces_random(Piece& piece, Shape & shape)
{
    vector<Piece> pieceList;
    set<int> pieceIDList;

    pieceIDList.insert(shape.getPieceID(piece));

    int centralPieceID = shape.getPieceID(piece);

    for (int i = 0; i < piece.pixels.size(); ++i)
    {
        vector<Pixel> currNeigbourPixels = getPossibleNeighbours(piece.pixels[i]);

        for (int j = 0; j < currNeigbourPixels.size(); ++j)
        {
            int pieceID;

            if (shape.hasPixel(piece.pixels[i]) && isPixelInPieces(currNeigbourPixels[j], shape.pieces, pieceID))
            {
                pieceIDList.insert(pieceID);
            }
        }
    }

    for (auto & pieceID : pieceIDList)
    {
        if (pieceID == centralPieceID)
        {
            pieceList.push_back(shape.pieces[pieceID]);
        }
        else
        {
            if (rand() % 2)
            {
                pieceList.push_back(shape.pieces[pieceID]);
            }
        }
    }

    return pieceList;  
}

// returns all pixels neighbouring a given piece 
// only returns pixels that are within the shape boundary and do not belong to other pieces in the shape
vector<Pixel> getNeighbouringPixels(const Piece& piece, const Shape& shape) {
    set<Pixel> neighbouring_pixels;

    for (const auto& px : piece.pixels) {
        vector<Pixel> neighbours = getNeighbouringPixels(px, shape);
        for (const  auto& pixel : neighbours) {
            neighbouring_pixels.insert(pixel);
        }
    }

    return vector<Pixel>(neighbouring_pixels.begin(), neighbouring_pixels.end());
}

// returns all neighbouring pixels that are within the boundary of the shape, and do not belong to other pieces
// assumes pixel belongs to one of the pieces in the shape
vector<Pixel> getNeighbouringPixels(const Pixel& pixel, const Shape& shape) {
    vector<Pixel> possible_pixels = getPossibleNeighbours(pixel);
    vector<Pixel> neighbouring_pixels;
    for (auto px : possible_pixels) {
        if (shape.hasPixel(px) && !isPixelInPieces(px, shape.pieces)) {
            neighbouring_pixels.push_back(px);
        } 
    }
    return neighbouring_pixels;
}

vector<Pixel> getNeighbouringPixels(const Pixel& pixel, const Shape& shape, int & distinctPieceNum)
{
    vector<Pixel> possible_pixels = getPossibleNeighbours(pixel);
    vector<Pixel> neighbouring_pixels;
    set<int> distinctPieceIDList;
    for (auto px : possible_pixels) {

        if (shape.hasPixel(px))
        {
            int currPieceID;
            if (isPixelInPieces(px, shape.pieces, currPieceID))
            {
                distinctPieceIDList.insert(currPieceID);
            }
            else
            {
                neighbouring_pixels.push_back(px);
            }
        } 
    }

    distinctPieceNum = distinctPieceIDList.size();
    return neighbouring_pixels;
}

vector<Pixel> getPossibleNeighbours(const Pixel& pixel) {
    return vector<Pixel>{
    Pixel(pixel.x + 1,pixel.y), Pixel(pixel.x, pixel.y - 1),
    Pixel(pixel.x - 1,pixel.y), Pixel(pixel.x,pixel.y + 1)
    };
}

vector<Pixel> getPossibleCircleNeigbours(const Pixel & pixel)
{
    return vector<Pixel>{
    Pixel(pixel.x + 1,pixel.y), Pixel(pixel.x, pixel.y - 1),
    Pixel(pixel.x - 1,pixel.y), Pixel(pixel.x,pixel.y + 1),
    Pixel(pixel.x - 1,pixel.y - 1), Pixel(pixel.x - 1,pixel.y + 1),
    Pixel(pixel.x + 1,pixel.y - 1), Pixel(pixel.x + 1,pixel.y + 1)
    };
}

vector<Pixel> getValidNeigbours(Pixel & pixel, Shape & shape)
{
    vector<Pixel> possible_pixels = getPossibleNeighbours(pixel);

    vector<Pixel> neighbouring_pixels;
    for (auto px : possible_pixels) {
        if (shape.hasPixel(px)) {
            neighbouring_pixels.push_back(px);
        } 
    }
    return neighbouring_pixels;
}

vector<Pixel> getUnoccupiedNeigbouringPixels(Pixel & pixel, Shape & shape)
{
    vector<Pixel> possible_pixels = getPossibleNeighbours(pixel);

    vector<Pixel> neighbouring_pixels;
    for (auto px : possible_pixels) {
        if (shape.hasPixel(px) and !shape.pixelInPieces(px)) {
            neighbouring_pixels.push_back(px);
        } 
    }
    return neighbouring_pixels;
}

vector<Pixel> getValidCircleNeigbours(Pixel & pixel, Shape & shape)
{
    vector<Pixel> possible_pixels = getPossibleCircleNeigbours(pixel);

    vector<Pixel> neighbouring_pixels;
    for (auto px : possible_pixels) {
        if (shape.hasPixel(px)) {
            neighbouring_pixels.push_back(px);
        } 
    }
    return neighbouring_pixels;
}


// returns the index of the piece in shape that pixel belongs to
// returns -1 if this pixel doesn't exist in the shape
int getPieceIndex(const Pixel& pixel, const Shape& shape) {
    for (int i = 0; i < shape.pieces.size(); i++) {
        if (shape.pieces[i].hasPixel(pixel)) {
            return i;
        }
    }

    return -1;
}

// checks if the given pixel belongs to any of the given pieces
bool isPixelInPieces(const Pixel& pixel, const vector<Piece> &pieces) {
    for (const auto& piece : pieces) {
        if (piece.hasPixel(pixel)) {
            return true;
        }
    }

    return false;
}

bool isPixelInPieces(const Pixel& pixel, const vector<Piece>& pieces, int & pieceID)
{
    pieceID = -1;

    int i = 0;
    for (const auto& piece : pieces) {
        if (piece.hasPixel(pixel)) {
            pieceID = i;
            return true;
        }

        ++i;
    }

    return false;
}

// returns random number inclusive of defined range
int getRandomNumber(int begin, int end) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(begin, end); // define the range
    return distr(gen);
}

// returns accesbility value of a pixel in a given shape
// right now, accessibility is just the number of 'free' neighbours the pixel has
// free neighbours are ones that do not not belong to any piece in the shape
int getAccessibilityValue(const Pixel& pixel, const Shape& shape) {
    vector<Pixel> possible_neighbours = getPossibleNeighbours(pixel);
    int free_neighbours_count = 0;

    for (const auto& px : possible_neighbours) {
        if (shape.hasPixel(px) && !isPixelInPieces(px, shape.pieces)) {
            free_neighbours_count++;
        }
    }

    return free_neighbours_count;
}

// returns the number of neigbouring pieces
int getBlockabilityValue(const Pixel& pixel,const Shape& shape)
{
    vector<Pixel> possible_neighbours = getPossibleNeighbours(pixel);
    set<int> pieceIDList;

    for (const auto& px : possible_neighbours) {
        if (shape.hasPixel(px) && isPixelInPieces(px, shape.pieces)) {
            pieceIDList.insert(getPieceIndex(px, shape));
        }
    }

    return pieceIDList.size();
}

void getTileInstanceList(const Pixel& pixel,const Shape& shape, set<int> & pieceIDList)
{
    vector<Pixel> possible_neighbours = getPossibleNeighbours(pixel);

    for (const auto& px : possible_neighbours) {
        if (shape.hasPixel(px) && isPixelInPieces(px, shape.pieces)) {
            pieceIDList.insert(getPieceIndex(px, shape));
        }
    }
}

float getAvgAccessibilityValue(Shape & shape)
{
    vector<Pixel> unfilled_pixels = shape.getUnfilledPixels();

    float totalAccess = 0;
    for (auto & pixel : unfilled_pixels)
    {
        vector<Pixel> possible_neighbours = getPossibleNeighbours(pixel);
        int free_neighbours_count = 0;

        for (const auto& px : possible_neighbours) 
            if (shape.hasPixel(px) && !isPixelInPieces(px, shape.pieces))
                free_neighbours_count++;
        

        totalAccess += free_neighbours_count;
    }

    return totalAccess / float(unfilled_pixels.size());
}


// returns the pixel with the least accessibility value among a given set of pixels in a shape
Pixel getPixelWithLeastAccessibilityValue(const vector<Pixel>& pixels, const Shape& shape) {
    int min_val = INT_MAX;
    int min_px_idx;
    for (int i = 0; i < pixels.size(); i++) {
        int accessibility_value = getAccessibilityValue(pixels[i], shape);
        if (accessibility_value < min_val) {
            min_val = accessibility_value;
            min_px_idx = i;
        }
    }

    return pixels[min_px_idx];
}

// similar to getPixelWithLeastAccessibilityValue, but chooses a pixel pseudo-randomly based on accesibility value
Pixel getRandomPixelBasedOnAccessibilityValue(const vector<Pixel>& pixels, const Shape& shape) {
    vector<Pixel> accessibility_val_px_array;

    for (const auto& px : pixels) {
        int accessibility_val = getAccessibilityValue(px, shape);
        int prob = 5 - accessibility_val; // highest possible access. val is 4, so scale the probability accordingly
        for (int i = 0; i < prob * prob; i++) {
            accessibility_val_px_array.push_back(px);
        }
    }

    return accessibility_val_px_array[getRandomNumber(0, accessibility_val_px_array.size() - 1)];
}

// used in a similar manner to getRandomPixelBasedOnAccessibilityValue, but this time purely returns a random pixel
Pixel getRandomPixel(const vector<Pixel>& pixels) {
    int idx = getRandomNumber(0, pixels.size() - 1);
    return pixels[idx];
}

bool isSameShape_old(const Piece& piece1, const Piece& piece2) {

    vector<Piece> all_possible_flips_and_rotations = piece1.get_all_possible_flips_and_rotations();

    for (const auto& piece : all_possible_flips_and_rotations) {
        // try all possible translations until there's a match
        for (int i = 0; i < piece1.pixels.size(); i++) {
            for (int j = 0; j < piece2.pixels.size(); j++) {
                // translation from pixel at i to pixel at j
                Eigen::Vector2i translation;
                translation << piece2.pixels[j].x - piece1.pixels[i].x, piece2.pixels[j].y - piece1.pixels[j].y;
                Eigen::MatrixXi translated_pixels = piece.pixelsAsMatrix.rowwise() + translation.transpose();
                Piece newPiece = Piece(translated_pixels, piece1.color);
                if (newPiece.isSamePiece(piece2)) {
                    return true;
                }
            }
        }
    }

    return false;
}

vector<vector<int>> rotate90(const vector<vector<int>>& original) {
    int rows = original.size();
    int cols = original[0].size();
    vector<vector<int>> rotated(cols, vector<int>(rows, 0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            rotated[j][rows - i - 1] = original[i][j];
            // rotated[j][i] = original[i][j];
        }
    }

    return rotated;
}

bool isSameShape(Piece& piece1, Piece& piece2)
{
    if (piece1.pixels.size() != piece2.pixels.size())
        return false;
    
    piece1.initializePieceShapeAsMatrix();
    piece2.initializePieceShapeAsMatrix();

    vector<vector<int>> A =  piece1.pieceShapeAsMatrix;
    vector<vector<int>> B =  piece2.pieceShapeAsMatrix;

    for (int i = 0; i < 5; ++i) {
        if (A == B || A == std::vector<std::vector<int>>(B.rbegin(), B.rend()) || A == std::vector<std::vector<int>>(B.begin(), B.end())) {
            return true;
        }

        B = rotate90(B);
    }

    return false;
}

bool isSubset(vector<vector<int>>& subset, vector<vector<int>>& superset)
{
    int rowsSubset = subset.size();
    int colsSubset = subset[0].size();
    int rowsSuperset = superset.size();
    int colsSuperset = superset[0].size();

    if (rowsSubset != rowsSuperset || colsSubset != colsSuperset) {
        return false;
    }

    for (int i = 0; i < rowsSubset; ++i)
    {
        for (int j = 0; j < colsSubset; ++j) 
        {
            if (subset[i][j] == 1 and superset[i][j] == 0)
                return false;
        }
    }

    return true;
}

// check if piece1 is the subset of piece2
// if piece1 == piece2, return false
bool isSubsetShape(Piece & piece1, Piece & piece2)
{
    piece1.initializePieceShapeAsMatrix();
    piece2.initializePieceShapeAsMatrix();

    vector<vector<int>> A =  piece1.pieceShapeAsMatrix;
    vector<vector<int>> B =  piece2.pieceShapeAsMatrix;

    for (int i = 0; i < 4; ++i)
    {
        if (isSubset(A, B) and A != B)
            return true;

        B = rotate90(B);
    }

    return false;
}

double distBetweenPx(const Pixel& px1, const Pixel& px2) {
    int dx = px2.x - px1.x;
    int dy = px2.y - px1.y;

    return sqrt(dx * dx + dy * dy);
}

int distBetweenPx(Pixel & px1, Pixel & px2, Shape & shape)
{
    int currDist = 0;

    vector<Pixel> visitedPixel;
    vector<Pixel> nextVisits;
    vector<int> currDepthList;

    nextVisits.push_back(px1);
    currDepthList.push_back(0);

    while (!nextVisits.empty())
    {
        int currDepth = currDepthList[0];
        Pixel currPixel = nextVisits[0];

        nextVisits.erase(nextVisits.begin());
        currDepthList.erase(currDepthList.begin());
        visitedPixel.push_back(currPixel);

        // vector<Pixel> newNeigbouringPixels = getValidCircleNeigbours(currPixel, shape);
        vector<Pixel> newNeigbouringPixels = getValidNeigbours(currPixel, shape);

        for (int i = 0; i < newNeigbouringPixels.size(); ++i)
        {
            if (newNeigbouringPixels[i] == px2)
                return (currDepth + 1);

            if(!isPixelInPixelList(newNeigbouringPixels[i], nextVisits) and !isPixelInPixelList(newNeigbouringPixels[i], visitedPixel))
            {
                nextVisits.push_back(newNeigbouringPixels[i]);
                currDepthList.push_back(currDepth + 1);
            }
        }
    }

    return 10000;
}

bool isBlocked(const Piece& piece, const Shape& shape) {
    return getNeighbouringPixels(piece, shape).size() == 0;
}

bool isPieceinTemplates(Piece& piece, vector<Piece>& templates) {
    for (auto& tmp : templates) {
        if (isSameShape(piece, tmp)) {
            return true;
        }
    }
    return false;
}

bool isPieceIsSubsetShapeinTemplates(Piece & piece, vector<Piece> & templates)
{
    for (auto & tmp : templates)
    {
        if (isSubsetShape(piece, tmp) and !isSameShape(piece, tmp))
            return true;
    }

    return false;
}

// finds total number of unmatched pixels between two clusters
int differenceInPixels(Cluster& cluster, Cluster& otherCluster) {
    std::set<std::pair<int, int>> clusterSet(cluster.pixels.begin(), cluster.pixels.end());
    std::set<std::pair<int, int>> otherClusterSet(otherCluster.pixels.begin(), otherCluster.pixels.end());

    std::set<std::pair<int, int>> union_set;
    std::set_union(clusterSet.begin(), clusterSet.end(),
        otherClusterSet.begin(), otherClusterSet.end(),
        std::inserter(union_set, union_set.begin()));

    std::set<std::pair<int, int>> intersection_set;
    std::set_intersection(clusterSet.begin(), clusterSet.end(),
        otherClusterSet.begin(), otherClusterSet.end(),
        std::inserter(intersection_set, intersection_set.begin()));

    return union_set.size() - intersection_set.size();
}


// function for calculating minimum distance between two clusters (after accounting for all valid rotations and flipping). 
// Stores the rotation matrix leading to the leading distance in bestRotationMatrix.
int calculateMinDistance(Cluster& cluster, Cluster& otherCluster, Eigen::MatrixXi* bestRotationMatrix) {
    Cluster origin_cluster = translateToOrigin(cluster);

    std::vector<Eigen::MatrixXi> all_valid_rotations_and_flips;

    Eigen::MatrixXi zero = Eigen::MatrixXi(2, 2);
    zero << 1, 0,
        0, 1;

    Eigen::MatrixXi ninety = Eigen::MatrixXi(2, 2);
    ninety << 0, -1,
        1, 0;

    Eigen::MatrixXi one_eighty = Eigen::MatrixXi(2, 2); //equivalent to a horizontal flip
    one_eighty << -1, 0,
        0, -1;

    Eigen::MatrixXi two_seventy = Eigen::MatrixXi(2, 2);
    two_seventy << 0, 1,
        -1, 0;

    Eigen::MatrixXi vertical_flip = Eigen::MatrixXi(2, 2);
    vertical_flip << -1, 0,
        0, 1;

    all_valid_rotations_and_flips.push_back(zero);
    all_valid_rotations_and_flips.push_back(ninety);
    all_valid_rotations_and_flips.push_back(one_eighty);
    all_valid_rotations_and_flips.push_back(two_seventy);
    all_valid_rotations_and_flips.push_back(vertical_flip);

    int minDist = INT_MAX;

    for (auto mat : all_valid_rotations_and_flips) {
        Eigen::MatrixXi rotated_matrix = (mat * origin_cluster.pixelsAsMatrix.transpose()).transpose();
        Cluster rotated_cluster = Cluster(rotated_matrix);
        Cluster rotated_and_translated_cluster = translateToCenterOfOtherCluster(rotated_cluster, otherCluster);
        int newDist = differenceInPixels(rotated_and_translated_cluster, otherCluster);
        if (newDist < minDist) {
            minDist = newDist;
            *bestRotationMatrix = mat;
        }
    }


    return minDist;
}

int getCurrentTemplateInstanceNum(Piece & current_template_piece, Shape & shape)
{
    int i = 0;

    for (auto & piece : shape.pieces)
    {
        if (isSameShape(piece, current_template_piece))
            i++;
    }

    return i;
}

Piece getTemplateWithLargestGroup(vector<Piece> & template_pieces, Shape & shape)
{
    Piece largestPieceGroupTemplate = template_pieces[0];
    int largestNum = getCurrentTemplateInstanceNum(template_pieces[0], shape);

    for (int i = 1; i < template_pieces.size(); ++i)
    {
        int currNum = getCurrentTemplateInstanceNum(template_pieces[i], shape);

        if (currNum > largestNum)
        {
            largestNum = currNum;
            largestPieceGroupTemplate = template_pieces[i];
        }
    }

    return largestPieceGroupTemplate;
}

vector<Piece> getTemplateInstances(Piece & current_template_piece, Shape & shape)
{
    vector<Piece> instances;

    for (auto & piece : shape.pieces)
    {
        if (isSameShape(piece, current_template_piece))
            instances.push_back(piece);
    }

    return instances;
}

vector<Piece> getTemplateInstances(Piece & current_template_piece, Shape & shape, vector<int> & pieceIDList)
{
    vector<Piece> instances;

    int i = 0;
    for (auto & piece : shape.pieces)
    {
        if (isSameShape(piece, current_template_piece))
        {
            instances.push_back(piece);
            pieceIDList.push_back(i);
        }
        ++i;
    }

    return instances;   
}

vector<Piece> getBlockedTemplateInstances(Piece & current_template_piece, Shape & shape)
{
    vector<Piece> instances;

    for (auto & piece : shape.pieces)
    {
        if (isSameShape(piece, current_template_piece) && isBlocked(piece, shape))
            instances.push_back(piece);
    }

    return instances;
}

vector<Piece> getUnblockedTemplateInstances(Piece & current_template_piece, Shape & shape)
{
    vector<Piece> instances;

    for (auto & piece : shape.pieces)
    {
        if (isSameShape(piece, current_template_piece) && !isBlocked(piece, shape))
            instances.push_back(piece);
    }

    return instances;
}

bool isPieceInPieceList(Piece & piece, vector<Piece> & pieces)
{
    for (auto & pe : pieces)
    {
        if (pe.isSamePiece(piece))
            return true;
    }

    return false;
}

bool isPixelInPixelList(Pixel & pixel, vector<Pixel> & pixels)
{
    for (auto & px : pixels)
    {
        if (pixel == px)
            return true;
    }

    return false;
}

bool isPixelInPixelList(Pixel & pixel, vector<pair<Pixel, int>> & pixels)
{
    for (auto & px : pixels)
    {
        if (pixel == px.first)
            return true;
    }

    return false;
}

Piece getPieceFromPixel(Pixel & pixel, vector<Piece> & pieces)
{
    Piece newPiece;

    for (const auto& piece : pieces) {
        if (piece.hasPixel(pixel)) {
            newPiece = piece;
            return newPiece;
        }
    }

    return newPiece;
}

vector<Pixel> getAlonePixels(Shape & shape, vector<int> & pieceIDs)
{
    vector<Pixel> unfilled_pixels = shape.getUnfilledPixels();
    // cout << "unfilled pixel number: " << unfilled_pixels.size() << endl;
    vector<Pixel> alonePixels;

    for (int k = 0; k < unfilled_pixels.size(); ++k)
    {
        vector<Pixel> visitedPixel;
        vector<Pixel> nextVisitPixels;

        nextVisitPixels.push_back(unfilled_pixels[k]);

        // vector<Pixel> neighbouringPixels = getNeighbouringPixels(pixel, shape);
        vector<int> neighbouringPieceIDs;

        while (!nextVisitPixels.empty())
        {
            Pixel currPixel = nextVisitPixels[0];
            nextVisitPixels.erase(nextVisitPixels.begin());

            visitedPixel.push_back(currPixel);

            if (isPixelInPieces(currPixel, shape.pieces))
            {
                int currID = getPieceIndex(currPixel, shape);

                if (std::find(neighbouringPieceIDs.begin(), neighbouringPieceIDs.end(), currID) == neighbouringPieceIDs.end())
                {
                    neighbouringPieceIDs.push_back(currID);

                    if (neighbouringPieceIDs.size() > 1)
                        break;
                }
            }
            else
            {
                vector<Pixel> newNeigbouringPixels = getValidNeigbours(currPixel, shape); 

                for (int i = 0; i < newNeigbouringPixels.size(); ++i)
                {
                    if(!isPixelInPixelList(newNeigbouringPixels[i], nextVisitPixels) and !isPixelInPixelList(newNeigbouringPixels[i], visitedPixel))
                    {
                        nextVisitPixels.push_back(newNeigbouringPixels[i]);
                    }
                } 
            }

            // cout << nextVisitPixels.size() << endl;
        }

        if (neighbouringPieceIDs.size() < 2)
        {
            alonePixels.push_back(unfilled_pixels[k]);
            pieceIDs.push_back(neighbouringPieceIDs[0]);
        }
    }

    return alonePixels;
}

int getPixelAccessibility(Pixel & pixel, Shape & shape, int depth)
{
    int pixelAccessibility = 0;
    vector<Pixel> visitedPixel;
    vector<Pixel> nextVisits;
    vector< pair<Pixel, int> > nextVisitPixels;

    nextVisitPixels.push_back(pair(pixel, 0));
    nextVisits.push_back(pixel);

    while (!nextVisitPixels.empty())
    {
        int currDepth = nextVisitPixels[0].second;
        Pixel currPixel = nextVisits[0];

        nextVisits.erase(nextVisits.begin());
        nextVisitPixels.erase(nextVisitPixels.begin());
        visitedPixel.push_back(currPixel);

        pixelAccessibility += getAccessibilityValue(currPixel, shape);

        vector<Pixel> newNeigbouringPixels = getValidNeigbours(currPixel, shape);

        for (int i = 0; i < newNeigbouringPixels.size(); ++i)
        {
            if (currDepth < depth)
            {
                if(!isPixelInPixelList(newNeigbouringPixels[i], nextVisits) and !isPixelInPixelList(newNeigbouringPixels[i], visitedPixel))
                {
                    nextVisitPixels.push_back(pair(newNeigbouringPixels[i], currDepth + 1));
                    nextVisits.push_back(newNeigbouringPixels[i]);
                }
            }
        } 
    }

    return pixelAccessibility;
}

int getPixelBlockability(Pixel & pixel, Shape & shape, int depth)
{
    int blockability = 0;

    vector<Pixel> visitedPixel;
    vector<Pixel> nextVisits;
    vector< pair<Pixel, int> > nextVisitPixels;

    nextVisitPixels.push_back(pair(pixel, 0));
    nextVisits.push_back(pixel);

    set<int> pieceIDList;

    while (!nextVisitPixels.empty())
    {
        int currDepth = nextVisitPixels[0].second;
        Pixel currPixel = nextVisits[0];

        nextVisits.erase(nextVisits.begin());
        nextVisitPixels.erase(nextVisitPixels.begin());
        visitedPixel.push_back(currPixel);

        // blockability += getBlockabilityValue(currPixel, shape);
        getTileInstanceList(currPixel, shape, pieceIDList);

        vector<Pixel> newNeigbouringPixels = getUnoccupiedNeigbouringPixels(currPixel, shape);

        for (int i = 0; i < newNeigbouringPixels.size(); ++i)
        {
            if (currDepth < depth)
            {
                if(!isPixelInPixelList(newNeigbouringPixels[i], nextVisits) and !isPixelInPixelList(newNeigbouringPixels[i], visitedPixel))
                {
                    nextVisitPixels.push_back(pair(newNeigbouringPixels[i], currDepth + 1));
                    nextVisits.push_back(newNeigbouringPixels[i]);
                }
            }
        } 
    }

    return pieceIDList.size();
}

void updateAccessibilityAndBlockability(Shape & shape)
{
    vector<float> accessibility(shape.volCol * shape.volRow, 0);
	vector<float> blockability(shape.volCol * shape.volRow, 0);

    vector<Pixel> unfilledPixels = shape.getUnfilledPixels();

    float minAccess = __FLT_MAX__;
    float minBlocka = __FLT_MAX__;
    float maxAccess = 0;
    float maxBlocka = 0;

	// compute
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
		float currAccessibility = getPixelAccessibility(unfilledPixels[i], shape, 3);
		float currBlockability = getPixelBlockability(unfilledPixels[i], shape, 2);

        if (currAccessibility < minAccess) minAccess = currAccessibility;
        if (currAccessibility > maxAccess) maxAccess = currAccessibility;
        if (currBlockability < minBlocka)  minBlocka = currBlockability;
        if (currBlockability > maxBlocka)  maxBlocka = currBlockability;
        
        int currPixelID = shape.getPixelID(unfilledPixels[i]);
		accessibility[currPixelID] = currAccessibility;
		blockability[currPixelID] = currBlockability;

        // unfilledPixels[i].print();
        // cout << currPixelID << " " << currAccessibility << " " << currBlockability << endl;
	}

    normalizeVector(accessibility);
    normalizeVector(blockability);

	// convert to possibility
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
        int currPixelID = shape.getPixelID(unfilledPixels[i]);

        // cout << accessibility[currPixelID] << endl;

        if (minAccess == maxAccess)
        {
            accessibility[currPixelID] = 10.0;
        }
        else 
        {
            accessibility[currPixelID] = round( 10.0 * (accessibility[currPixelID] - minAccess) / (maxAccess - minAccess) );
        }

        if (accessibility[currPixelID] < 1.0)
            accessibility[currPixelID] = 1.0;

        if (minBlocka == maxBlocka)
        {
            blockability[currPixelID] = 10.0;
        }
        else
        {
            blockability[currPixelID] = round( 10.0 - 10.0 * (blockability[currPixelID] - minBlocka) / (maxBlocka - minBlocka));
        }

        if (blockability[currPixelID] < 1.0)
            blockability[currPixelID] = 1.0;
	}

	// store
	shape.accessibilityList = accessibility;
	shape.blockabilityList = blockability;
}

void updateAccessibility(Shape & shape)
{
    vector<float> accessibility(shape.volCol * shape.volRow, 0);
    vector<Pixel> unfilledPixels = shape.getUnfilledPixels();

    float minAccess = __FLT_MAX__;
    float minBlocka = __FLT_MAX__;
    float maxAccess = 0;
    float maxBlocka = 0;

	// compute
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
		float currAccessibility = getPixelAccessibility(unfilledPixels[i], shape, 2);

        if (currAccessibility < minAccess) minAccess = currAccessibility;
        if (currAccessibility > maxAccess) maxAccess = currAccessibility;
        
        int currPixelID = shape.getPixelID(unfilledPixels[i]);
		accessibility[currPixelID] = currAccessibility;

        // unfilledPixels[i].print();
        // cout << currPixelID << " " << currAccessibility << " " << currBlockability << endl;
	}

    normalizeVector(accessibility);

	// convert to possibility
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
        int currPixelID = shape.getPixelID(unfilledPixels[i]);

        // cout << accessibility[currPixelID] << endl;

        // accessibility[currPixelID] = 1 - accessibility[currPixelID];

        if (minAccess == maxAccess)
        {
            accessibility[currPixelID] = 10.0;
        }
        else 
        {
            accessibility[currPixelID] = round( 10.0 * (accessibility[currPixelID] - minAccess) / (maxAccess - minAccess) );
        }

        if (accessibility[currPixelID] < 1.0)
            accessibility[currPixelID] = 1.0;
	}

	// store
	shape.accessibilityList = accessibility;
}

void updateBlockability(Shape & shape)
{
	vector<float> blockability(shape.volCol * shape.volRow, 0);

    vector<Pixel> unfilledPixels = shape.getUnfilledPixels();

    float minAccess = __FLT_MAX__;
    float minBlocka = __FLT_MAX__;
    float maxAccess = 0;
    float maxBlocka = 0;

	// compute
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
		float currBlockability = getPixelBlockability(unfilledPixels[i], shape, 2);

        if (currBlockability < minBlocka)  minBlocka = currBlockability;
        if (currBlockability > maxBlocka)  maxBlocka = currBlockability;
        
        int currPixelID = shape.getPixelID(unfilledPixels[i]);
		blockability[currPixelID] = currBlockability;

        // unfilledPixels[i].print();
        // cout << currPixelID << " " << currAccessibility << " " << currBlockability << endl;
	}

    normalizeVector(blockability);

	// convert to possibility
	for (int i = 0; i < unfilledPixels.size(); ++i)
	{
        int currPixelID = shape.getPixelID(unfilledPixels[i]);

        // cout << accessibility[currPixelID] << endl;

        if (minBlocka == maxBlocka)
        {
            blockability[currPixelID] = 10.0;
        }
        else
        {
            blockability[currPixelID] = round( 10.0 - 10.0 * (blockability[currPixelID] - minBlocka) / (maxBlocka - minBlocka));
        }

        if (blockability[currPixelID] < 1.0)
            blockability[currPixelID] = 1.0;
	}

	// store
	shape.blockabilityList = blockability; 
}

void saveShape2File(Shape & shape, string folderPath, string filePath, double runningTime)
{   
    if (!fs::exists(folderPath)) 
    {
        if (fs::create_directory(folderPath)) {
            // std::cout << "Create folder for saving puzzle successfully." << std::endl;
        } else {
            std::cerr << "Cannot create folder." << std::endl;
            return;
        }
    } 

    // else 
    // {
    //     std::cout << "Folder existing." << std::endl;

    //     return;
    // }

    std::ofstream outputFile(filePath);

    if (!outputFile.is_open()) {
        std::cerr << "Cannot create a new file." << std::endl;
        return;
    }

    outputFile << shape.volRow << " " << shape.volCol << " " << 1 << "\n" << "\n";

    outputFile << 0.5 << " " << 0.5 << " " << 0.1 << "\n" << "\n";

    std::vector<std::vector<int>> puzzleMatrix = shape.get2DPuzzleMatrix();

    for (const auto &row : puzzleMatrix) {
        for (const int &num : row) {
            outputFile << num << " ";
        }
        outputFile << "\n";  
    }

    outputFile << "\n"; 

    int i = 1;
    for (auto & pieceList : shape.pieceListsByCluster)
    {
        outputFile << "#class_" << i << ": " ;

        for (auto & piece : pieceList)
        {
            outputFile << piece.pieceID << " ";
        }

        outputFile << "\n"; 

        ++i;
    }

    outputFile << "\n"; 

    outputFile << "running time: " << runningTime << " min." << "\n\n";

    outputFile.close();

    std::cout << "Write result into file successfully." << std::endl;
}


Shape readShapeFromFile(string filePath)
{
    set<Pixel> pixels;

    if (filePath.size() < 4)
    {
        std::cerr << "Cannot open file: this is not a .dom input file." << std::endl;
        return Shape(pixels);
    }

    if (filePath.substr(filePath.size() - 4) != ".dom") {
        std::cerr << "Cannot open file: this is not a .dom input file." << std::endl;
        return Shape(pixels);
    }

    std::ifstream inputFile(filePath);

    if (!inputFile.is_open()) {
        std::cerr << "Cannot open file: the file is invalid." << std::endl;
        return Shape(pixels);
    }

    if (filePath.substr(filePath.size() - 4) == ".domain")
    {
        int numRows, numCols, numHeight;
        inputFile >> numRows >> numCols >> numHeight;

        double l, w, h;
        inputFile >> l >> w >> h;

        // 2D matrix to store the coord
        std::vector<std::vector<int>> matrix;

        int row = 0;
        int col = 0;
        int num;

        while (inputFile >> num) {
            
            matrix.push_back({num, row, col});

            col++;

            if (col == numCols) {
                col = 0;
                row++;
            }
        }

        inputFile.close();

        cout << "Initializing region..." << endl;

        // initialize the shape 
        for (const auto &entry : matrix) {
            if (entry[0] == 1) {
                pixels.insert(Pixel(entry[1], entry[2]));
            }
        }

        cout << "Initializing region done." << endl;

        return Shape(pixels, numRows, numCols);
    }

    else{
        int numRows, numCols, numHeight;
        inputFile >> numRows >> numCols >> numHeight;

        double l, w, h;
        inputFile >> l >> w >> h;

        // 2D matrix to store the coord
        std::vector<std::vector<int>> matrix;

        int row = 0;
        int col = 0;
        int num;
        int maxValue = 0;

        while (inputFile >> num) {
            
            matrix.push_back({num, row, col});

            if (num > maxValue)
                maxValue = num;

            col++;

            if (col == numCols) {
                col = 0;
                row++;
            }
        }

        inputFile.close();

        // initialize each piece
        vector<Piece> currPieceList;
        for (int i = 0; i < maxValue; ++i)
        {
            vector<Pixel> currPixelList;

            for (const auto &entry : matrix) {
                if (entry[0] == i + 1) {
                    currPixelList.push_back(Pixel(entry[1], entry[2]));
                }
            }

            currPieceList.push_back(Piece(currPixelList));
        }

        // calculate K
        vector<Piece> templates;
        for (auto & piece : currPieceList)
        {
            bool isExist = false;
            for (auto & curr : templates)
            {
                if (isSameShape(curr, piece))
                {
                    isExist = true;
                    break;
                }
            }

            if (!isExist)
            {
                templates.push_back(piece);
            }
        }

        cout << "The existing tiling solution includes " << templates.size() << " types of tiles\n" << endl;

        Shape currShape = Shape(currPieceList, numRows, numCols);
        currShape.updateTemplates(false);

        currShape.assignColors();
        currShape.assignPieceID();
        currShape.assignPieceClusterID();

        return currShape;
    }
}

Shape readSeedFromFile(string filePath)
{
    set<Pixel> pixels;

    if (filePath.size() < 4)
    {
        std::cerr << "Cannot open file: this is not a .vol or .puz input file." << std::endl;
        return Shape(pixels);
    }

    if (filePath.substr(filePath.size() - 4) != ".vol" && filePath.substr(filePath.size() - 4) != ".puz") {
        std::cerr << "Cannot open file: this is not a .vol or .puz input file." << std::endl;
        return Shape(pixels);
    }

    std::ifstream inputFile(filePath);

    if (!inputFile.is_open()) {
        std::cerr << "Cannot open file: the file is invalid." << std::endl;
        return Shape(pixels);
    }

    int numRows, numCols, numHeight;
    inputFile >> numRows >> numCols >> numHeight;

    double l, w, h;
    inputFile >> l >> w >> h;

    // 2D matrix to store the coord
    std::vector<std::vector<int>> matrix;

    int row = 0;
    int col = 0;
    int num;

    while (inputFile >> num) {
            
        matrix.push_back({num, row, col});

        col++;

        if (col == numCols) {
            col = 0;
            row++;
        }
    }

    inputFile.close();

    // initialize the shape 
    for (const auto &entry : matrix) {
        if (entry[0] >= 1) {
            pixels.insert(Pixel(entry[1], entry[2]));
        }
    }

    Shape currShape = Shape(pixels, numRows, numCols);

    vector<Pixel> seed_pixels;
    vector<Piece> pieces;

    for (const auto &entry : matrix) {
        if (entry[0] > 1) {
            seed_pixels.push_back(Pixel(entry[1], entry[2]));
        }
    }

    for (const auto& px : seed_pixels) 
    {
        pieces.push_back(Piece(px));
    }

    for (auto piece : pieces) 
    {
        currShape.addPiece(piece);
    }

    currShape.assignColors();
    currShape.valid = true;

    Piece newPiece = pieces[0];

    currShape.current_templates.clear();
    currShape.current_templates.push_back(newPiece);

    updateAccessibilityAndBlockability(currShape);

    return currShape;
}


vector<int> GetRandomObjIndexList(vector<float> possibList, float alpha, int objNum)
{
	vector<int> objIndexList;
	if ( objNum > possibList.size() || objNum <= 0 )
	{
		printf("Warning: The input objNum is not correct %d ! \n\n", objNum);
		return objIndexList;
	}

	while( objIndexList.size() < objNum )
	{
		int tempIndex = GetRandomObjIndex(possibList, alpha);
		if ( std::find(objIndexList.begin(), objIndexList.end(), tempIndex) == objIndexList.end() )
		{
			objIndexList.push_back(tempIndex);
		}
	}

	return objIndexList;
}


// Input:  the selection possibility value of each object
// Output: randomly selected object index
int GetRandomObjIndex(vector<float> possibList, float alpha)
{
	if ( possibList.size() == 0)
		return -1;

	if ( possibList.size() == 1)
		return 0;

	vector<float> possibMapList = PossibExpMapping(possibList, alpha);

	// Compute possibility regions for objects with possibility value [P0, P1, P2, ..., P(k-1), Pk]
	// which is [P0, P0+P1, P0+P1+P2, ..., P0+P1+..+P(k-1), P0+P1+..+P(k-1)+Pk]
	vector<float> possibRegions;
	for (int i=0; i<possibMapList.size(); i++)
	{
		float possibSum = 0;
		for (int j=0; j<=i; j++)
		{
			possibSum += possibMapList[j];
		}
		possibRegions.push_back(possibSum);
	}

	// Generate a random value in range [0, P0+P1+..+P(k-1)+Pk)
	int lastObjIndex = possibMapList.size() - 1;
	float randValue = (rand()/(RAND_MAX+1.0)) * possibRegions[lastObjIndex];

	// Return object index by finding which region the random value falls into
	for (int i=0; i<possibRegions.size(); i++)
	{
		// Find each object's possibility range
		float regionMinValue, regionMaxValue;
		if ( i == 0 ) 
		{
			regionMinValue = 0;
			regionMaxValue = possibRegions[0];
		}
		else           
		{
			regionMinValue = possibRegions[i-1]; 
			regionMaxValue = possibRegions[i];
		}

		// Return the randomly selected object index
		if ( randValue >= regionMinValue &&
			 randValue <= regionMaxValue )
		{
			return i;
		}
	}

	printf("Warning: the return random value may not be correct. \n");

	return 0;
}


vector<float> PossibExpMapping(vector<float> possibList, float alpha)
{
	// Check if the input possibility values are correct
	for (int i=0; i<possibList.size(); i++)
	{
		if ( possibList[i] < 0 || possibList[i] > INT_MAX )
		{
			printf("Warning: The input possibility values [i=%2d  P=%.2f] are not correct! \n\n", i, possibList[i]);
		}
	}

	// Possibility value exponential mapping 
	vector<float> possibMapList;
	for (int i=0; i<possibList.size(); i++)
	{
		float tempPossib = pow(possibList[i], alpha);
		possibMapList.push_back(tempPossib);
		//printf("i=%d  P1 %.2f  P2 %.2f \n", i, possibList[i], possibMapList[i]);
	}

	// Calculate the sum of all the possibility values in each list
	//float totalPossib1 = 0;
	//float totalPossib2 = 0;
	//for (int i=0; i<possibList.size(); i++)
	//{
	//	totalPossib1 += possibList[i];
	//	totalPossib2 += possibMapList[i];
	//}
	////printf("Total Possib 1: %.2f \n", totalPossib1);
	////printf("Total Possib 2: %.2f \n", totalPossib2);

	//// Normalize the possibility values in each list
	//for (int i=0; i<possibList.size(); i++)
	//{
	//	possibList[i]    /= totalPossib1;
	//	possibMapList[i] /= totalPossib2;
	//	//printf("i=%d  P1 %.2f  P2 %.2f \n", i, possibList[i], possibMapList[i]);
	//}

	return possibMapList;
}


// Return a float value in [0 1)
float GetRandomNumber(int seedIndex)
{
	if ( seedIndex < 1 )
		printf("Warning: seedIndex cannot be smaller than 1! \n\n");

	float randValue = seedIndex * rand()/(RAND_MAX+1.0);

	if( randValue > 1.0 )
		randValue = randValue - floor(randValue);
	else if( randValue == 1.0 )
		randValue = 0.999999;

	return randValue;
}

void normalizeVector(std::vector<float>& vec) {

    if (vec.size() == 1)
    {
        vec[0] = 1;
        return;
    }
    
    // find min value
    float minValue = *std::min_element(vec.begin(), vec.end());
    float maxValue = *std::max_element(vec.begin(), vec.end());

    if (minValue == maxValue)
    {
        for (float& value : vec) 
            value = 1;
    
        return;
    }

    for (float& value : vec) {
        value = (value - minValue) / float(maxValue - minValue);
    }
}

// Function to generate a random unsaturated color
// Function to generate a random unsaturated color
void generateUnsaturatedColor(Color_eigen& color, double minSaturation, double maxSaturation, double minValue, double maxValue, const std::vector<Color_eigen>& existingColors) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, 359);

    // Generate a random hue value (0-360) for an unsaturated color
    int h = distribution(gen);

    // Generate random saturation and value within the specified range
    std::uniform_real_distribution<double> satDistribution(minSaturation, maxSaturation);
    std::uniform_real_distribution<double> valDistribution(minValue, maxValue);
    double s = satDistribution(gen);
    double v = valDistribution(gen);

    // Convert HSV to RGB
    double c = s * v;
    double x = c * (1 - std::abs(std::fmod(h / 60.0, 2) - 1));
    double m = v - c;

    double r1, g1, b1;
    if (h < 60) {
        r1 = c;
        g1 = x;
        b1 = 0;
    } else if (h < 120) {
        r1 = x;
        g1 = c;
        b1 = 0;
    } else if (h < 180) {
        r1 = 0;
        g1 = c;
        b1 = x;
    } else if (h < 240) {
        r1 = 0;
        g1 = x;
        b1 = c;
    } else if (h < 300) {
        r1 = x;
        g1 = 0;
        b1 = c;
    } else {
        r1 = c;
        g1 = 0;
        b1 = x;
    }

    // Scale to 0-255 range and set color values
    color << (r1 + m) * 255, (g1 + m) * 255, (b1 + m) * 255;

    // Calculate the squared Euclidean distance between the new color and existing colors
    double minDistanceSquared = std::numeric_limits<double>::max();
    for (const auto& existingColor : existingColors) {
        double dr = color[0] - existingColor[0];
        double dg = color[1] - existingColor[1];
        double db = color[2] - existingColor[2];
        double distSquared = dr * dr + dg * dg + db * db;
        minDistanceSquared = std::min(minDistanceSquared, distSquared);
    }

    // If the minimum distance is below a threshold, regenerate the color
    if (minDistanceSquared < 300) { // Adjust threshold as needed
        generateUnsaturatedColor(color, minSaturation, maxSaturation, minValue, maxValue, existingColors);
    }
}