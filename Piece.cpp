#include "Piece.h"
#include <iostream>

// Function to find the index of a row in an Eigen matrix
int findRow(const Eigen::MatrixXd& matrix, const Eigen::RowVectorXd& targetRow) {
    for (int i = 0; i < matrix.rows(); ++i) {
        if (matrix.row(i) == targetRow) {
            return i;
        }
    }
    return -1;  // Return -1 if the row is not found
}

void Pixel::print()
{
    cout << "Pixel( " << x << ", " << y << ")" << endl;
}

Piece::Piece() {

}


//constructs a piece from a single pixel
Piece::Piece(const Pixel& pixel) {
    this->pixels = vector<Pixel>{ pixel };
    initializeMatrix();
    this->color = Color{ 255,255,255 };
}

Piece::Piece(const vector<Pixel> &pixels,const Color &color) {
	this->pixels = pixels;
    initializeMatrix();
	this->color = color;
}

Piece::Piece(const vector<Pixel> &pixels,const Color &color, const int & pieceID) {
	this->pixels = pixels;
    initializeMatrix();
	this->color = color;
    this->pieceID = pieceID;
}

Piece::Piece(const Eigen::MatrixXi& pixelsAsMatrix, const Color& color) {
    this->pixelsAsMatrix = pixelsAsMatrix;
    for (int i = 0; i < pixelsAsMatrix.rows(); i++) {
        this->pixels.push_back(Pixel(pixelsAsMatrix.coeff(i, 0), pixelsAsMatrix.coeff(i, 1)));
    }
    this->color = color;
}

Piece::Piece(const vector<Pixel> &pixels)
{
    this->pixels = pixels;
    initializeMatrix(); 
    this->color = Color{ 255,255,255 };
}

void Piece::initializeMatrix() {
    pixelsAsMatrix = Eigen::MatrixXi(pixels.size(), 2);
    int row = 0;
    for (const auto& pixel : pixels) {
        pixelsAsMatrix.row(row) << pixel.x, pixel.y;
        row++;
    }
}

void Piece::initializePieceShapeAsMatrix()
{
    // 找到最小和最大的行和列坐标
    int minRow = INT_MAX, maxRow = INT_MIN;
    int minCol = INT_MAX, maxCol = INT_MIN;

    for (const auto & pixel : pixels) {
        minRow = min(minRow, pixel.x);
        maxRow = max(maxRow, pixel.x);
        minCol = min(minCol, pixel.y);
        maxCol = max(maxCol, pixel.y);
    }

    // 计算裁剪后的行数和列数
    int numRows = maxRow - minRow + 1;
    int numCols = maxCol - minCol + 1;

    // 创建二维数组并初始化为0
    vector<vector<int>> croppedPiece(numRows, vector<int>(numCols, 0));

    // 将像素坐标映射到二维数组中
    for (const Pixel& pixel : pixels) {
        int newRow = pixel.x - minRow;
        int newCol = pixel.y - minCol;
        croppedPiece[newRow][newCol] = 1; // 1表示该位置有像素
    }

    pieceShapeAsMatrix = croppedPiece;
}

void Piece::generateRenderMeshes() {
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
        V.row(Vrow) << pixel.x - 0.5, pixel.y - 0.5;
        V1.row(Virow) << pixel.x - 0.5, pixel.y - 0.5;

        V.row(Vrow + 1) << pixel.x - 0.5, pixel.y + 0.5;
        V2.row(Virow) << pixel.x - 0.5, pixel.y + 0.5;

        V.row(Vrow + 2) << pixel.x + 0.5, pixel.y + 0.5;
        V3.row(Virow) << pixel.x + 0.5, pixel.y + 0.5;

        V.row(Vrow + 3) << pixel.x + 0.5, pixel.y - 0.5;
        V4.row(Virow) << pixel.x + 0.5, pixel.y - 0.5;

        F.row(Frow) << Vrow, Vrow + 1, Vrow + 2;
        F.row(Frow + 1) << Vrow, Vrow + 3, Vrow + 2;

        Vrow += 4;
        Virow++;
        Frow += 2;
            
    }

    int boundaryEdgeNum = 0;
    vector<Pixel> currPixels = pixels;
    vector< pair<Eigen::Vector2d, Eigen::Vector2d>  > edgeEndPointsList;

    for (int k = 0; k < pixels.size(); ++k)
    {
        Pixel currPixel = pixels[k];

        vector<Pixel> possibleNeigbours;
        possibleNeigbours.push_back(Pixel(currPixel.x + 1,currPixel.y)); // right
        possibleNeigbours.push_back(Pixel(currPixel.x,currPixel.y - 1)); // down
        possibleNeigbours.push_back(Pixel(currPixel.x - 1,currPixel.y)); // left
        possibleNeigbours.push_back(Pixel(currPixel.x,currPixel.y + 1)); // top

        int m = 0;
        for (auto neighbour : possibleNeigbours)
        {
            bool isInThisPiece = false;

            for (int i = 0; i < currPixels.size(); ++i)
            {
                if (neighbour == currPixels[i])
                {
                    isInThisPiece = true;
                    break;
                }
            }
            if (!isInThisPiece)
            {
                // cout << "find a bounadry edge. " << endl;
                boundaryEdgeNum++;

                // right
                if (m == 0)
                    edgeEndPointsList.push_back(pair<Eigen::Vector2d, Eigen::Vector2d>(Eigen::Vector2d(currPixel.x + 0.5, currPixel.y + 0.5), Eigen::Vector2d(currPixel.x + 0.5, currPixel.y - 0.5)));

                // down
                if (m == 1)
                    edgeEndPointsList.push_back(pair<Eigen::Vector2d, Eigen::Vector2d>(Eigen::Vector2d(currPixel.x + 0.5, currPixel.y - 0.5), Eigen::Vector2d(currPixel.x - 0.5, currPixel.y - 0.5)));
                
                // left
                if (m == 2)
                    edgeEndPointsList.push_back(pair<Eigen::Vector2d, Eigen::Vector2d>(Eigen::Vector2d(currPixel.x - 0.5, currPixel.y - 0.5), Eigen::Vector2d(currPixel.x - 0.5, currPixel.y + 0.5)));

                // top
                if (m == 3)
                    edgeEndPointsList.push_back(pair<Eigen::Vector2d, Eigen::Vector2d>(Eigen::Vector2d(currPixel.x + 0.5, currPixel.y + 0.5), Eigen::Vector2d(currPixel.x - 0.5, currPixel.y + 0.5))); 
            }

            m++;
        }
    }

    V_boundary = Eigen::MatrixXd(boundaryEdgeNum * 4, 3);
    F_boundary = Eigen::MatrixXi(boundaryEdgeNum * 2, 3);

    Vrow = 0;
    Frow = 0;
    double shortSide = 0.05;

    for (int i = 0; i < edgeEndPointsList.size(); ++i)
    {
        Eigen::Vector2d topLeft, topRight, bottomLeft, bottomRight;

        Eigen::Vector2d p1 = edgeEndPointsList[i].first;
        Eigen::Vector2d p2 = edgeEndPointsList[i].second;

        p1 = p1 + 0.5 * shortSide * (p1 - p2);
        p2 = p2 + 0.5 * shortSide * (p2 - p1);
        
        // 计算矩形的两个边向量
        Eigen::Vector2d midPoint = (p1 + p2) / 2.0;
        Eigen::Vector2d e1 = (p2 - p1) / 2.0;
        Eigen::Vector2d e2(-e1.y(), e1.x());  // 旋转90度，得到垂直边向量

        e2.normalize();
        e2 = (0.5 * shortSide) * e2;

        // 计算矩形的四个顶点
        topLeft = midPoint + e2 - e1;
        topRight = midPoint + e2 + e1;
        bottomLeft = midPoint - e2 - e1;
        bottomRight = midPoint - e2 + e1;

        V_boundary.row(Vrow)     << topLeft(0), topLeft(1), 0.1;
        V_boundary.row(Vrow + 1) << bottomLeft(0), bottomLeft(1), 0.1;
        V_boundary.row(Vrow + 2) << bottomRight(0), bottomRight(1), 0.1;
        V_boundary.row(Vrow + 3) << topRight(0), topRight(1), 0.1;

        F_boundary.row(Frow) << Vrow, Vrow + 1, Vrow + 2;
        F_boundary.row(Frow + 1) << Vrow, Vrow + 3, Vrow + 2;

        Vrow += 4;
        Frow += 2;
    }

    // color matrix
    C = Eigen::MatrixXd(pixels.size() * 2, 3);
    for (int i = 0; i < pixels.size() * 2; i++) {
        C.row(i) << color.r, color.g, color.b;
    }
}

void Piece::generateSavingMeshes()
{
    V_obj = Eigen::MatrixXd(pixels.size() * 4, 3);
    F_obj = Eigen::MatrixXi(pixels.size() * 2, 3);

    int Vrow = 0;
    int Virow = 0;
    int Frow = 0;

    for (auto pixel : pixels) {
        V_obj.row(Vrow) << pixel.x - 0.5, pixel.y - 0.5, 0;

        V_obj.row(Vrow + 1) << pixel.x - 0.5, pixel.y + 0.5, 0;

        V_obj.row(Vrow + 2) << pixel.x + 0.5, pixel.y + 0.5, 0;

        V_obj.row(Vrow + 3) << pixel.x + 0.5, pixel.y - 0.5, 0;

        F_obj.row(Frow) << Vrow, Vrow + 2, Vrow + 1;
        F_obj.row(Frow + 1) << Vrow, Vrow + 3, Vrow + 2;

        Vrow += 4;
        Virow++;
        Frow += 2;
    }

    // Convert Eigen matrix to std::vector<std::vector<double>>
    std::vector<std::vector<double>> matrixVec(V_obj.rows(), std::vector<double>(V_obj.cols()));
    for (int i = 0; i < V_obj.rows(); ++i) {
        for (int j = 0; j < V_obj.cols(); ++j) {
            matrixVec[i][j] = V_obj(i, j);
        }
    }

    // Remove duplicate rows from the std::vector
    std::sort(matrixVec.begin(), matrixVec.end());
    auto last = std::unique(matrixVec.begin(), matrixVec.end());
    matrixVec.erase(last, matrixVec.end());

    // Convert std::vector<std::vector<double>> back to Eigen matrix
    Eigen::MatrixXd V_obj_result(matrixVec.size(), matrixVec[0].size());
    for (int i = 0; i < V_obj_result.rows(); ++i) {
        for (int j = 0; j < V_obj_result.cols(); ++j) {
            V_obj_result(i, j) = matrixVec[i][j];
        }
    }

    // cout << endl << V_obj_result << endl << endl;

    Eigen::MatrixXi F_obj_result = Eigen::MatrixXi(pixels.size() * 2, 3);

    for (int i = 0; i < F_obj.rows(); ++i)
    {
        F_obj_result.row(i) << findRow(V_obj_result, V_obj.row(F_obj(i, 0))), findRow(V_obj_result, V_obj.row(F_obj(i, 1))), findRow(V_obj_result, V_obj.row(F_obj(i, 2)));
    }

    // cout << endl << F_obj_result << endl << endl;

    V_obj = V_obj_result;
    F_obj = F_obj_result;
}

bool Piece::hasPixel(const Pixel& pixel) const {
    for (const auto& px : pixels) {
        if (pixel == px) {
            return true;
        }
    }
    return false;
}

void Piece::setColor(int r, int g, int b) {
    this->color = Color{ r,g,b };
}

void Piece::setColor(const Color &c) {
    this->color = c;
}

void Piece::addPixel(const Pixel& px) {
    pixels.push_back(px);
    initializeMatrix();
}

void Piece::deletePixel(const Pixel & px)
{
    for (int i = 0; i < pixels.size(); ++i)
    {
        if (px == pixels[i])
        {
            pixels.erase(pixels.begin() + i);
        }
    }
}


// returns a rotated piece given a rotation matrix and the pixel to rotate about
Piece Piece::rotate(const Eigen::MatrixXi& rotation_matrix, int px_index) const {
    Eigen::Vector2i translation;
    translation << -pixels[px_index].x, -pixels[px_index].y;
    
    // translate to origin
    Eigen::MatrixXi translated_pixels = pixelsAsMatrix.rowwise() + translation.transpose();
    // rotate
    Eigen::MatrixXi rotated_pixels = (rotation_matrix * translated_pixels.transpose()).transpose();
    // translate back to original position
    Eigen::Vector2i translate_back;
    translate_back << pixels[px_index].x, pixels[px_index].y;
    rotated_pixels = rotated_pixels.rowwise() + translate_back.transpose();
    
    Piece newPiece = Piece(rotated_pixels, this->color);
    
    return newPiece;

}

vector<Piece> Piece::get_all_possible_flips_and_rotations() const {
    vector<Eigen::MatrixXi> all_rotations;

    Eigen::MatrixXi zero = Eigen::MatrixXi(2, 2);
    zero << 1, 0,
        0, 1;

    Eigen::MatrixXi ninety = Eigen::MatrixXi(2, 2);
    ninety << 0, -1,
        1, 0;

    Eigen::MatrixXi one_eighty = Eigen::MatrixXi(2, 2);
    one_eighty << -1, 0,
        0, -1;

    Eigen::MatrixXi two_seventy = Eigen::MatrixXi(2, 2);
    two_seventy << 0, 1,
        -1, 0;

    all_rotations.push_back(zero);
    all_rotations.push_back(ninety);
    all_rotations.push_back(one_eighty);
    all_rotations.push_back(two_seventy);

    vector<Eigen::MatrixXi> all_flips;

    Eigen::MatrixXi no_flip = Eigen::MatrixXi(2, 2);
    no_flip << 1, 1,
        1, 1;


    Eigen::MatrixXi vertical_flip = Eigen::MatrixXi(2, 2);
    vertical_flip << -1, 0,
        0, 1;

    Eigen::MatrixXi horizontal_flip = Eigen::MatrixXi(2, 2);
    horizontal_flip << 1, 0,
        0, -1;

    all_flips.push_back(no_flip);
    all_flips.push_back(vertical_flip);
    all_flips.push_back(horizontal_flip);

    vector<Piece> all_possible_pieces;

    for (int i = 0; i < pixels.size();i++) {
        for (const auto& rotation : all_rotations) {
            for (const auto& flip : all_flips) {
                Piece flipped = this->rotate(flip, i);
                Piece rotated_and_flipped = flipped.rotate(rotation, i);
                all_possible_pieces.push_back(rotated_and_flipped);
            }
        }
    }

    return all_possible_pieces;
}

bool Piece::isSamePiece(const Piece& other) {
    set<Pixel> this_pixels(this->pixels.begin(), this->pixels.end());
    set<Pixel> other_pixels(other.pixels.begin(), other.pixels.end());
    return this_pixels == other_pixels;
}


string Piece::printPixels() {
    string out = "";
    out.append("Color : { " + to_string(color.r) + ',' + to_string(color.g) + ',' + to_string(color.b) + string(" }\n"));
    for (auto px : pixels) {
        out.append("Pixel: " + to_string(px.x) + ',' + to_string(px.y) + '\n');
    }

    return out;
}

void Piece::printPieceMatrix()
{
    initializePieceShapeAsMatrix();

    // Get the number of rows and columns in the matrix
    int numRows = pieceShapeAsMatrix.size();
    int numCols = pieceShapeAsMatrix[0].size(); // Assuming all rows have the same number of columns

    // Print the 2D matrix
    // std::cout << "Piece:" << std::endl;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << pieceShapeAsMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}