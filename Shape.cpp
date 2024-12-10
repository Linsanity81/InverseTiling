#include "Shape.h"
#include "Util.h"
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/closeness_centrality.hpp>
#include <boost/graph/breadth_first_search.hpp>

// Define an undirected graph using the Boost Graph Library
struct vertex_info { 
    int x; 
    int y; 
    
	float betweenessCentrality; 
	float closenessCentrality;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_info> Graph;

// Function to export Eigen matrices V and F to an OBJ file
void exportToObj(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::string& outputPath) {
    std::ofstream objFile(outputPath);
    if (!objFile.is_open()) {
        std::cerr << "Error opening file: " << outputPath << std::endl;
        return;
    }

    // Write vertices to the OBJ file
    for (int i = 0; i < V.rows(); ++i) {
        objFile << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
    }

    // Write faces to the OBJ file
    for (int i = 0; i < F.rows(); ++i) {
        objFile << "f";
        for (int j = 0; j < F.cols(); ++j) {
            objFile << " " << F(i, j) + 1; // OBJ indices start from 1
        }
        objFile << std::endl;
    }

    objFile.close();
    // std::cout << "OBJ file saved at: " << outputPath << std::endl;
}

Shape::Shape() {

}

Shape::Shape(const set<Pixel>& all_pixels) {
	this->all_pixels = all_pixels;
	this->backtrack_count = 0;
}

Shape::Shape(const set<Pixel> &all_pixels, const int row, const int col) {
	this->all_pixels = all_pixels;
	this->backtrack_count = 0;

	this->volRow = row;
	this->volCol = col;
}

Shape::Shape(const vector<Piece>& pieces) {
	this->backtrack_count = 0;
	this->pieces = pieces;

	for (const auto& piece : pieces) {
		for (auto px : piece.pixels) {
			all_pixels.insert(px);
		}
	}
}

Shape::Shape(const vector<Piece>& pieces, const int row, const int col)
{
	this->backtrack_count = 0;
	this->pieces = pieces;

	for (const auto& piece : pieces) {
		for (auto px : piece.pixels) {
			all_pixels.insert(px);
		}
	}

	this->volRow = row;
	this->volCol = col;
}

void Shape::clear()
{
	pieces.clear();

	current_templates.clear();
	blocked_templates.clear();

	pieceListsByCluster.clear();

	growthLog = "";

	backtrack_count = 0;
}

void Shape::addPiece(const Piece& piece) {
	pieces.push_back(piece);
	for (auto px : piece.pixels) {
		all_pixels.insert(px);
	}
}

bool Shape::hasPixel(const Pixel& pixel) const {
	return all_pixels.count(pixel);
}

bool Shape::pixelInPieces(const Pixel & pixel)
{
	for (auto piece : pieces)
	{
		if (piece.hasPixel(pixel))
			return true;
	}

	return false;
}

// assigns a unique color to each piece in the shape
void Shape::assignColors() {
	vector<Color> colors = { Color{150,150,150},
							Color{255,0,0}, 
	                     Color{0,255,0}, 
                         Color{0,0,255}, 
                         Color{255, 255, 0}, 
                         Color{255,0,255}, 
                         Color{255,102,255},
	                     Color{0,153,153}, 
                         Color{255,128,0}, 
                         Color{20,20,20},
						 Color{236, 82, 21},
						 Color{203,13,120}
						  };

	for (int i = 0; i < pieces.size(); i++) {
		pieces[i].setColor(colors[i + 1]);
	}
}



// returns the current size of each piece in the shape
// returns -1 if the pieces are of different sizes
int Shape::getAllPiecesSize() {
	if (pieces.size() == 0) {
		return -1;
	}
	int firstPieceSize = pieces[0].pixels.size();

	for (int i = 1; i < pieces.size(); i++) {
		if (pieces[i].pixels.size() != firstPieceSize) {
			return -1;
		}
	}

	return firstPieceSize;
}

vector<Pixel> Shape::getUnfilledPixels() {
	vector<Pixel> unfilled_px;

	for (const auto& px : all_pixels) {
		if (!isPixelInPieces(px, pieces)) {
			unfilled_px.push_back(px);
		}
	}

	return unfilled_px;
}

int Shape::getNumUnfilledPx() {
	int filled = 0;
	for (const auto& piece : pieces) {
		filled += piece.pixels.size();
	}

	return all_pixels.size() - filled;
}

void Shape::addLogGrowPiece(int piece_idx, Color piece_color) {
	// growthLog.append("\n");
	// growthLog.append("Grew piece: " + piece_idx + string(" Color = (") + to_string(int(piece_color.r)) + "," + to_string(int(piece_color.g)) + "," + to_string(int(piece_color.b)) + ")\n");
	// growthLog.append("\n");
	growthLog.append(pieces[piece_idx].printPixels());
}

void Shape::addLogCurrentShapeState() {
	growthLog.append("\nCurrent Shape: \n");

	if (!valid) {
		growthLog.append("Shape is now Invalid.\n");
	}

	growthLog.append("Current pieces: \n");
	for (int i = 0; i < pieces.size(); i++) {
		growthLog.append("\nPiece " + to_string(i) + string(":\n"));
		growthLog.append(pieces[i].printPixels());
	}

	growthLog.append("\nUnfilled pixels: \n");
	for (auto px : getUnfilledPixels()) {
		growthLog.append("Pixel: " + to_string(px.x) + ',' + to_string(px.y) + '\n');
	}
}

void Shape::printTemplates()
{
	cout << "-------------------------------------" << endl << endl;

	for (int i = 0; i < current_templates.size(); ++i)
	{
        cout << "current Template " << i << ": " << endl;
        current_templates[i].printPieceMatrix();
        cout << endl;
    }

    for (int i = 0; i < blocked_templates.size(); ++i)
    {
        cout << "blocked Template " << i << ": " << endl;
        blocked_templates[i].printPieceMatrix();
        cout << endl;
    }

	cout << "-------------------------------------" << endl << endl;
}

string Shape::printAllPixels()
{
    string out = "";
    for (auto px : all_pixels) {
        out.append("Pixel: " + to_string(px.x) + ',' + to_string(px.y) + '\n');
    }

    return out;
}

void Shape::assignPieceID()
{
	int currentID = 1;
	for (auto & piece : pieces)
	{
		piece.pieceID = currentID;
		currentID++;
	}
}

void Shape::assignPieceClusterID()
{
	pieceListsByCluster.clear();

	for (auto & piece : current_templates)
	{
		vector<Piece> currPieceList;
		currPieceList.push_back(piece);
		pieceListsByCluster.push_back(currPieceList);
	}

	for (auto & piece : blocked_templates)
	{
		vector<Piece> currPieceList;
		currPieceList.push_back(piece);
		pieceListsByCluster.push_back(currPieceList);
	}

	for (auto & piece : pieces)
	{
		for (auto & pieceList : pieceListsByCluster)
		{
			if (isSameShape(piece, pieceList[0]))
			{
				pieceList.push_back(piece);
				// std::cout << "Push_back successfully." << std::endl;
				break;
			}
		}
	}

	for (auto & pieceList : pieceListsByCluster)
	{
		pieceList.erase(pieceList.begin());
	}

	// 使用 lambda 函数作为排序的比较函数
    std::sort(pieceListsByCluster.begin(), pieceListsByCluster.end(),
              [](const std::vector<Piece>& a, const std::vector<Piece>& b) {
                  return a.size() > b.size();  // 按元素数量降序排序
              });

	// assign clusterID to each piece
	for (auto & piece : pieces)
	{
		int i = 1;
		for (auto & pieceList : pieceListsByCluster)
		{
			if (isSameShape(piece, pieceList[0]))
			{
				piece.pieceClusterID = i;
				break;
			}
			i++;
		}
	}
}

std::vector<std::vector<int>> Shape::get2DPuzzleMatrix()
{
	assignPieceID();
	assignPieceClusterID();

	std::vector<std::vector<int>> matrix(volRow, std::vector<int>(volCol, 0));

	for (const auto & piece : pieces)
	{
		int currentPieceID = piece.pieceID;
		for (const auto & pixel : piece.pixels)
		{
			matrix[pixel.x][pixel.y] = currentPieceID;
		}
	}

	return matrix;
}

string Shape::getOutputFolderPath(string prefix)
{
	string extension = "_N" + to_string(pieces.size()) + "_K" + to_string(score);

	// cout << prefix + extension << endl;

	return prefix + extension;
}

string Shape::getOutputFileFullName(string prefix)
{
	string extension = "_N" + to_string(pieces.size()) + "_K" + to_string(score);

	// cout << prefix + extension + "/" + extension + ".puz" << endl;

	return prefix + extension + ".puz";
}

int Shape::getSmallestPieceSize()
{
	int minSize = INT_MAX;

	for (auto & piece : pieces)
	{
		if (piece.pixels.size() < minSize)
		{
			minSize = piece.pixels.size();
		}
	}

	return minSize;
}

int Shape::getLargestPieceSize()
{
	int maxSize = 0;

	for (auto & piece : pieces)
	{
		if (piece.pixels.size() > maxSize)
		{
			maxSize = piece.pixels.size();
		}
	}

	return maxSize;
}

double Shape::getInstanceNumDeviation()
{
	assignPieceClusterID();

	vector<int> instanceNumList(current_templates.size(), 0);

	for (auto & piece : pieces)
	{
		instanceNumList[piece.pieceClusterID - 1] ++;
	}

	cout << "intanceNumList: " << endl;
	for (auto & i : instanceNumList)
	{
		cout << i << " ";
	}
	cout << endl;

	double averageInstanceNum = double(pieces.size()) / double(current_templates.size());
	double deviation = 0;

	for (auto & i : instanceNumList)
	{
		deviation += (fabs(averageInstanceNum - double(i)) * fabs(averageInstanceNum - double(i)));
	}

	cout << "currentInstanceNumDeviation: " << deviation << endl << endl;

	return deviation;
}

bool Shape::isSatisfyMinSizeRequirement()
{
	double target_template_size = ceil((double)all_pixels.size() / (double)pieces.size());

	// cout << "smallest size: " << getSmallestPieceSize() << endl;

	// if (getSmallestPieceSize() >= ceil(target_template_size * minAvgSizeRatio))
	if (getSmallestPieceSize() >= 4)
		return true;
	
	else
		return false;
}

int Shape::getPixelID(const Pixel & pixel)
{
	return (pixel.x * volCol + pixel.y);
}

int Shape::getPieceID(const Piece & piece)
{
	for (int i = 0; i < pieces.size(); ++i)
	{
		if (pieces[i].isSamePiece(piece))
			return i;
	}

	return -1;
}

int Shape::getSmallSizePiecesNum(int currMinTileSize)
{
	int count = 0;
	for (int i = 0; i < pieces.size(); ++i)
	{
		if (pieces[i].pixels.size() < currMinTileSize)
		{
			count++;	
		}
	}

	return count;
}

int Shape::getLargeSizePiecesNum(int currMaxTileSize)
{
	int count = 0;
	for (int i = 0; i < pieces.size(); ++i)
	{
		if (pieces[i].pixels.size() > currMaxTileSize)
		{
			count++;	
		}
	}

	return count;
}

double Shape::getAverageEnlargeability()
{
	float count = 0;
	for (auto & piece : pieces)
	{
		vector<Pixel> neighbourPixels = getNeighbouringPixels(piece, *this);

		for (auto & neighourPixel : neighbourPixels)
        {
            int currPixelID = getPixelID(neighourPixel);

            // count += round(0.6 * shape.accessibilityList[currPixelID] + 0.4 * shape.blockabilityList[currPixelID]);
            count += round(blockabilityList[currPixelID]);
        }
	}

	return count / double(pieces.size());
}

void Shape::updateTemplates(bool isLooseSizeRequirement)
{
    current_templates.clear();
    blocked_templates.clear();

    for (int i = 0; i < pieces.size(); ++i)
    {
        pieces[i].isBlocked = isBlocked(pieces[i], *this);

        if (!isPieceinTemplates(pieces[i], current_templates))
            current_templates.push_back(pieces[i]);
    }

	int currMaxSize = maxTileSize;
    int currMinSize = minTileSize;

	if (getNumUnfilledPx() <= 10 and isLooseSizeRequirement)
    {
        currMaxSize += 1;
        currMinSize -= 1;
    }

    for (auto & tmp : current_templates)
    {
        bool isBlockedTmp = true;
        vector<Piece> instances = getTemplateInstances(tmp, *this);

        for (auto & instance : instances)
        {
            if (!instance.isBlocked)
            {
                isBlockedTmp = false;
                break;
            }
        }

        if (isBlockedTmp or tmp.pixels.size() >= currMaxSize)
            blocked_templates.push_back(tmp);
    }

	score = current_templates.size();
}

void Shape::updateCentrality()
{
	// 0）build the graph and compute accessibility
	Graph g;
	for (int i = 0; i < volRow; ++i)
	{
		for (int j = 0; j < volCol; ++j)
		{
			boost::add_vertex(g);
		}
	}

	std::vector<float> accessibility(num_vertices(g), 0);

	int i = 0;
	vector<Pixel> unfilledPixels = getUnfilledPixels();
	for (auto & pixel : unfilledPixels)
	{
		int distinctPieceNum;
		vector<Pixel> unfilledNeigbourPixels = getNeighbouringPixels(pixel, *this, distinctPieceNum);
		int currPixelID = getPixelID(pixel);

		accessibility[currPixelID] = distinctPieceNum;

		for (auto & px : unfilledNeigbourPixels)
		{
			int currNeigbourPixelID = getPixelID(px);
			boost::add_edge(currPixelID, currNeigbourPixelID, g);
		}

 		++i;
	}

	normalizeVector(accessibility);

	// 1) compute closeness centrality
	std::vector<float> closenessCentrality(num_vertices(g), 0);
	for (auto & pixel : unfilledPixels)
	{
		int currPixelID = getPixelID(pixel);
		
		// vector to store distances for each vertex
    	std::vector<int> distances(boost::num_vertices(g));

		// bfs
    	boost::breadth_first_search(g, currPixelID, boost::visitor(boost::make_bfs_visitor(boost::record_distances(distances.data(), boost::on_tree_edge()))));

		closenessCentrality[currPixelID] = 0;

		for (auto & px : unfilledPixels)
		{
			int _currPixelID = getPixelID(px);
			closenessCentrality[currPixelID] += distances[_currPixelID];
		}

		// closenessCentrality[currPixelID] = float(all_pixels.size() - 1.0) / closenessCentrality[currPixelID];

		if (closenessCentrality[currPixelID] != 0)
			closenessCentrality[currPixelID] = ( 1.0 / closenessCentrality[currPixelID] ) / float(unfilledPixels.size() - 1);
	}

	normalizeVector(closenessCentrality);

	// 2) compute betweeness centrality
	std::vector<float> betweenessCentrality(num_vertices(g), 0);
	boost::brandes_betweenness_centrality(g, boost::make_iterator_property_map(betweenessCentrality.begin(), get(boost::vertex_index, g)));

	normalizeVector(betweenessCentrality);

	// 3) store the centrality to the pixel list
	accessibilityList = accessibility;
	closenessCentralityList = closenessCentrality;
	betweenessCentralityList = betweenessCentrality;

	remainingGrowingEvaluation = 0;

	for (auto & pixel : all_pixels)
	{
		Pixel currPixel = pixel;

		if (isPixelInPixelList(currPixel, unfilledPixels))
		{
			int currPixelID = getPixelID(pixel);

			remainingGrowingEvaluation += (round (2.0 * (1.0 - closenessCentralityList[currPixelID]) ) 
                                    + round (2.0 * (1.0 - betweenessCentralityList[currPixelID]) ) 
                                    + round (1.0 * (1.0 - accessibilityList[currPixelID]) ) );

			// cout << 1.0 - closenessCentrality[currPixelID] << " " << 1.0 - betweenessCentrality[currPixelID] << " " << 1.0 - accessibility[currPixelID] << endl;
		}
	}

	// cout << endl;

	// for (auto & pixel : all_pixels)
	// {
	// 	Pixel currPixel = pixel;

	// 	currPixel.print();

	// 	cout << currPixel.closenessCentrality << " " << currPixel.betweenessCentrality << endl;
	// }

	// print for debug
	// boost::graph_traits<Graph>::vertex_iterator vi, vend;
	// std::cout << "Betweenness Centrality:\n";
    // for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	// {
	// 	std::cout << "Vertex " << *vi << ": " << betweenessCentrality[*vi] << std::endl;
    //     // std::cout << "Vertex " << *vi << ": (" << g[*vi].x << ", " << g[*vi].y << ") " << betweenessCentrality[*vi] << std::endl;
    // }

	// std::cout << "Closeness Centrality:\n";
    // for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi) 
	// {
	// 	std::cout << "Vertex " << *vi << ": " << closenessCentrality[*vi] << std::endl;
    //     // std::cout << "Vertex " << *vi << ": (" << g[*vi].x << ", " << g[*vi].y << ") " << closenessCentrality[*vi] << std::endl;
    // }
}

void Shape::saveEachPieceOBJs(string prefix)
{
	int i = 1;
	for (auto & piece : pieces)
	{
		piece.generateSavingMeshes();

		exportToObj(piece.V_obj, piece.F_obj, prefix + "/piece" + to_string(i) + "_class" + to_string(piece.pieceClusterID) + ".obj");

		++i;
	}
}

void Shape::saveEachClassTemplate(string prefix)
{
	if (!fs::exists(prefix + "/class_templates")) {
        if (fs::create_directory(prefix + "/class_templates")) {
            std::cout << "Create folder for saving puzzle class template successfully." << std::endl;
        } else {
            std::cerr << "Cannot create folder." << std::endl;
        }
    } else {
        std::cout << "Folder existing." << std::endl;

        return;
    }

	int i = 1;

	for (auto & pieceList : pieceListsByCluster)
	{
		// cout << "pieceList size: " << pieceList.size() << endl;

		for (auto & piece : pieceList)
		{
			piece.generateSavingMeshes();

			exportToObj(piece.V_obj, piece.F_obj, prefix + "/class_templates/class_template_" + to_string(i) + ".obj");

			break;
			
		}
		++i;
	}
}