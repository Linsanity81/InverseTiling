#include "Viewer.h"
#include "Util.h"
#include <iostream>
#include <igl/stb/write_image.h>
#include <igl/stb/read_image.h>
#include <random>

// Color scheme table
Eigen::RowVector3d colorTable[32] = {

		// Eigen::RowVector3d(0.85, 0.85, 0.85),	  //  1: light green
		Eigen::RowVector3d(0.63, 0.9, 0.63),	  //  1: light green 
		Eigen::RowVector3d(0.95, 0.95, 0.66),     //  2: light yellow
		Eigen::RowVector3d(0.81, 0.63, 0.9),   	  //  3: light purple
        Eigen::RowVector3d(0.9, 0.74, 0.63),      //  4: light brown
        Eigen::RowVector3d(0.63, 0.63, 0.9),      //  5: light blue

		Eigen::RowVector3d(0.64, 0.9, 0.9),       //  6: light cyan
		Eigen::RowVector3d(0.9, 0.63, 0.63),      //  7: light red
		Eigen::RowVector3d(0.9, 0.63, 0.8),       //  8: light pink
		Eigen::RowVector3d(0.8, 1.0, 0.7),        //  9: light pink
		Eigen::RowVector3d(0.7, 0.85, 1.0),       //  10: light pink

		Eigen::RowVector3d(0.392, 0.784, 0.706),  //  11: light pink
		Eigen::RowVector3d(0.8, 0.8, 0.56),       //  12: light pink
		Eigen::RowVector3d(0.95, 0.78, 0.67),     //  13: light pink
		Eigen::RowVector3d(0.582, 0.701, 0.9),    //  14: light pink
		Eigen::RowVector3d(0.714, 0.522, 0.8),    //  15: Red

		Eigen::RowVector3d(0.627, 0.949, 0.102),  //  16: Red
		Eigen::RowVector3d(0.927, 0.588, 0.948),  //  17: Red
		Eigen::RowVector3d(0.424, 0.965, 0.796),  //  18: light pink
		Eigen::RowVector3d(0.902, 0.510, 0.745),  //  19: light pink
		Eigen::RowVector3d(0.941, 0.941, 0.475),  //  20: Red

		Eigen::RowVector3d(0.900, 0.565, 0.549),  //  21: Red
		Eigen::RowVector3d(0.902, 0.668, 0.497),  //  22: Red
		Eigen::RowVector3d(0.608, 0.863, 0.341),  //  23: light pink
		Eigen::RowVector3d(0.957, 0.941, 0.380),  //  24: light pink
		Eigen::RowVector3d(0.59, 0.85, 0.79),     //  25: Red

		Eigen::RowVector3d(0.392, 0.902, 0.902),  //  26: Red
		Eigen::RowVector3d(0.565, 0.733, 0.804),  //  27: Red
		Eigen::RowVector3d(1.000, 0.604, 0.830),  //  28: light pink
		Eigen::RowVector3d(0.950, 0.740, 0.100),  //  29: light pink
		Eigen::RowVector3d(0.118, 0.929, 0.365),  //  30: Red

		Eigen::RowVector3d(0.565, 0.804, 0.659),  //  31: light pink
        Eigen::RowVector3d(1.000, 0.915, 0.500)   //  32: light pink

};

void generateUnsaturatedColor(int& r, int& g, int& b, float minSaturation, float maxSaturation, float minValue, float maxValue) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, 359);

    // Generate a random hue value (0-360) for an unsaturated color
    int h = distribution(gen);

    // Generate random saturation and value within the specified range
    std::uniform_real_distribution<float> satDistribution(minSaturation, maxSaturation);
    std::uniform_real_distribution<float> valDistribution(minValue, maxValue);
    float s = satDistribution(gen);
    float v = valDistribution(gen);

    // Convert HSV to RGB
    float c = s * v;
    float x = c * (1 - std::abs(std::fmod(h / 60.0f, 2) - 1));
    float m = v - c;

    float r1, g1, b1;
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

    // Scale to 0-255 range and convert to integers
    r = static_cast<int>((r1 + m) * 255);
    g = static_cast<int>((g1 + m) * 255);
    b = static_cast<int>((b1 + m) * 255);
}

Window::Window(int originX, int originY, int width, int height) {
	this->originX = originX;
	this->originY = originY;
	this->width = width;
	this->height = height;
}

void Window::addPiece(const Piece &piece) {
	pieces.push_back(piece);
}

// this just adds the individual pieces of the shape to the window
// and then adds remaining pixels as another piece
void Window::addShape(Shape &shape) 
{
	pieces.clear();
	for (auto piece : shape.pieces) {
		pieces.push_back(piece);
	}

	vector<Pixel> remaining_pixels;
	accessibilityList.clear();
	for (const auto &px : shape.all_pixels) {
		if (!isPixelInPieces(px, shape.pieces)) {
			remaining_pixels.push_back(px);
			int currPixelID = shape.getPixelID(px);
			accessibilityList.push_back(shape.accessibilityList[currPixelID]);
		} 
	}

	pieces.push_back(Piece(remaining_pixels, Color{150,150,150}));

	solution = shape;
}

void Window::addLabel(const string& label) {
	this->label = label;
}

void Window::clear()
{
	pieces.clear();
}

Viewer::Viewer() {
	viewer.plugins.push_back(&plugin);
	plugin.widgets.push_back(&menu);
}

void Viewer::addWindow(Window window) {
	windows.push_back(window);
}

void Viewer::clearViewerMeshes()
{
	int dataListMeshNum = viewer.data_list.size();

	for (int i = dataListMeshNum - 1; i>0; i--)
	{
		viewer.data_list[i].meshgl.free();
		viewer.data_list.erase(viewer.data_list.begin() + i);
	}
}

void Viewer::updateWindows()
{
	Eigen::Matrix2d rotationMatrix;
    rotationMatrix << 0, -1,
                      1,  0;

	Eigen::Matrix3d rotationMatrix_3d;
    rotationMatrix_3d << 0, -1, 0,
                         1,  0, 0,
						 0, 0, 1;

	int maxPieceClusterID = 0;

	srand(3);

	for (int i = 0; i < windows.size(); i++)
	{
		for (int j = 0; j < windows[i].pieces.size(); j++)
		{
			Piece* piece = &(windows[i].pieces[j]);

			if (piece->pixels.size() == 0)
				continue;

			if (piece->pieceClusterID > maxPieceClusterID) 
				maxPieceClusterID = piece->pieceClusterID;
		}
	}

	vector<Eigen::RowVector3d> colorList;
	if (maxPieceClusterID > 32)
	{
		int r, g, b;
		float minSaturation = 0.25f; // Minimum saturation
    	float maxSaturation = 0.5f; // Maximum saturation
    	float minValue = 0.85f;      // Minimum value
    	float maxValue = 1.0f;      // Maximum value

		for (int i = 0; i < 200; i++)
		{
			generateUnsaturatedColor(r, g, b, minSaturation, maxSaturation, minValue, maxValue);

			Eigen::RowVector3d currColor;
			currColor << double(r) / 255.0, double (g) / 255.0, double(b) / 255.0; 

			colorList.push_back(currColor);
		}
	}

	vector<int> templateMeshIds;

	// set mesh data
	for (int i = 0; i < windows.size(); i++) {

		vector<int> templateBucket(maxPieceClusterID, 0);
		vector<Eigen::MatrixXd> templateVerticesList;
		vector<Eigen::MatrixXi> templateFacesList;
		vector<Eigen::MatrixXd> templateBoundaryVerticesList;
		vector<Eigen::MatrixXi> templateBoundaryFacesList;
		vector<tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>> templateEdgesList;
		vector<int> templateIDList;

		double currMinX = 10, currMinY = 10, currMaxX = -10, currMaxY = -10;

		for (int j = 0; j < windows[i].pieces.size(); j++) 
		{
			// first generate render mesh for this shape
			Piece* piece = &(windows[i].pieces[j]);

			if (piece->pixels.size() == 0)
				continue;
				
			piece->generateRenderMeshes();

			// then set mesh data in the viewer, and add edges for each pixel
			piece->viewer_mesh_id = viewer.append_mesh(false);

			piece->V = piece->V * rotationMatrix;
			piece->V1 = piece->V1 * rotationMatrix;
			piece->V2 = piece->V2 * rotationMatrix;
			piece->V3 = piece->V3 * rotationMatrix;
			piece->V4 = piece->V4 * rotationMatrix;

			viewer.data().set_mesh(piece->V, piece->F);
			viewer.data().show_lines = false;
			viewer.data().face_based = true;
			viewer.data().double_sided = true;
			viewer.data().shininess = 0.0f;

			if (windows[i].pieces.size() == 2)
			{
				viewer.data().set_colors(Eigen::RowVector3d(0.85, 0.85, 0.85));
			}
			else
			{
				if (piece->pieceClusterID <= 32)
				{
					viewer.data().set_colors(colorTable[(piece->pieceClusterID - 1) % 32]);
				}

				else
				{
					viewer.data().set_colors(colorList[(piece->pieceClusterID - 33) % 32]);
				}
			}

			viewer.data().add_edges(piece->V1, piece->V2, Eigen::RowVector3d(0.6, 0.6, 0.6));
			viewer.data().add_edges(piece->V2, piece->V3, Eigen::RowVector3d(0.6, 0.6, 0.6));
			viewer.data().add_edges(piece->V3, piece->V4, Eigen::RowVector3d(0.6, 0.6, 0.6));
			viewer.data().add_edges(piece->V4, piece->V1, Eigen::RowVector3d(0.6, 0.6, 0.6));
			viewer.data().double_sided = true;
			viewer.data().line_width = 0.25;

			// compute bounding box
			if (piece->V.rows() > 0)
			{
				if (piece->V.col(0).minCoeff() < currMinX)
					currMinX = piece->V.col(0).minCoeff();

				if (piece->V.col(1).minCoeff() < currMinY)
					currMinY = piece->V.col(1).minCoeff();

				if (piece->V.col(0).maxCoeff() > currMaxX)
					currMaxX = piece->V.col(0).maxCoeff();

				if (piece->V.col(1).maxCoeff() > currMaxY)
					currMaxY = piece->V.col(1).maxCoeff();
			}

			// set the boundary mesh
			piece->viewer_boundary_mesh_id = viewer.append_mesh(false);
			piece->V_boundary = piece->V_boundary * rotationMatrix_3d;

			viewer.data().set_mesh(piece->V_boundary, piece->F_boundary);
			viewer.data().show_lines = false;
			viewer.data().set_colors(Eigen::RowVector3d(0, 0, 0));
			viewer.data().double_sided = true;

			// save to template list 
			if (templateBucket[piece->pieceClusterID - 1] == 0)
			{
				templateIDList.push_back(piece->pieceClusterID);

				templateVerticesList.push_back(piece->V);
				templateFacesList.push_back(piece->F);

				templateEdgesList.push_back(tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>(piece->V1, piece->V2, piece->V3, piece->V4));

				templateBoundaryVerticesList.push_back(piece->V_boundary);
				templateBoundaryFacesList.push_back(piece->F_boundary);

				templateBucket[piece->pieceClusterID - 1] = 1;
			}
		}

		double currMinX_ = 10000, currMinY_ = 10000, currMaxX_ = -10000, currMaxY_ = -10000;

		// display templates below tiling solution
		if (templateIDList.size() > 1)
		{
			for (int k = 0; k < maxPieceClusterID; ++k)
			{
				int currClusterID = templateIDList[k];

				Eigen::MatrixXd currPieceV = templateVerticesList[k];
				Eigen::MatrixXi currPieceF = templateFacesList[k];

				Eigen::MatrixXd currV_1 = std::get<0>(templateEdgesList[k]);
				Eigen::MatrixXd currV_2 = std::get<1>(templateEdgesList[k]);
				Eigen::MatrixXd currV_3 = std::get<2>(templateEdgesList[k]);
				Eigen::MatrixXd currV_4 = std::get<3>(templateEdgesList[k]);

				Eigen::MatrixXd currBoundaryV = templateBoundaryVerticesList[k];
				Eigen::MatrixXi currBoundaryF = templateBoundaryFacesList[k];

				double centerY = currMinY - 6 * (floor((currClusterID - 1) / 10) + 1);
				double centerX = currMinX + 6 * ((currClusterID - 1) % 10);

				Eigen::RowVectorXd currCenter(2);
				currCenter << centerX , centerY;

				Eigen::RowVectorXd currCenter_3d(3);
				currCenter_3d << centerX , centerY, 0 ;

				Eigen::RowVectorXd centroid = currPieceV.colwise().mean();
				Eigen::RowVectorXd centroid3d = Eigen::RowVector3d(centroid(0), centroid(1), 0);

				// set shape mesh
				{
					currPieceV.rowwise() -= centroid;
					currPieceV.rowwise() += currCenter;

					currV_1.rowwise() -= centroid;
					currV_1.rowwise() += currCenter;

					currV_2.rowwise() -= centroid;
					currV_2.rowwise() += currCenter;

					currV_3.rowwise() -= centroid;
					currV_3.rowwise() += currCenter;

					currV_4.rowwise() -= centroid;
					currV_4.rowwise() += currCenter;

					int templateMeshId = viewer.append_mesh(false);
					templateMeshIds.push_back(templateMeshId);

					viewer.data().clear();
					viewer.data().set_mesh(currPieceV, currPieceF);
					viewer.data().show_lines = false;
					viewer.data().face_based = true;
					viewer.data().double_sided = true;
					viewer.data().shininess = 0;

					if (currClusterID <= 32)
					{
						viewer.data().set_colors(colorTable[(currClusterID - 1) % 32]);
					}

					else
					{
						viewer.data().set_colors(colorList[(currClusterID - 33) % 32]);
					}

					viewer.data().add_edges(currV_1, currV_2, Eigen::RowVector3d(0, 0, 0));
					viewer.data().add_edges(currV_2, currV_3, Eigen::RowVector3d(0, 0, 0));
					viewer.data().add_edges(currV_3, currV_4, Eigen::RowVector3d(0, 0, 0));
					viewer.data().add_edges(currV_4, currV_1, Eigen::RowVector3d(0, 0, 0));

					// compute bounding box
					if (currPieceV.rows() > 0)
					{
						if (currPieceV.col(0).minCoeff() < currMinX_)
							currMinX_ = currPieceV.col(0).minCoeff();

						if (currPieceV.col(1).minCoeff() < currMinY_)
							currMinY_ = currPieceV.col(1).minCoeff();

						if (currPieceV.col(0).maxCoeff() > currMaxX_)
							currMaxX_ = currPieceV.col(0).maxCoeff();

						if (currPieceV.col(1).maxCoeff() > currMaxY_)
							currMaxY_ = currPieceV.col(1).maxCoeff();
					}
				}

				// set boundary mesh
				{
					currBoundaryV.rowwise() -= centroid3d;
					currBoundaryV.rowwise() += currCenter_3d;

					int templateMeshId = viewer.append_mesh(false);
					templateMeshIds.push_back(templateMeshId);

					viewer.data().clear();
					viewer.data().set_mesh(currBoundaryV, currBoundaryF);
					viewer.data().show_lines = false;
					viewer.data().set_colors(Eigen::RowVector3d(0, 0, 0));
					viewer.data().double_sided = true;
					viewer.data().shininess = 0;
				}
			}
		}
		
		// add label for this window
		viewer.data().add_label(Eigen::Vector2d(4, 12), windows[i].label);
		viewer.data().show_custom_labels = true;;

		// set best camera view
		Eigen::MatrixXd bBox(4, 2);
		bBox << 	currMinX, currMinY,
					currMinX, currMaxY,	
					currMaxX, currMaxY,
					currMaxX, currMinY;

		Eigen::MatrixXi bBoxF(2, 3);
		bBoxF << 	0, 2, 1,
					0, 3, 2;

		viewer.core_list[i].align_camera_center(bBox, bBoxF);
		viewer.core_list[i].camera_zoom = 0.5;

		// set best camera view
		Eigen::MatrixXd bBox_(4, 2);
		bBox_ << 	currMinX_, currMinY_,
					currMinX_, currMaxY_,	
					currMaxX_, currMaxY_,
					currMaxX_, currMinY_;

		Eigen::MatrixXi bBoxF_(2, 3);
		bBoxF_ << 	0, 2, 1,
					0, 3, 2;

		if (i == 2)
		{
			viewer.core_list[i].align_camera_center(bBox_, bBoxF_);
			viewer.core_list[i].camera_zoom = 2.0;
		}
	}

	// set visibility of meshes in each window
	for (int i = 0; i < 2; i++) {
		int window_id = viewer.core_list[i].id;
		for (int j = 0; j < windows[i].pieces.size(); j++) {
			viewer.data(windows[i].pieces[j].viewer_mesh_id).set_visible(true, window_id);
			viewer.data(windows[i].pieces[j].viewer_boundary_mesh_id).set_visible(true, window_id);
		}
	}

	for (int i = 0; i < templateMeshIds.size(); ++i)
	{
		viewer.data(templateMeshIds[i]).set_visible(true, viewer.core_list[2].id);
	}
}

void Viewer::launch(int w, int h) 
{
	// set up views
	Window* firstWindow = &(windows[0]);
	viewer.core().viewport = Eigen::Vector4f(firstWindow->originX, firstWindow->originY, firstWindow->width, firstWindow->height);
	
	for (int i = 1; i < windows.size(); i++) viewer.append_core(Eigen::Vector4f(windows[i].originX, windows[i].originY, windows[i].width, windows[i].height), true);

	viewer.core_list[0].orthographic = true;
	viewer.core_list[1].orthographic = true;
	viewer.core_list[2].orthographic = true;

	viewer.core_list[0].background_color.setOnes();
	viewer.core_list[1].background_color.setOnes();
	viewer.core_list[2].background_color.setOnes();

	viewer.core_list[0].light_position = Eigen::Vector3f(100,100,100);
	viewer.core_list[1].light_position = Eigen::Vector3f(100,100,100);
	viewer.core_list[2].light_position = Eigen::Vector3f(100,100,100);

	viewer.core_list[0].rotation_type = igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION;
	viewer.core_list[1].rotation_type = igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION;
	viewer.core_list[2].rotation_type = igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION;

	menu.callback_draw_viewer_window = []() {};

	//launch 
	viewer.launch(false, "Inverse Tiling", w, h);

}

void Viewer::updateViewerBufferAndSave(vector<Piece> all_templates, Shape state, int pieceNum, string filepath, int stateID, vector<Eigen::RowVector3d> additionalColorList)
{
	Eigen::Matrix2d rotationMatrix;
    rotationMatrix << 0, -1,
                      1,  0;

	Eigen::Matrix3d rotationMatrix_3d;
    rotationMatrix_3d << 0, -1, 0,
                         1,  0, 0,
						 0, 0, 1;

	// for template display
	vector<Eigen::MatrixXd> templateVerticesList;
	vector<Eigen::MatrixXi> templateFacesList;
	vector<Eigen::MatrixXd> templateBoundaryVerticesList;
	vector<Eigen::MatrixXi> templateBoundaryFacesList;
	vector<tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>> templateEdgesList;
	vector<int> templateIDList;

	double currMinX = 10, currMinY = 10, currMaxX = -10, currMaxY = -10;

	// join matrix
	vector<Eigen::MatrixXd> joinedVMatrixList;
	vector<Eigen::MatrixXi> joinedFMatrixList;
	vector<Eigen::MatrixXd> joinedCMatrixList;

	int existingFaceNum = 0;

	vector<Eigen::MatrixXd> v1_list;
	vector<Eigen::MatrixXd> v2_list;

	for (int j = 0; j < windows[0].pieces.size(); j++)
	{
		// first generate render mesh for this shape
		Piece* piece = &(windows[0].pieces[j]);

		if (piece->pixels.size() == 0)
			continue;

		piece->generateRenderMeshes();

		// set the boundary mesh
		piece->V_boundary = piece->V_boundary * rotationMatrix_3d;

		Eigen::MatrixXd V_boundary_2d = piece->V_boundary.block(0, 0, piece->V_boundary.rows(), 2);
		Eigen::MatrixXi F_boundary_2d = piece->F_boundary;

		joinedVMatrixList.push_back(V_boundary_2d);

		Eigen::MatrixXi reorderedF_boundary = F_boundary_2d;
		reorderedF_boundary.array() += existingFaceNum;

		joinedFMatrixList.push_back(reorderedF_boundary);

		Eigen::RowVectorXd rowValues_boundary(3);
		rowValues_boundary << 0.0, 0.0, 0.0;

		Eigen::MatrixXd currColorMat_boundary(piece->V_boundary.rows(), 3);
		currColorMat_boundary = rowValues_boundary.replicate(piece->V_boundary.rows(), 1);
		joinedCMatrixList.push_back(currColorMat_boundary);

		existingFaceNum += piece->V_boundary.rows();

		piece->V1 = piece->V1 * rotationMatrix;
		piece->V2 = piece->V2 * rotationMatrix;
		piece->V3 = piece->V3 * rotationMatrix;
		piece->V4 = piece->V4 * rotationMatrix;

		v1_list.push_back(piece->V1); 
		v1_list.push_back(piece->V2);
		v1_list.push_back(piece->V3);
		v1_list.push_back(piece->V4);

		v2_list.push_back(piece->V2); 
		v2_list.push_back(piece->V3);
		v2_list.push_back(piece->V4);
		v2_list.push_back(piece->V1);
	}

	int totalGridRows = 0;
	for (auto & grid : v1_list)
	{
		totalGridRows += grid.rows();
	}

	cout << "totalGridRows: " << totalGridRows << endl;

	Eigen::MatrixXd V_grid = Eigen::MatrixXd(totalGridRows * 4, 2);
    Eigen::MatrixXi F_grid = Eigen::MatrixXi(totalGridRows * 2, 3);

	int Vrow = 0;
    int Frow = 0;
    double shortSide = 0.015;

	for (int j = 0; j < v1_list.size(); ++j)
	{
		for (int i = 0; i < v1_list[j].rows(); ++i)
		{
			Eigen::Vector2d topLeft, topRight, bottomLeft, bottomRight;

			Eigen::Vector2d p1;
			p1 << v1_list[j](i, 0), v1_list[j](i, 1); 

			Eigen::Vector2d p2;
			p2 << v2_list[j](i, 0), v2_list[j](i, 1);

			// cout << p1 << endl << p2 << endl;
			// cout << "Load 2 edge end points. " << endl;

			p1 = p1 + 0.5 * shortSide * (p1 - p2);
			p2 = p2 + 0.5 * shortSide * (p2 - p1);
			
			Eigen::Vector2d midPoint = (p1 + p2) / 2.0;
			Eigen::Vector2d e1 = (p2 - p1) / 2.0;
			Eigen::Vector2d e2(-e1.y(), e1.x());

			e2.normalize();
			e2 = (0.5 * shortSide) * e2;

			topLeft = midPoint + e2 - e1;
			topRight = midPoint + e2 + e1;
			bottomLeft = midPoint - e2 - e1;
			bottomRight = midPoint - e2 + e1;

			V_grid.row(Vrow)     << topLeft(0), topLeft(1);
			V_grid.row(Vrow + 1) << bottomLeft(0), bottomLeft(1);
			V_grid.row(Vrow + 2) << bottomRight(0), bottomRight(1);
			V_grid.row(Vrow + 3) << topRight(0), topRight(1);

			F_grid.row(Frow) << Vrow, Vrow + 1, Vrow + 2;
			F_grid.row(Frow + 1) << Vrow, Vrow + 3, Vrow + 2;

			Vrow += 4;
			Frow += 2;
		}
	}

	joinedVMatrixList.push_back(V_grid);

	Eigen::MatrixXi F_grid_reordered = F_grid;
	F_grid_reordered.array() += existingFaceNum;
	joinedFMatrixList.push_back(F_grid_reordered);

	Eigen::RowVectorXd rowValues_boundary(3);
	rowValues_boundary << 0.6, 0.6, 0.6;

	Eigen::MatrixXd currColorMat_grid(V_grid.rows(), 3);
	currColorMat_grid = rowValues_boundary.replicate(V_grid.rows(), 1);
	joinedCMatrixList.push_back(currColorMat_grid);

	existingFaceNum += V_grid.rows();

	for (int j = 0; j < windows[0].pieces.size(); j++) 
	{
		// first generate render mesh for this shape
		Piece* piece = &(windows[0].pieces[j]);

		if (piece->pixels.size() == 0)
			continue;
				
		piece->generateRenderMeshes();

		// then set mesh data in the viewer, and add edges for each pixel
		piece->V = piece->V * rotationMatrix;

		joinedVMatrixList.push_back(piece->V);

		Eigen::MatrixXi reorderedF = piece->F;
		reorderedF.array() += existingFaceNum;
		joinedFMatrixList.push_back(reorderedF);

		Eigen::MatrixXd currColorMat(piece->V.rows(), 3);
		Eigen::RowVectorXd rowValues(3);
		rowValues << 0.85, 0.85, 0.85;

		if (windows[0].pieces.size() - 1 == 1 or j >= pieceNum or stateID == 0)
		{
			rowValues << 0.85, 0.85, 0.85;
		}
		else
		{
			int currPieceClusterID = -1;
			int clusterCount = 0;
			for (auto & currTemplate : all_templates)
			{
				if (isSameShape(windows[0].pieces[j], currTemplate))
				{
					currPieceClusterID = clusterCount;
					break;
				}	
				clusterCount ++;
			}

			cout << "currPieceClusterID: " << currPieceClusterID << endl;

			if (currPieceClusterID == -1)
			{
				rowValues << 0.85, 0.85, 0.85;
			}
			else
			{
				// if (currPieceClusterID <= stateID - 1)
				// {
				// 	rowValues << additionalColorList[clusterCount][0], additionalColorList[clusterCount][1], additionalColorList[clusterCount][2];
				// }
				// else
				// {
				// 	rowValues << 0.85, 0.85, 0.85;
				// }

				rowValues << additionalColorList[clusterCount][0], additionalColorList[clusterCount][1], additionalColorList[clusterCount][2];
				
				if (currPieceClusterID < 32)
					rowValues << colorTable[currPieceClusterID % 32][0], colorTable[currPieceClusterID % 32][1], colorTable[currPieceClusterID % 32][2];

				else
				{
					rowValues << additionalColorList[(currPieceClusterID - 32) % 200][0], additionalColorList[(currPieceClusterID - 32) % 200][1], additionalColorList[(currPieceClusterID - 32) % 200][2];
				}
			}
		}

		currColorMat = rowValues.replicate(piece->V.rows(), 1);
		joinedCMatrixList.push_back(currColorMat);

		existingFaceNum += piece->V.rows();

		// compute bounding box
		if (piece->V.rows() > 0)
		{
			if (piece->V.col(0).minCoeff() < currMinX)
				currMinX = piece->V.col(0).minCoeff();

			if (piece->V.col(1).minCoeff() < currMinY)
				currMinY = piece->V.col(1).minCoeff();

			if (piece->V.col(0).maxCoeff() > currMaxX)
				currMaxX = piece->V.col(0).maxCoeff();

			if (piece->V.col(1).maxCoeff() > currMaxY)
				currMaxY = piece->V.col(1).maxCoeff();
		}
	}

	int totalVertices = 0;
	int totalFaces = 0;
	int totalColors = 0;
	for (const auto& v : joinedVMatrixList) 
		totalVertices += v.rows();
	for (const auto& f : joinedFMatrixList) 
		totalFaces += f.rows();
	for (const auto& c : joinedCMatrixList) 
		totalColors += c.rows();

	cout << "totalVertices: " << totalVertices << endl;
	cout << "totalFaces: " << totalFaces << endl;
	cout << "totalColors: " << totalColors << endl;

	// 合并顶点矩阵
	Eigen::MatrixXd combinedVertices(totalVertices, 2);
	int rowOffset = 0;
	for (const auto& v : joinedVMatrixList) {
		combinedVertices.block(rowOffset, 0, v.rows(), 2) = v;
		rowOffset += v.rows();
	}

	// 合并面矩阵
	Eigen::MatrixXi combinedFaces(totalFaces, 3);
	int faceRowOffset = 0;
	for (const auto& f : joinedFMatrixList) {
		combinedFaces.block(faceRowOffset, 0, f.rows(), 3) = f;
		faceRowOffset += f.rows();
	}

	// 合并颜色矩阵
	Eigen::MatrixXd combinedColors(totalColors, 3);
	int colorRowOffset = 0;
	for (const auto& c : joinedCMatrixList) {
		combinedColors.block(colorRowOffset, 0, c.rows(), 3) = c;
		colorRowOffset += c.rows();
	}

	cout << "combine matrix done. " << endl;

	viewer.data().clear();
	viewer.data().set_mesh(combinedVertices, combinedFaces);
	viewer.data().double_sided = true;
	viewer.data().show_lines = false;
	viewer.data().set_colors(combinedColors);
	viewer.data().is_visible = true;

	// viewer.data().clear_edges();
	// for (int i = 0; i < v1_list.size(); ++i)
	// {
	// 	viewer.data().add_edges(v1_list[i], v2_list[i], Eigen::RowVector3d(0.6, 0.6, 0.6));
	// }

	// set best camera view
	if (stateID == 0)
	{
		Eigen::MatrixXd bBox(4, 2);
		bBox << 	currMinX, currMinY,
					currMinX, currMaxY,	
					currMaxX, currMaxY,
					currMaxX, currMinY;

		Eigen::MatrixXi bBoxF(2, 3);
		bBoxF << 	0, 2, 1,
					0, 3, 2;

		viewer.core_list[0].align_camera_center(bBox, bBoxF);
		viewer.core_list[0].camera_zoom = 0.5;
	}
	
}
