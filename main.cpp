#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include "Piece.h"
#include "Shape.h"
#include "Viewer.h"
#include "Dissection.h"
#include "Util.h"
#include <iostream>
#include <filesystem>
#include <cstdio>
#include <chrono>
#include <thread>
#include <igl/stb/write_image.h>
#include <igl/stb/read_image.h>
#include <queue>

using namespace std::chrono;
using namespace std;
namespace fs = std::filesystem;

// Color scheme table
Eigen::RowVector3d helperColorList[32] = {
		// Eigen::RowVector3d(0.64, 0.9, 0.9),       //  6: light cyan
		// Eigen::RowVector3d(0.9, 0.63, 0.63),      //  7: light red

		Eigen::RowVector3d(0.63, 0.9, 0.63),	  //  1: light green 
		// Eigen::RowVector3d(0.9, 0.63, 0.63),      //  7: light red
		Eigen::RowVector3d(0.95, 0.95, 0.66),     //  2: light yellow
		Eigen::RowVector3d(0.81, 0.63, 0.9),   	  //  3: light purple
        Eigen::RowVector3d(0.9, 0.74, 0.63),      //  4: light brown
        Eigen::RowVector3d(0.63, 0.63, 0.9),      //  5: light blue

		Eigen::RowVector3d(0.64, 0.9, 0.9),       //  6: light cyan
		Eigen::RowVector3d(0.9, 0.63, 0.63),      //  7: light red
		Eigen::RowVector3d(0.9, 0.63, 0.8),       //  8: light pink
		Eigen::RowVector3d(0.8, 1.0, 0.7),        //  9: light pink
		// Eigen::RowVector3d(0.85, 0.85, 0.85),     //  10: light pink
		Eigen::RowVector3d(0.7, 0.85, 1.0),    //  10: light pink

		Eigen::RowVector3d(0.392, 0.784, 0.706),     //  11: light pink
		Eigen::RowVector3d(0.8, 0.8, 0.56),       //  12: light pink
		Eigen::RowVector3d(0.95, 0.78, 0.67),     //  13: light pink
		Eigen::RowVector3d(0.582, 0.701, 0.9),    //  14: light pink
		Eigen::RowVector3d(0.714, 0.522, 0.8),    //  15: Red

		Eigen::RowVector3d(0.627, 0.949, 0.102),  //  16: Red
		Eigen::RowVector3d(0.927, 0.588, 0.948),  //  17: Red
		Eigen::RowVector3d(0.424, 0.965, 0.796),  //  18: light pink
		Eigen::RowVector3d(0.902, 0.510, 0.745),  //  19: light pink
		Eigen::RowVector3d(0.941, 0.941, 0.157),  //  20: Red

		Eigen::RowVector3d(0.900, 0.565, 0.549),  //  21: Red
		Eigen::RowVector3d(0.902, 0.668, 0.497),  //  22: Red
		Eigen::RowVector3d(0.608, 0.863, 0.341),  //  23: light pink
		Eigen::RowVector3d(0.957, 0.941, 0.380),  //  24: light pink
		Eigen::RowVector3d(0.59, 0.85, 0.79),  //  25: Red

		Eigen::RowVector3d(0.392, 0.902, 0.902),  //  26: Red
		Eigen::RowVector3d(0.565, 0.733, 0.804),  //  27: Red
		Eigen::RowVector3d(1.000, 0.604, 0.830),  //  28: light pink
		Eigen::RowVector3d(0.950, 0.740, 0.100),  //  29: light pink
		Eigen::RowVector3d(0.118, 0.929, 0.365),  //  30: RedEigen::RowVector3d(0.118, 0.929, 0.365)

		Eigen::RowVector3d(0.565, 0.804, 0.659),
        Eigen::RowVector3d(1.000, 0.915, 0.500)  //  31: light pink

};

// Parameters
int N = 15; // num of template pieces to solve for
int G = 10; // number of times to try growing process for each seed
int T = 120; // amount of time to run algorithm for in minutes
int K = 100; // number of discreate equivalence classes
int seedNum = 50;
bool isUniformDistribution = true;
int minSeedDist = 2;
int pieceNumThreshold = 2;
int shapeCandisNum = 10;
float minAvgSizeRatio = 0.5;
float maxAvgSizeRatio = 1.5;
int minTileSize = 4;
int maxTileSize = 6;
bool autoSave = true;
bool isLooseSizeRequirement = true;
// bool DEBUG = true;

// debug mode parameters
int growing_round = 0;
bool isDebugShapeDisconnected = false;
Shape debug_seed;
Shape debug_shape;

// For control panel
string inputFileName;  // full path
string fileName;       // file name only
string pureFileName;   // file name without extension

// for rendering tiling states
string renderInputFileName;
int stateID = 0;
int totalPieceNum;
vector<Piece> all_templates;
Shape terminalState;

// Input & Output & Data
Shape input;
Shape best_solution;
double runningTime;
int bestScore;

// Init the viewer
Viewer viewer;
Window input_window = Window(0, 0, 1600, 1600);
Window output_window = Window(1600, 0, 1600, 1600);

// Allocate temporary buffers
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> _R(1600,1600);
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> _G(1600,1600);
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> _B(1600,1600);
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> _A(1600,1600);

vector<Eigen::RowVector3d> currentColorList;

// Define directions for movement (up, down, left, right)
vector<pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

// Function to check if the piece is connected
bool isPieceConnected(const vector<Pixel>& currPixelList) {
    if (currPixelList.empty()) return true;  // An empty piece is trivially connected

    set<pair<int, int>> visited;
    queue<Pixel> toVisit;

    // Start BFS/DFS from the first pixel in the list
    toVisit.push(currPixelList[0]);
    visited.insert({currPixelList[0].x, currPixelList[0].y});

    while (!toVisit.empty()) {
        Pixel current = toVisit.front();
        toVisit.pop();

        // Explore all possible 4-connected neighbors (up, down, left, right)
        for (const auto& dir : directions) {
            int newX = current.x + dir.first;
            int newY = current.y + dir.second;

            // Check if this neighbor is part of the piece
            for (const auto& pixel : currPixelList) {
                if (pixel.x == newX && pixel.y == newY) {
                    if (visited.find({newX, newY}) == visited.end()) {
                        toVisit.push(Pixel(newX, newY));
                        visited.insert({newX, newY});
                    }
                }
            }
        }
    }

    // If the number of visited pixels matches the number of pixels in the piece, it's connected
    return (visited.size() == currPixelList.size()) || (visited.size() == currPixelList.size() - 1) || (visited.size() == currPixelList.size() - 2) || (visited.size() == currPixelList.size() - 3) || (visited.size() == currPixelList.size() - 4);
}

void statusBar()
{
    fileName = inputFileName;

    for (int i = fileName.size() - 1; i > 0; --i)
    {
        if (fileName[i] == '/')
        {
            fileName.erase(fileName.begin(), fileName.begin() + i + 1);
            break;
        }
    }

    size_t dotPosition = fileName.find_last_of('.');
    pureFileName = fileName.substr(0, dotPosition);

    ImGui::Text("Shape :     %s", fileName.c_str());
    ImGui::Text("Pixels :     %zu", input.all_pixels.size());
    ImGui::Text("Avg size:      %zu", input.all_pixels.size() / N);
    ImGui::Text("N         :     %lu", best_solution.pieces.size());
    ImGui::Text("K         :     %d", bestScore);
    ImGui::Text("minSize         :     %d", best_solution.getSmallestPieceSize());
    ImGui::Text("maxSize         :     %d", best_solution.getLargestPieceSize());
    ImGui::Text("Time   :     %f s", runningTime);
}

void updateViewer()
{
    viewer.windows[0].clear();
    viewer.windows[0].addShape(input);

    viewer.windows[1].clear();
    viewer.windows[1].addShape(best_solution);

    viewer.clearViewerMeshes();
    viewer.updateWindows();
}

void addControlPanel()
{
    // Draw additional windows
    viewer.menu.callback_draw_custom_window = [&]()
    {
        //// color preset:
        ImGui::GetStyle().Colors[ImGuiCol_WindowBg] = ImVec4(0.9f, 0.9f, 0.9f, 1.00f);
        ImGui::GetStyle().Colors[ImGuiCol_TitleBg] = ImVec4(0.9f, 0.9f, 0.9f, 1.00f);
        ImGui::GetStyle().Colors[ImGuiCol_TitleBgActive] = ImVec4(0.9f, 0.9f, 0.9f, 1.00f);
        ImGui::GetStyle().Colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.9f, 0.9f, 0.9f, 1.00f);

        ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);
        ImGui::GetStyle().Colors[ImGuiCol_FrameBgHovered] = ImVec4(1.0f, 1.0f, 1.0f, 0.50f);
        ImGui::GetStyle().Colors[ImGuiCol_FrameBgActive] = ImVec4(1.0f, 1.0f, 1.0f, 0.50f);

        ImGui::GetStyle().Colors[ImGuiCol_MenuBarBg] = ImVec4(0.9f, 0.9f, 0.9f, 1.00f);
        ImGui::GetStyle().Colors[ImGuiCol_Text] = ImVec4(0.0f, 0.0f, 0.0f, 1.00f);

        ImGui::GetStyle().Colors[ImGuiCol_PopupBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.00f);

        // Define next window position + size
        ImGui::SetNextWindowPos(ImVec2(1600, 0), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 800), ImGuiCond_FirstUseEver);
        ImGui::Begin(
            "Control Panel", nullptr,
            ImGuiWindowFlags_NoSavedSettings
        );

        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        statusBar();

        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        double p = 2;
        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        ImGui::Text("Seed num");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Seed num", &seedNum);

        ImGui::Text("N piece num");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##N piece num", &N);

        ImGui::Text("Target K");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Target K", &K);

        ImGui::Text("Min seed dist");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Min seed dist", &minSeedDist);

        ImGui::Text("G grow attempts per seed");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##G growing attempts", &G);

        ImGui::Text("T max mins");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##T max mins", &T);

        ImGui::Text("Shape Candis Num");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Shape Candis Num", &shapeCandisNum);

        // ImGui::Text("Min Size Ratio");
        // ImGui::SameLine(100, p);
        // ImGui::DragFloat("##Min Size Ratio", &minAvgSizeRatio);

        // ImGui::Text("Max Size Ratio");
        // ImGui::SameLine(100, p);
        // ImGui::DragFloat("##Max Size Ratio", &maxAvgSizeRatio);

        ImGui::Text("Min Tile Size");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Min Tile Size", &minTileSize);

        ImGui::Text("Max Tile Size");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Max Tile Size", &maxTileSize);

        ImGui::Text("Min Piece Threshold");
        ImGui::SameLine(100, p);
        ImGui::DragInt("##Min Piece Threshold", &pieceNumThreshold);

        ImGui::Text("Uniform Seed Distribution");
        ImGui::SameLine(150, p);
        ImGui::Checkbox("##Uniform Seed Distribution", &isUniformDistribution);

        ImGui::Text("Loose Size Requirement");
        ImGui::SameLine(150, p);
        ImGui::Checkbox("##Loose Size Requirement", &isLooseSizeRequirement);

        ImGui::Text("Auto Save");
        ImGui::SameLine(150, p);
        ImGui::Checkbox("##Auto Save", &autoSave);

        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        // Add a button
        if (ImGui::Button("Read", ImVec2(-1,0)))
        {
            inputFileName = igl::file_dialog_open();
            if( !inputFileName.empty() )
            {
                input = readShapeFromFile(inputFileName);
                // updateAccessibilityAndBlockability(input);

                updateAccessibility(input);
                updateBlockability(input);

                best_solution = input;
                bestScore = input.current_templates.size();
                growing_round = 0;

                if (input.pieces.size() != 0)
                    N = input.pieces.size();
                else    
                    N = 15;

                updateViewer();
            }
        }

        if (ImGui::Button("Optimize -- Brute-force", ImVec2(-1,0)))
        {
            auto start = high_resolution_clock::now();
            int minK = 100000;
            std::vector<std::vector<int>> puzzleMat = input.get2DPuzzleMatrix(); 

            while (1)
            {
                bool isPieceValid = true;

                for (int i = 0; i < puzzleMat.size(); ++i)
                {
                    for (int j = 0; j < puzzleMat[i].size(); ++j)
                    {
                        if (puzzleMat[i][j] > 0)
                        {
                            puzzleMat[i][j] = rand() % N + 1;
                        }
                    }
                }

                // initialize each piece
                vector<Piece> currPieceList;
                for (int i = 0; i < N; ++i)
                {
                    vector<Pixel> currPixelList;

                    for (int j = 0; j < puzzleMat.size(); ++j)
                    {
                        for (int k = 0; k < puzzleMat[j].size(); ++k)
                        {
                            if (puzzleMat[j][k] == i + 1)
                            {
                                currPixelList.push_back(Pixel(j, k));
                            }
                        }
                    }

                    isPieceValid = isPieceConnected(currPixelList);

                    if (!isPieceValid)
                        break;

                    Piece currPiece = Piece(currPixelList);
                    currPieceList.push_back(currPiece);
                }

                if (!isPieceValid)
                    continue;

                // cout << "pieceListSize: " << currPieceList.size() << endl;

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

                auto time_now = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(time_now - start);
                double time_elapsed = (duration.count() / (double)1000) / 60;

                if (templates.size() < minK)
                {
                    minK = templates.size();

                    cout << "currMinK: " << minK << endl;   
                    cout << "time: " << time_elapsed << " mins" << endl;
                }

                // break;
            }
        }

        // Add a button
        if (ImGui::Button("Optimize -- Global", ImVec2(-1,0)))
        {
            vector<Piece> existingTemplates;
            bool isFound = false;
            best_solution = generatePieces(input, N, K, G, T, autoSave, pureFileName, 
                                           seedNum, minAvgSizeRatio, maxAvgSizeRatio, minTileSize, maxTileSize,
                                           minSeedDist, isUniformDistribution, shapeCandisNum,
                                           runningTime, bestScore, existingTemplates, isFound, isLooseSizeRequirement);

            K = best_solution.current_templates.size();

            // update viewport
            updateViewer();
        }

        // ImGui::Dummy(ImVec2(0.0f, 2.0f));

        // if (ImGui::Button("Read Shape with Seeds", ImVec2(-1,0)))
        // {
        //     inputFileName = igl::file_dialog_open();

        //     if( !inputFileName.empty() )
        //     {
        //         input = readSeedFromFile(inputFileName);
        //         updateAccessibilityAndBlockability(input);

        //         best_solution = input;
        //         bestScore = input.current_templates.size();
        //         growing_round = 0;

        //         if (input.pieces.size() != 0)
        //             N = input.pieces.size();
        //         else    
        //             N = 15;

        //         // updateViewer();
        //     }
        // }

        // if (ImGui::Button("Optimize with Seeds -- Global", ImVec2(-1,0)))
        // {
        //     vector<Piece> existingTemplates;
        //     bool isFound = false;
        //     best_solution = generatePiecesWithGivenSeeds(input, N, K, G, T, autoSave, pureFileName, 
        //                                    seedNum, minAvgSizeRatio, maxAvgSizeRatio, minTileSize, maxTileSize,
        //                                    minSeedDist, isUniformDistribution, shapeCandisNum,
        //                                    runningTime, bestScore, existingTemplates, isFound, isLooseSizeRequirement);

        //     K = best_solution.current_templates.size();

        //     // update viewport
        //     updateViewer();
        // }

        // ImGui::Dummy(ImVec2(0.0f, 2.0f));

        // if (ImGui::Button("Optimize Size -- Local", ImVec2(-1,0)))
        // {
        //     input = best_solution;
        //     best_solution = optimizePiecesSizeLocally(input, minTileSize, maxTileSize);

        //     // update viewport
        //     updateViewer();
        // }
        

        if (ImGui::Button("Optimize K -- Local", ImVec2(-1,0)))
        {
            input = best_solution;
            best_solution = generatePiecesLocally(input, pieceNumThreshold);

            // update viewport
            updateViewer();
        }

        if (ImGui::Button("Optimize Distribution -- Local", ImVec2(-1,0)))
        {
            input = best_solution;
            best_solution = generatePiecesLocally_distribution(input, pieceNumThreshold);

            // update viewport
            updateViewer();
        }

        if (ImGui::Button("Save Tiling Image", ImVec2(-1,0)))
        {
            for (int i = 0; i < viewer.viewer.data_list.size(); ++i)
                viewer.viewer.core().draw_buffer(viewer.viewer.data(),false,_R,_G,_B,_A);

            _A.setConstant(255);
            
            igl::stb::write_image("/home/yishin/CLionProjects/InverseTiling/Rendering/temp.png",_R,_G,_B,_A);
        }


        // ImGui::Dummy(ImVec2(0.0f, 2.0f));

        // // Add a button
        // if (ImGui::Button("Generate Seeds", ImVec2(-1,0)))
        // {
        //     srand(seedNum);

        //     input.clear();
        //     Shape seed = generateNSeedPieces(N, input, minSeedDist, isUniformDistribution);

        //     N = seed.pieces.size();
        //     best_solution = seed;

        //     // update viewport
        //     updateViewer();
        // }

        // // Add a button
        // if (ImGui::Button("Grow Pieces", ImVec2(-1,0)))
        // {
        //     if (growing_round == 0)
        //     {
        //         srand(seedNum);

        //         input.clear();

        //         debug_shape = generateNSeedPieces(N, input, minSeedDist, isUniformDistribution);
        //         // debug_shape.updateCentrality();

        //         best_solution = debug_shape;

        //         isDebugShapeDisconnected = false;

        //         growing_round++;

        //         // update viewport
        //         updateViewer(); 

        //     }
        //     else
        //     {
        //         vector<Piece> existingTemplates;
        //         debug_shape = growPiecesByOnePixel(debug_shape, shapeCandisNum, isDebugShapeDisconnected, existingTemplates, minTileSize, maxTileSize, isLooseSizeRequirement);

        //         if (isDebugShapeDisconnected)
        //         {
        //             cout << "Find a disconnected result. " << endl;

        //             isDebugShapeDisconnected = true;

        //             best_solution = debug_shape;

        //             growing_round = 0;
        //         }
        //         else
        //         {
        //             cout << "Grew all pieces by one pixel. " << debug_shape.getNumUnfilledPx() << " pixels remaining. \n";

        //             best_solution = debug_shape;

        //             growing_round++;
        //         }

        //         // update viewport
        //         updateViewer(); 
        //     }
        // }

        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        // Add a button
        if (ImGui::Button("Write", ImVec2(-1,0)))
        {
            string outputFolderPath = igl::file_dialog_save();
            if( !outputFolderPath.empty())
            {
                cout << outputFolderPath << endl;

                size_t lastSlashPosition = outputFolderPath.find_last_of('/');
                string savingPath = outputFolderPath.substr(0, lastSlashPosition + 1);
                string pureFolderName = best_solution.getOutputFolderPath(savingPath + pureFileName) + "_size" + to_string(best_solution.getSmallestPieceSize()) + "To" + to_string(best_solution.getLargestPieceSize());

                // save .puz file
                saveShape2File(best_solution, pureFolderName, best_solution.getOutputFileFullName(pureFolderName + "/" + pureFileName), runningTime);
            }
        }

        ImGui::Dummy(ImVec2(0.0f, 2.0f));

        if (ImGui::Button("Render states", ImVec2(-1,0)))
        {
            srand(1);

            renderInputFileName = igl::file_dialog_open();
            if( !renderInputFileName.empty() )
            {
                fs::path fullPath(renderInputFileName);

                fs::path directoryPath = fullPath.parent_path();

                std::cout << "Directory path: " << directoryPath.string() << std::endl;

                fs::create_directory(directoryPath / "state_images");

                stateID = 0;

                terminalState = readShapeFromFile(renderInputFileName);
                totalPieceNum = terminalState.pieces.size();

                all_templates.clear();
                currentColorList.clear();

                terminalState.updateTemplates(false);
                terminalState.assignPieceClusterID();

                for (int k = 1; k <= 20; ++k)
                {
                    for (auto & piece : terminalState.pieces)
                    {
                        if (piece.pieceClusterID == k)
                        {
                            all_templates.push_back(piece);
                            currentColorList.push_back(helperColorList[k - 1]);

                            cout << "find a template in first 20 templates. " << endl;

                            break;
                        }
                    }
                }

                int countColor = 20;

                while(fs::exists(directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")))
                {
                    Shape currState = readShapeFromFile((directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")).string());

                    int countPiece = 0;
                    for (auto & piece : currState.pieces)
                    {
                        if (currState.pieces.size() < totalPieceNum)
                            break;

                        countPiece++;

                        if (!isPieceinTemplates(piece, all_templates) and countPiece <= totalPieceNum)
                        {
                            all_templates.push_back(piece);

                            if (countColor < 32)
                            {
                                currentColorList.push_back(helperColorList[countColor]);    

                                countColor ++;
                            }

                            else
                            {
                                Color_eigen color;

                                float minSaturation = 0.25f; // Minimum saturation
                                float maxSaturation = 0.5f; // Maximum saturation
                                float minValue = 0.85f;      // Minimum value
                                float maxValue = 1.00f;      // Maximum value

                                generateUnsaturatedColor(color, minSaturation, maxSaturation, minValue, maxValue, currentColorList);

                                cout << "successfully to generate a new color. " << endl;

                                currentColorList.push_back(Eigen::RowVector3d(color[0] / 255.0, color[1] / 255.0, color[2] / 255.0));
                            }
                        }
                    }

                    stateID++;
                }

                cout << stateID << " states in total. " << endl;
                cout << all_templates.size() << " templates in total. " << endl;

                stateID = 0;

                // load the input domain first
                while(fs::exists(directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")))
                {
                    Shape currState = readShapeFromFile((directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")).string());

                    viewer.windows[0].clear();
                    viewer.windows[0].addShape(currState);

                    viewer.windows[1].clear();
                    viewer.windows[1].addShape(currState);

                    viewer.clearViewerMeshes();
                    viewer.updateViewerBufferAndSave(all_templates, currState, totalPieceNum, (directoryPath / "state_images" / (to_string(stateID) + ".png")).string(), stateID, currentColorList);

                    break;
                }
            }
        }

        if (ImGui::Button("Prev State", ImVec2(-1,0)))
        {
            fs::path fullPath(renderInputFileName);

            fs::path directoryPath = fullPath.parent_path();

            while(fs::exists(directoryPath / "enlarging_states" / (to_string(stateID - 1) + ".puz")))
            {
                stateID--;

                Shape currState = readShapeFromFile((directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")).string());

                viewer.windows[0].clear();
                viewer.windows[0].addShape(currState);

                viewer.windows[1].clear();
                viewer.windows[1].addShape(currState);

                viewer.clearViewerMeshes();
                viewer.updateViewerBufferAndSave(all_templates, currState, totalPieceNum, (directoryPath / "state_images" / (to_string(stateID) + ".png")).string(), stateID, currentColorList);

                break;
            }
        }

        if (ImGui::Button("Next State", ImVec2(-1,0)))
        {
            fs::path fullPath(renderInputFileName);

            fs::path directoryPath = fullPath.parent_path();

            while(fs::exists(directoryPath / "enlarging_states" / (to_string(stateID + 1) + ".puz")))
            {
                stateID++;

                Shape currState = readShapeFromFile((directoryPath / "enlarging_states" / (to_string(stateID) + ".puz")).string());

                viewer.windows[0].clear();
                viewer.windows[0].addShape(currState);

                viewer.windows[1].clear();
                viewer.windows[1].addShape(currState);

                viewer.clearViewerMeshes();
                viewer.updateViewerBufferAndSave(all_templates, currState, totalPieceNum, (directoryPath / "state_images" / (to_string(stateID) + ".png")).string(), stateID, currentColorList);

                break;
            }
        }

        if (ImGui::Button("Save Image", ImVec2(-1,0)))
        {
            fs::path fullPath(renderInputFileName);

            fs::path directoryPath = fullPath.parent_path();

            for (int i = 0; i < viewer.viewer.data_list.size(); ++i)
                viewer.viewer.core().draw_buffer(viewer.viewer.data(),false,_R,_G,_B,_A);

            _A.setConstant(255);

            cout << (directoryPath / "state_images" / (to_string(stateID) + ".png")).string() << endl;
            
            igl::stb::write_image(((directoryPath / "state_images" / (to_string(stateID) + ".png")).string()),_R,_G,_B,_A);
        }

        ImGui::End();
    };
}

void initWindows()
{   
    addControlPanel();

    viewer.addWindow(input_window);
    viewer.addWindow(output_window);

    // std::thread updateThread(checkFrameUpdate);

    viewer.launch(1900, 800);

    // 等待并行线程结束
    // updateThread.join();
}

int main() {

    initWindows();

    return 0;
}