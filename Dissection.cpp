#include "Dissection.h"
#include "Util.h"
#include <random>
#include <iostream>

using namespace std::chrono;
using namespace std;

Shape generatePiecesLocally_distribution(Shape & shape, int pieceNumThreshold)
{
    Shape newShape = shape;
    bool isSuccess = false;

    int minTileSize = newShape.getSmallestPieceSize();
    int maxTileSize = newShape.getLargestPieceSize();

    for(int k = 0; k < 30; ++k)
    {
        Shape currShape = newShape;
        bool noPossiblePieceTemplate;

        Piece eliminationPieceTemplate = selectEliminationPieceTemplate(currShape.current_templates, currShape, pieceNumThreshold, noPossiblePieceTemplate);
        cout << endl << eliminationPieceTemplate.printPixels() << endl;

        if (noPossiblePieceTemplate) 
        {
            cout << "no Possible Piece Template. " << endl;
            break;
        }

        vector<Piece> eliminationPieces = getTemplateInstances(eliminationPieceTemplate, currShape);

        // cout << "eliminationPieces: " << eliminationPieces.size() << endl;

        // select one piece from them
        vector<float> possibList(eliminationPieces.size(), 1);

        vector<Piece> oneEliminationPieceList;
        oneEliminationPieceList.push_back(eliminationPieces[GetRandomObjIndex(possibList, 1)]);

        int N, K;
        vector<int> eliminationPieceIDList;
    
        Shape subShape = getSubShape(oneEliminationPieceList, currShape, N, K, eliminationPieceIDList);
        subShape.updateTemplates(false);
        updateAccessibilityAndBlockability(subShape);

        // cout << "Subshape: N = " << N << ", K = " << K << endl;  
        // cout << endl << subShape.printAllPixels() << endl;

        int original_K = currShape.current_templates.size();
        double original_var = currShape.getInstanceNumDeviation();

        // cout << "before delete K: " << currShape.current_templates.size() << endl;
        // cout << "before delete N: " << currShape.pieces.size() << endl;

        sort(eliminationPieceIDList.begin(), eliminationPieceIDList.end());
        for (int i = eliminationPieceIDList.size() - 1; i >= 0; --i)
        {
            currShape.pieces.erase(currShape.pieces.begin() + eliminationPieceIDList[i]);
        }

        currShape.updateTemplates(false);

        // cout << "after delete K: " << currShape.current_templates.size() << endl;
        // cout << "after delete N: " << currShape.pieces.size() << endl;

        double timeLimitForEachSubshape = 10.0;
        double runningTime;
        int bestScore;

        auto start = high_resolution_clock::now();

        int seed = 0;
        while (1)
        {
            bool isFound = false;
            Shape newSubShape = generatePieces(subShape, N, K, 15, 2, false, " ",
                                           seed, 0.3, 3.0, minTileSize, maxTileSize, 2, true, 10, runningTime, bestScore, currShape.current_templates, isFound, false);

            if (isFound)
            {
                vector<Piece> updatedTemplates = currShape.current_templates;

                for (auto & piece : newSubShape.current_templates)
                {
                    if (!isPieceinTemplates(piece, updatedTemplates))
                    {
                        updatedTemplates.push_back(piece);
                    }
                }

                cout << "original_K: " << original_K << endl;
                cout << "after_K: " << updatedTemplates.size() << endl;

                if (updatedTemplates.size() == original_K)
                // if (1)
                {
                    // cout << "successfully reduce " << original_K - updatedTemplates.size() << " templates. " << endl;

                    for (auto & piece : newSubShape.pieces)
                    {
                        currShape.pieces.push_back(piece);
                    }

                    Shape tempShape;
                    tempShape = currShape;

                    tempShape.updateTemplates(false);
                    tempShape.assignPieceID();
                    tempShape.assignPieceClusterID();

                    double after_var = tempShape.getInstanceNumDeviation();

                    cout << "original_var: " << original_var << endl;
                    cout << "after_var: " << after_var << endl;

                    if (after_var < original_var)
                    {                      
                        cout << "original_var: " << original_var << endl;
                        cout << "after_var: " << after_var << endl;

                        newShape = currShape;
                        newShape.updateTemplates(false);
                        newShape.assignPieceID();
                        newShape.assignPieceClusterID();

                        isSuccess = true;

                    break;
                    }
                }

                // if (updatedTemplates.size() == original_K)
                // {
                //     if (newSubShape.getSmallestPieceSize() > minTileSize or newSubShape.getLargestPieceSize() < maxTileSize)
                //     {
                //         newShape = currShape;
                //         newShape.updateTemplates();
                //         newShape.assignPieceID();
                //         newShape.assignPieceClusterID();
                //     }

                //     isSuccess = true;

                //     break;
                // }
            }

            auto time_now = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(time_now - start);
            double time_elapsed = (duration.count() / (double)1000) / 60;

            if (time_elapsed > timeLimitForEachSubshape)
            {
                cout << "unsuccessfuly to reduce K. " << endl;
                break;
            }

            ++seed;         
        }

        if (isSuccess)
            break;
    }

    return newShape;
}

Shape generatePiecesLocally(Shape & shape, int pieceNumThreshold)
{
    Shape newShape = shape;
    bool isSuccess = false;

    int minTileSize = newShape.getSmallestPieceSize();
    int maxTileSize = newShape.getLargestPieceSize();

    for(int k = 0; k < 20; ++k)
    {
        Shape currShape = newShape;
        bool noPossiblePieceTemplate;

        Piece eliminationPieceTemplate = selectEliminationPieceTemplate(currShape.current_templates, currShape, pieceNumThreshold, noPossiblePieceTemplate);
        cout << endl << eliminationPieceTemplate.printPixels() << endl;

        if (noPossiblePieceTemplate) 
        {
            cout << "no Possible Piece Template. " << endl;
            break;
        }

        vector<Piece> eliminationPieces = getTemplateInstances(eliminationPieceTemplate, currShape);

        cout << "eliminationPieces: " << eliminationPieces.size() << endl;

        for (auto & piece : eliminationPieces)
        {
            cout << endl << piece.printPixels() << endl;
        }

        int N, K;
        vector<int> eliminationPieceIDList;
    
        Shape subShape = getSubShape(eliminationPieces, currShape, N, K, eliminationPieceIDList);
        subShape.updateTemplates(false);
        updateAccessibilityAndBlockability(subShape);

        // cout << "Subshape: N = " << N << ", K = " << K << endl;  
        // cout << endl << subShape.printAllPixels() << endl;

        int original_K = currShape.current_templates.size();

        // cout << "before delete K: " << currShape.current_templates.size() << endl;
        // cout << "before delete N: " << currShape.pieces.size() << endl;

        sort(eliminationPieceIDList.begin(), eliminationPieceIDList.end());
        for (int i = eliminationPieceIDList.size() - 1; i >= 0; --i)
        {
            currShape.pieces.erase(currShape.pieces.begin() + eliminationPieceIDList[i]);
        }

        currShape.updateTemplates(false);

        // cout << "after delete K: " << currShape.current_templates.size() << endl;
        // cout << "after delete N: " << currShape.pieces.size() << endl;

        double timeLimitForEachSubshape = 5.0;
        double runningTime;
        int bestScore;

        auto start = high_resolution_clock::now();

        int seed = 0;
        while (1)
        {
            bool isFound = false;
            Shape newSubShape = generatePieces(subShape, N, K, 15, 5, false, " ",
                                           seed, 0.3, 3.0, minTileSize, maxTileSize, 2, true, 10, runningTime, bestScore, currShape.current_templates, isFound, false);

            if (isFound)
            {
                vector<Piece> updatedTemplates = currShape.current_templates;

                for (auto & piece : newSubShape.current_templates)
                {
                    if (!isPieceinTemplates(piece, updatedTemplates))
                    {
                        updatedTemplates.push_back(piece);
                    }
                }

                cout << "original_K: " << original_K << endl;
                cout << "after_K: " << updatedTemplates.size() << endl;

                if (updatedTemplates.size() < original_K)
                {
                    cout << "successfully reduce " << original_K - updatedTemplates.size() << " templates. " << endl;

                    for (auto & piece : newSubShape.pieces)
                    {
                        currShape.pieces.push_back(piece);
                    }

                    newShape = currShape;
                    newShape.updateTemplates(false);
                    newShape.assignPieceID();
                    newShape.assignPieceClusterID();

                    isSuccess = true;

                    break;
                }
            }

            auto time_now = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(time_now - start);
            double time_elapsed = (duration.count() / (double)1000) / 60;

            if (time_elapsed > timeLimitForEachSubshape)
            {
                cout << "unsuccessfully to reduce K. " << endl;
                break;
            }

            ++seed;         
        }

        if (isSuccess)
            break;
    }

    return newShape;
}

Shape optimizePiecesSizeLocally(Shape & shape, int minTileSize, int maxTileSize)
{
    Shape newShape = shape;

    cout << "There are " <<  newShape.getSmallSizePiecesNum(minTileSize) << " small pieces and " << newShape.getLargeSizePiecesNum(maxTileSize) << " large pieces." << endl;

    for(int k = 0; k < 1; ++k)
    {
        Shape currShape = newShape;

        vector<Piece> eliminationPieces;
        for (int i = 0; i < currShape.pieces.size(); ++i)
        {
            if (currShape.pieces[i].pixels.size() < minTileSize or currShape.pieces[i].pixels.size() > maxTileSize)
            {
                eliminationPieces.push_back(currShape.pieces[i]);
                break;
            }
        }

        int N, K;
        vector<int> eliminationPieceIDList;
    
        Shape subShape = getSubShape(eliminationPieces, currShape, N, K, eliminationPieceIDList);
        subShape.updateTemplates(false);
        updateAccessibilityAndBlockability(subShape);

        int originalInvalidPieceNum = subShape.getSmallSizePiecesNum(minTileSize) + subShape.getLargeSizePiecesNum(maxTileSize);

        sort(eliminationPieceIDList.begin(), eliminationPieceIDList.end());
        for (int i = eliminationPieceIDList.size() - 1; i >= 0; --i)
        {
            currShape.pieces.erase(currShape.pieces.begin() + eliminationPieceIDList[i]);
        }

        currShape.updateTemplates(false);

        double timeLimitForEachSubshape = 1.0;
        double runningTime;
        int bestScore;

        auto start = high_resolution_clock::now();

        int seed = 0;
        while (1)
        {
            bool isFound = false;
            Shape newSubShape = generatePieces(subShape, N, N, 10, 1, false, " ",
                                           seed, 0.3, 3.0, minTileSize, maxTileSize, 2, true, 5, runningTime, bestScore, currShape.current_templates, isFound, false);

            if (isFound)
            {
                for (auto & piece : newSubShape.pieces)
                {
                    currShape.pieces.push_back(piece);
                }

                cout << "successfully reduce " << originalInvalidPieceNum << " invalid pieces. " << endl;

                newShape = currShape;
                newShape.updateTemplates(false);
                newShape.assignPieceID();
                newShape.assignPieceClusterID();

                break;
            }

            auto time_now = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(time_now - start);
            double time_elapsed = (duration.count() / (double)1000) / 60;

            if (time_elapsed > timeLimitForEachSubshape)
            {
                cout << "unsuccessfuly to reduce K. " << endl;
                break;
            }

            ++seed;         
        }

    }

    cout << "There are " <<  newShape.getSmallSizePiecesNum(minTileSize) << " small pieces and " << newShape.getLargeSizePiecesNum(maxTileSize) << " large pieces." << endl;

    return newShape;
}

Shape getSubShape(vector<Piece> & eliminationPieces, Shape & shape, int & pieceNum, int & K, vector<int> & eliminationPieceIDList)
{
    Shape newShape = shape;
    newShape.all_pixels.clear();
    newShape.pieces.clear();

    set<int> pieceIDList;
    vector<Piece> currentTemplates;

    for (int i = 0; i < eliminationPieces.size(); ++i)
    {
        vector<Piece> pieceList = getNeighbouringAndCurrentPieces(eliminationPieces[i], shape);
        for (auto & piece : pieceList)
        {
            if (!isPieceinTemplates(piece, currentTemplates))
            {
                currentTemplates.push_back(piece);
            }

            for (auto & px : piece.pixels)
            {
                newShape.all_pixels.insert(px);
            }

            pieceIDList.insert(shape.getPieceID(piece));
        }
    }

    K = currentTemplates.size();

    for (auto & i : pieceIDList)
    {
        eliminationPieceIDList.push_back(i);
    }

    pieceNum = eliminationPieceIDList.size();

    return newShape;
}

Shape generatePieces(Shape & shape, int N, int K, int G, int T, bool autoSave, string pureFileName, 
                     int seedNum, float minAvgSizeRatio, float maxAvgSizeRatio, int minTileSize, int maxTileSize,
                     int minSeedDist, bool isUniformDistribution, int shapeCandisNum,
                     double & runningTime, int & bestScore, vector<Piece> & existingTemplates, bool & isFound, bool isLooseSizeRequirement)
{
    auto start = high_resolution_clock::now();

    int best_score = INT_MAX;
    srand(seedNum);

    Shape best_solution;    
    Shape input = shape;
    input.clear();

    updateAccessibility(input);
    updateBlockability(input);

    while(1) 
    {
        Shape seed = generateNSeedPieces(N, input, minSeedDist, isUniformDistribution);

        cout << "Generate " << seed.pieces.size() << " seeds. " << endl;

        seed.maxAvgSizeRatio = maxAvgSizeRatio;
        seed.minAvgSizeRatio = minAvgSizeRatio;
        seed.minTileSize = minTileSize;
        seed.maxTileSize = maxTileSize;

        auto time_now = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_now - start);
        double time_elapsed = (duration.count() / (double)1000) / 60;

        if (time_elapsed > T) {
            cout << "Exceeded runtime....exiting.\n";
            break;
        }
                
        // repeat growing process G times and store each solution for testing
        for (int j = 0; j < G; j++) {
            Shape newShape = seed;
            bool isDisconnected = false;

            vector<Shape> growingStates;
            growingStates.push_back(input);
            growingStates.push_back(newShape);

            while (true) {        

                newShape = growPiecesByOnePixel(newShape, shapeCandisNum, isDisconnected, existingTemplates, minTileSize, maxTileSize, isLooseSizeRequirement, growingStates);

                if (isDisconnected)
                {
                    cout << "Find unreachable pixels or fail to satisfy tile size requirement. " << endl;
                    break;
                }

                cout << "Grew all pieces by one pixel. " << newShape.getNumUnfilledPx() << " pixels remaining. \n";

                // complete
                if (newShape.getNumUnfilledPx() <= 0) 
                { 
                    if (newShape.score <= best_score && newShape.getNumUnfilledPx() == 0 && newShape.score > 0) 
                    {
                        cout << "Find a N = " << newShape.pieces.size() << " K = " << newShape.score << " result. " << endl;

                        best_score = newShape.score;
                        best_solution = newShape;

                        auto now = high_resolution_clock::now();
                        auto dura = duration_cast<milliseconds>(time_now - start);
                        double elapsed = (duration.count() / (double)1000) / 60;

                        best_solution.assignPieceID();
                        best_solution.assignPieceClusterID();

                        // save current best result
                        if (autoSave)
                        {
                            // string savingPath =  "./Result/auto_save/";
                            string savingPath =  "../";

                            string pureFolderName = best_solution.getOutputFolderPath(savingPath + pureFileName) + "_size" + to_string(best_solution.getSmallestPieceSize()) + "To" + to_string(best_solution.getLargestPieceSize()) + "_" + to_string(elapsed) + "min";
                            saveShape2File(best_solution, pureFolderName, best_solution.getOutputFileFullName(pureFolderName + "/" + pureFileName), elapsed);

                            string pureGrowingStateSavingFolderName = pureFolderName + "/enlarging_states";
                            for (int stateID = growingStates.size() - 1; stateID >= 0; --stateID)
                            {
                                if (growingStates[stateID].getNumUnfilledPx() > 0)
                                {
                                    vector<Pixel> remaining_pixels;
                                    for (const auto &px : growingStates[stateID].all_pixels) 
                                        if (!isPixelInPieces(px, growingStates[stateID].pieces)) 
                                            remaining_pixels.push_back(px);
                                    
                                    if (!remaining_pixels.empty())
                                        growingStates[stateID].pieces.push_back(Piece(remaining_pixels, Color{150,150,150}));
                                }

                                // saveShape2File(growingStates[stateID], pureGrowingStateSavingFolderName, pureGrowingStateSavingFolderName + "/" + to_string(stateID) + ".puz", elapsed);
                            }
                        }

                    }

                    if (best_score <= K) {
                        isFound = true;
                    }

                    break;
                }
            }

            // cout << "Time Elapsed: " << time_elapsed  << " minutes.\n";

            if (isFound)
                break;    
        }

        if (isFound)
            break;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    runningTime = duration.count()/(double)1000/60;
    bestScore = best_solution.score;

    cout << endl << "Best solution with K = " << bestScore << "\n";
    cout << "Total runtime: " << runningTime << "min.\n";

    return best_solution;
}

Shape generatePiecesWithGivenSeeds(Shape & shape, int N, int K, int G, int T, bool autoSave, string pureFileName, 
                     int seedNum, float minAvgSizeRatio, float maxAvgSizeRatio, int minTileSize, int maxTileSize,
                     int minSeedDist, bool isUniformDistribution, int shapeCandisNum,
                     double & runningTime, int & bestScore, vector<Piece> & existingTemplates, bool & isFound, bool isLooseSizeRequirement)
{
    auto start = high_resolution_clock::now();

    int best_score = INT_MAX;
    srand(seedNum);

    Shape best_solution;    
    Shape input = shape;
    input.clear();

    while(1) 
    {
        Shape seed = input;

        seed.maxAvgSizeRatio = maxAvgSizeRatio;
        seed.minAvgSizeRatio = minAvgSizeRatio;
        seed.minTileSize = minTileSize;
        seed.maxTileSize = maxTileSize;

        auto time_now = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_now - start);
        double time_elapsed = (duration.count() / (double)1000) / 60;

        if (time_elapsed > T) {
            cout << "Exceeded runtime....exiting.\n";
            break;
        }
                
        // repeat growing process G times and store each solution for testing
        for (int j = 0; j < G; j++) {
            Shape newShape = seed;
            bool isDisconnected = false;
            vector<Shape> growingStates;

            while (true) {        

                newShape = growPiecesByOnePixel(newShape, shapeCandisNum, isDisconnected, existingTemplates, minTileSize, maxTileSize, isLooseSizeRequirement, growingStates);

                if (isDisconnected)
                {
                    cout << "Find unreachable pixels or fail to satisfy tile size requirement. " << endl;
                    break;
                }

                cout << "Grew all pieces by one pixel. " << newShape.getNumUnfilledPx() << " pixels remaining. \n";

                if (newShape.getNumUnfilledPx() <= 0) { 
                    // complete

                    if (newShape.getLargestPieceSize() > maxTileSize or newShape.getSmallestPieceSize() < minTileSize)
                    {
                        // cout << "Fail to satisfy tile size requirement." << endl;
                        break;
                    }

                    if (newShape.score <= best_score && newShape.getNumUnfilledPx() == 0 && newShape.score > 0)  {

                        cout << "Find a N = " << newShape.pieces.size() << " K = " << newShape.score << " result. " << endl;

                        best_score = newShape.score;
                        best_solution = newShape;

                        auto now = high_resolution_clock::now();
                        auto dura = duration_cast<milliseconds>(time_now - start);
                        double elapsed = (duration.count() / (double)1000) / 60;

                        // save current best result
                        if (autoSave)
                        {
                            string savingPath =  "../Result/auto_save/";
                            string pureFolderName = best_solution.getOutputFolderPath(savingPath + pureFileName) + "_size" + to_string(best_solution.getSmallestPieceSize()) + "To" + to_string(best_solution.getLargestPieceSize()) + "_" + to_string(elapsed) + "min";
                            saveShape2File(best_solution, pureFolderName, 
                            best_solution.getOutputFileFullName(pureFolderName + "/" + pureFileName), elapsed);
                        } 

                    }

                    if (best_score <= K) {
                        isFound = true;
                    }
                    break;
                }
            }

            // cout << "Time Elapsed: " << time_elapsed  << " minutes.\n";

            if (isFound)
                break;    
        }

        if (isFound)
            break;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    runningTime = duration.count()/(double)1000/60;
    bestScore = best_solution.score;

    cout << endl << "Best solution with K = " << bestScore << "\n";
    cout << "Total runtime: " << runningTime << "min.\n";

    return best_solution;   
}

// Assumes N < the number of pixels in that piece
Shape generateNSeedPieces(int N, Shape & shape, int minSeedDist, bool isUniformDistribution) {

    vector<float> possibList = shape.accessibilityList;
    vector<Piece> pieces;
    vector<Pixel> seed_pixels;
    vector<Pixel> pixelsInShape = vector<Pixel>(shape.all_pixels.begin(), shape.all_pixels.end());
    int count = 0;

    while(1)
    {
        seed_pixels.clear();

        possibList = shape.accessibilityList;

        if (isUniformDistribution)
        {
            for (int i = 0; i < possibList.size(); ++i)
            {
                if (possibList[i] > 0)
                    possibList[i] = 1;
            }
        }

        while(1)
        {
            int random_index = GetRandomObjIndex(possibList, 3.0);
            Pixel seed_pixel = Pixel(int(random_index / shape.volCol), random_index % shape.volCol);

            cout << "currSelectedSeedID: " << random_index << endl;
            cout << "(x,y): " <<  seed_pixel.x << " " << seed_pixel.y << endl << endl;

            seed_pixels.push_back(seed_pixel);
            count ++;

            if (count == N)
                break;

            vector<Pixel> visitedPixels;
            vector<Pixel> nextVisits;
            vector< pair<Pixel,int> > nextVisitPixels;

            nextVisitPixels.push_back(pair(seed_pixel,0));
            nextVisits.push_back(seed_pixel);

            while (!nextVisitPixels.empty())
            {
                int currDepth = nextVisitPixels[0].second;
                Pixel currPixel = nextVisits[0];

                nextVisits.erase(nextVisits.begin());
                nextVisitPixels.erase(nextVisitPixels.begin());
                visitedPixels.push_back(currPixel);

                cout << "currPixelID: " << currPixel.y * shape.volCol + currPixel.x << endl;
                cout << "(x,y): " <<  currPixel.x << " " << currPixel.y << endl << endl;

                possibList[currPixel.x * shape.volCol + currPixel.y] = 0;

                // vector<Pixel> newNeigbouringPixels = getValidNeigbours(currPixel, shape);
                vector<Pixel> newNeigbouringPixels = getValidCircleNeigbours(currPixel, shape);

                for (int i = 0; i < newNeigbouringPixels.size(); ++i)
                {
                    if (currDepth < minSeedDist - 1)
                    {
                        if(!isPixelInPixelList(newNeigbouringPixels[i], nextVisits) and !isPixelInPixelList(newNeigbouringPixels[i], visitedPixels))
                        {
                            nextVisitPixels.push_back(pair(newNeigbouringPixels[i], currDepth + 1));
                            nextVisits.push_back(newNeigbouringPixels[i]);
                        }
                    }
                } 
            }

            if (*max_element(possibList.begin(), possibList.end()) == 0)
                break;
        }

        if (seed_pixels.size() == N)
            break;
    }

    // cout << "Find " << count << " seeds with minimal dist " << minSeedDist << endl;

    for (const auto& px : seed_pixels) 
    {
        pieces.push_back(Piece(px));
    }

    Shape seed = shape;
    
    for (auto piece : pieces) 
    {
        seed.addPiece(piece);
    }

    seed.assignColors();
    seed.valid = true;

    Piece newPiece = pieces[0];

    seed.current_templates.clear();
    seed.current_templates.push_back(newPiece);

    updateAccessibilityAndBlockability(seed);

    return seed;
}

// grows each piece in the shape by one pixel according to the algorithm
// if use_access_val is true, it will pick pixels based on accessibility value, otherwise: pick randomly
// returns the result as a new shape
Shape growAllPiecesByOnePixel(const Shape& shape, bool use_access_val, vector<Shape> & shape_state) {

    Shape newShape = shape; //copy shape
    newShape.current_templates = vector<Piece>(); // clear current templates

    // 1) randomly choose a starting piece (that is not blocked)
    int random_piece_idx = getRandomNumber(0, newShape.pieces.size() - 1);
    while (isBlocked(newShape.pieces[random_piece_idx], newShape)) {
        random_piece_idx = getRandomNumber(0, newShape.pieces.size() - 1);
    }
    
    // 2) grow that piece by one pixel
    Piece newPiece = growPiece(newShape, random_piece_idx, use_access_val);
    newShape.pieces[random_piece_idx] = newPiece;
    newShape.current_templates.push_back(newPiece);

    // 3) grow remaining pieces
    for (int i = 0; i < newShape.pieces.size(); i++) {
        if (i != random_piece_idx) {
            if (isBlocked(newShape.pieces[i], newShape)) {
                if (isPieceinTemplates(newShape.pieces[i], newShape.blocked_templates)) {
                    // already a blocked template, leave this piece be and continue
                    continue;
                } else {
                    double target_template_size = ((double)newShape.all_pixels.size() / (double)newShape.pieces.size());
                    if ((double)newShape.pieces[i].pixels.size() >= 0.8*target_template_size) { // if this piece is already big enough, add it to the blocked template list and continue
                        // add to blocked template list in this shape, and continue
                        newShape.blocked_templates.push_back(newShape.pieces[i]);
                        continue;
                    } else { 

                        // remove backtracking temporarily 

                        // backtrack randomly 1-3 rounds
                        // int backtrack_amt = getRandomNumber(1, 5);
                        // if (backtrack_amt >= shape_state.size()) {
                        //     backtrack_amt = shape_state.size() - 1;
                        // }
                        // cout << "Backtacking " << backtrack_amt << " rounds....\n";
                        // while (backtrack_amt != 0) {
                        //     shape_state.pop_back();
                        //     backtrack_amt--;
                        // }

                        // int backtrack_count = newShape.backtrack_count;
                        // newShape = shape_state[shape_state.size() - 1];
                        // newShape.backtrack_count = backtrack_count + 1;

                        // shape_state.pop_back();
                        return newShape;

                    }
                }
            }
            
            vector<Pixel> next_possible_pixels = getPossiblePixelsToGrow(newShape.pieces[i], newShape);
            if (next_possible_pixels.size() == 0) {
                // no possible pixels that match existing templates, so we create a new template
                Piece newPiece = growPiece(newShape, i, use_access_val);
                newShape.pieces[i] = newPiece;
                newShape.current_templates.push_back(newPiece);
            } else {
                // randomly choose one of the possible pixels to grow
                int rand_idx = getRandomNumber(0, next_possible_pixels.size() - 1);
                Pixel next_px = next_possible_pixels[rand_idx];
                if (next_px.x == -1) { // this represents a blocked template match: no need to grow the piece
                    continue;
                } else {
                    newShape.pieces[i].addPixel(next_px);
                }
                
            }
        }
    }

    newShape.score = newShape.current_templates.size() + newShape.blocked_templates.size();
    return newShape;
}

Shape growPiecesByOnePixel(Shape & shape, int backtrackNum, bool & isDisconnected, vector<Piece> & existingTemplates, int minTileSize, int maxTileSize, bool isLooseSizeRequirement, vector<Shape> & growingStates)
{
    vector<Shape> shapeCandis;
    vector<vector<Shape>> intermediateShapeCandis;
    vector<float> distinctPieceNumList;
    vector<float> growingSpaceEvaluationList;
    vector<float> possibList;

    // growing pieces in the same piece group by adding one pixel
    for(int k = 0; k < backtrackNum; ++k)
    {
        Shape newShape = shape; // new shape to grow

        // 1) select the guiding piece class template according to:
        // i. the number of instances
        // ii. piece size
        Piece guidingPieceClassTemplate = selectGuidingPieceClassTemplate(shape.current_templates, newShape, 5.0, isDisconnected);

        if (isDisconnected)
        {
            // updateTemplates(newShape);
            // cout << "Find one or more unreachable pixels." << endl;
            return shape;
        }

        // 2) select a piece in guiding piece class according to:
        // i. the less growing space

        vector<Piece> unblocked_instances = getUnblockedTemplateInstances(guidingPieceClassTemplate, newShape);
        // cout << "unblocked_instances: " << unblocked_instances.size() << endl;

        Piece firstSelectedPiece = selectFirstGrowingPieceInGuidingPieceClass(unblocked_instances, newShape, 3.0);
        // cout << "find first selected piece. " << endl;

        // 3) grow the first selected piece by adding one pixel according to: 
        // i. match one of the existing current templates
        // ii. be the subset of one of the existing current templates
        // iii. accessibility

        vector<Piece> combinedTemplates;
        for (int i = 0; i < shape.current_templates.size(); ++i)
        {
            combinedTemplates.push_back(shape.current_templates[i]);
        }

        for (int i = 0; i < existingTemplates.size(); ++i)
        {
            combinedTemplates.push_back(existingTemplates[i]);
        }

        Piece firstGrowedPiece = growFirstSelectedPiece(newShape, combinedTemplates, firstSelectedPiece, 3.0);

        int firstGrowedPieceID = newShape.getPieceID(firstSelectedPiece);
        newShape.pieces[firstGrowedPieceID] = firstGrowedPiece;

        // 4) grow the remaining piece to match this piece

        // for each remaining piece that may be grow to the same shape:
        //      if blocked:
        //          add to blocked_templates<>;
        //      else:
        //          if can match the growed template:
        //              grow to match;
        //          else:
        //              leave it to be to the next iteration;
        //              isEraseSelectedTemplate = false;

        vector<Pixel> addedPixels;
        vector< pair<int, vector<Pixel>> > piecePairs;

        for (int i = 0; i < newShape.pieces.size(); ++i)
        {
            if (i == firstGrowedPieceID)
            {
                continue;
            }

            else
            {
                // if the pieces[i] and firstSelectedPiece has the same shape 
                if (isSameShape(newShape.pieces[i], firstSelectedPiece))
                {
                    if (isBlocked(newShape.pieces[i], newShape))
                    {
                        continue;
                    }
                    else
                    {
                        vector<Pixel> currPossibMatchPixels = getPossiblePixelsGrow2Match(newShape, newShape.pieces[i], firstGrowedPiece);

                        if (!currPossibMatchPixels.empty())
                        {
                            piecePairs.push_back(pair(i, currPossibMatchPixels));
                        }
                    }

                    continue;
                }
            }
        }

        std::sort(piecePairs.begin(), piecePairs.end(),
              [](pair<int, vector<Pixel>> & a, pair<int, vector<Pixel>> & b) {
                  return a.second.size() < b.second.size();
              });

        while (!piecePairs.empty())
        {
            int currID = piecePairs[0].first;
            Piece currPiece = newShape.pieces[currID];
            vector<Pixel> currPossibPixels = piecePairs[0].second;
            vector<float> currPossibList;

            for (int i = 0; i < currPossibPixels.size(); ++i)
            {
                if (!isPixelInPixelList(currPossibPixels[i], addedPixels))
                {
                    int currPixelID = newShape.getPixelID(currPossibPixels[i]);
                    currPossibList.push_back(round(5.0 * newShape.blockabilityList[currPixelID]));
                }
                else
                {
                    currPossibList.push_back(0);
                }
            }

            if (*max_element(currPossibList.begin(), currPossibList.end()) > 0)
            {
                int growPixelID = GetRandomObjIndex(currPossibList, 3.0);
                Pixel growPixel = currPossibPixels[growPixelID];

                currPiece.addPixel(growPixel);
                newShape.pieces[currID] = currPiece;

                addedPixels.push_back(growPixel);
            }

            piecePairs.erase(piecePairs.begin());
        }

        // 5) fill the disconnected pixels
        vector<Shape> currentCandis;
        currentCandis.push_back(newShape);

        int prevUnfilledPxNum = newShape.getNumUnfilledPx();
        fillDisconnectedPixels(newShape);
        int afterUnfilledPxNum = newShape.getNumUnfilledPx();

        if (afterUnfilledPxNum < prevUnfilledPxNum)
        {
            currentCandis.push_back(newShape);
        }

        // currentCandis.push_back(newShape);

        // cout << "Fill " << prevUnfilledPxNum - afterUnfilledPxNum << " dangling pixels. " << newShape.getNumUnfilledPx() << " pixels remaining\n" ;

        // 6) update templates, accessibility and score
        newShape.updateTemplates(isLooseSizeRequirement);
        updateBlockability(newShape);

        // 7) add current shape to the candidate list if valid
        bool isValid = true;
        int currMaxSize = maxTileSize;
        int currMinSize = minTileSize;

        if (newShape.getNumUnfilledPx() <= 10 and isLooseSizeRequirement)
        {
            currMaxSize += 1;
            currMinSize -= 1;
        }

        for (int i = 0; i < newShape.pieces.size(); ++i)
        {
            if (newShape.pieces[i].pixels.size() > currMaxSize)
            {
                isValid = false;
                break;
            }

            if (isBlocked(newShape.pieces[i], newShape) and newShape.pieces[i].pixels.size() < currMinSize)
            {
                isValid = false;
                break;
            }
        }

        if (isValid)
        {
            shapeCandis.push_back(newShape);
            intermediateShapeCandis.push_back(currentCandis); 
        }

        if (shapeCandis.size() >= 2)
            break;
    }

    if (shapeCandis.empty())
    {
        isDisconnected = true;
        return shape;
    }
    else
    {
        Shape finalShape = shapeCandis[0];
        int selectedID = 0;

        for (int i = 1; i < shapeCandis.size(); ++i)
        {
            if (shapeCandis[i].current_templates.size() < finalShape.current_templates.size())
            {
                finalShape =  shapeCandis[i];
                selectedID = i;
            }

            // if (shapeCandis[i].getAverageEnlargeability() > finalShape.getAverageEnlargeability())
            // {
            //     finalShape =  shapeCandis[i];
            //     selectedID = i;
            // }
        }

        // for (int i = 1; i < shapeCandis.size(); ++i)
        // {
        //     if (shapeCandis[i].getInstanceNumDeviation() < finalShape.getInstanceNumDeviation())
        //        finalShape =  shapeCandis[i];
        // }

        for (int i = 0; i < intermediateShapeCandis[selectedID].size(); ++i)
        {
            growingStates.push_back(intermediateShapeCandis[selectedID][i]);
        }

        // cout << "enlargeability: " << finalShape.getAverageEnlargeability() << endl;

        return finalShape;
    }
}


// given a piece in a shape, determines which pixels this piece can grow into next that match the set of existing templates in the shape
vector<Pixel> getPossiblePixelsToGrow(Piece & piece, Shape & shape) {
    vector<Pixel> possible_pixels;

    if (isPieceinTemplates(piece, shape.blocked_templates)) {
        possible_pixels.push_back(Pixel(-1, -1));
    }

    vector<Pixel> neighbouring_px = getNeighbouringPixels(piece, shape);
    for (auto& px : neighbouring_px) {
        Piece newPiece = piece; //copy piece
        newPiece.addPixel(px);

        if (isPieceinTemplates(newPiece, shape.current_templates)) {
            possible_pixels.push_back(px);
        }
    }

    return possible_pixels;

}

// given a piece in a shape, grows it by one pixel
// returns the same piece if this is not possible
Piece growPiece(const Shape& shape, int piece_index, bool use_access_val) {
    // get pixels neighbouring the piece
    vector<Pixel> neighbouring_pixels = getNeighbouringPixels(shape.pieces[piece_index], shape);
    if (neighbouring_pixels.size() == 0) {
        return shape.pieces[piece_index];
    }
    else {
        // choose the next pixel
        Pixel next_px(0,0);
        if (use_access_val) {
            next_px = getRandomPixelBasedOnAccessibilityValue(neighbouring_pixels, shape);
        }
        else {
            next_px = getRandomPixel(neighbouring_pixels);
        }
        
        // grow piece with that pixel
        Piece newPiece = shape.pieces[piece_index];
        newPiece.addPixel(next_px);
        return newPiece;
    }
}

// given a piece in a shape, grows it by one pixel
// returns the same piece if this is not possible
Piece growPiece(const Shape& shape, Piece & piece) {
    // get pixels neighbouring the piece
    vector<Pixel> neighbouring_pixels = getNeighbouringPixels(piece, shape);
    if (neighbouring_pixels.size() == 0) {
        return piece;
    }
    else {
        // choose the next pixel
        Pixel next_px(0,0);

        next_px = getRandomPixelBasedOnAccessibilityValue(neighbouring_pixels, shape);
        
        // grow piece with that pixel
        Piece newPiece = piece;
        newPiece.addPixel(next_px);
        return newPiece;
    }
}

Piece growFirstSelectedPiece(Shape& shape, vector<Piece> & current_templates , Piece & piece, float alpha)
{
    // get pixels neighbouring the piece
    vector<Pixel> neighbouring_pixels = getNeighbouringPixels(piece, shape);

    if (neighbouring_pixels.size() == 0) {
        return piece;
    }
    else
    {
        vector<float> possibList;
        vector<Piece> candidate_pieces;

        for (auto & pixel : neighbouring_pixels)
        {
            Pixel currPixel = pixel;
            int currPixelID = shape.getPixelID(currPixel);

            // grow piece with that pixel
            Piece newPiece = piece;
            newPiece.addPixel(currPixel);

            if (isPieceinTemplates(newPiece, current_templates))
            {
                possibList.push_back(30 + round(shape.blockabilityList[currPixelID]));
            }
            else if (isPieceIsSubsetShapeinTemplates(newPiece, current_templates))
            {
                // cout << "PieceIsSubsetShapeinTemplates. " << endl;
                possibList.push_back(30 + round(shape.blockabilityList[currPixelID]));
            }
            else
            {
                possibList.push_back(round(shape.blockabilityList[currPixelID]));
            }

            candidate_pieces.push_back(newPiece);
        }

        return candidate_pieces[GetRandomObjIndex(possibList, alpha)];
    }
}

vector<Pixel> getPossiblePixelsGrow2Match(Shape & shape, Piece & piece, Piece & piece_template)
{
    vector<Pixel> possiblePixels;

    // get pixels neighbouring the piece
    vector<Pixel> neighbouring_pixels = getNeighbouringPixels(piece, shape);

    if (neighbouring_pixels.size() == 0) {
        return possiblePixels;
    }

    else
    {
        for (auto & pixel : neighbouring_pixels)
        {
            Piece newPiece = piece;
            newPiece.addPixel(pixel);
            int currPixelID = shape.getPixelID(pixel);

            if (isSameShape(newPiece, piece_template))
            {
                possiblePixels.push_back(pixel);
            }
        }

        return possiblePixels;
    }
}

Piece growPiece2Match(Shape & shape, Piece & piece, Piece & piece_template, bool & isMatched, float alpha)
{
    // get pixels neighbouring the piece
    vector<Pixel> neighbouring_pixels = getNeighbouringPixels(piece, shape);

    if (neighbouring_pixels.size() == 0) {
        isMatched = false;
        return piece;
    }

    else
    {
        isMatched = false;
        vector<Piece> candidates;
        vector<float> possibList;

        for (auto & pixel : neighbouring_pixels)
        {
            Piece newPiece = piece;
            newPiece.addPixel(pixel);
            int currPixelID = shape.getPixelID(pixel);

            if (isSameShape(newPiece, piece_template))
            {
                isMatched = true;

                candidates.push_back(newPiece);
                // possibList.push_back(round(0.6 * shape.accessibilityList[currPixelID] + 0.4 * shape.blockabilityList[currPixelID]));
                possibList.push_back(round(shape.blockabilityList[currPixelID]));
            }
        }

        if (!isMatched)
        {
            return piece;
        }
        else
        {
            return candidates[GetRandomObjIndex(possibList, alpha)];
        }
    }
}

// tries to grow the a given piece in a shape in the same orientation as a given template
// returns the same piece if this is not possible
Piece growPieceWithSameShape(const Shape& shape, int piece_index, vector<Pixel> neighbouring_pixels, Piece& templatePiece, bool use_access_val) {
    if (neighbouring_pixels.size() == 0) {
        return shape.pieces[piece_index];
    } else {
        // choose the next pixel 
        Pixel next_px(0, 0);
        if (use_access_val) {
            next_px = getRandomPixelBasedOnAccessibilityValue(neighbouring_pixels, shape);
        } else {
            next_px = getRandomPixel(neighbouring_pixels);
        }
        // grow piece with that pixel
        Piece newPiece = shape.pieces[piece_index];
        newPiece.addPixel(next_px);
        // check shape
        if (isSameShape(newPiece, templatePiece)) {
            return newPiece;
        } else {
            neighbouring_pixels.erase(remove(neighbouring_pixels.begin(), neighbouring_pixels.end(), next_px), neighbouring_pixels.end());
            return growPieceWithSameShape(shape, piece_index, neighbouring_pixels, templatePiece, use_access_val);
        }
    }
}

Piece selectGuidingPieceClassTemplate(vector<Piece> & piece_templates, Shape & shape, float alpha, bool & isDisconnected)
{
    vector<float> possibList;
    vector<float> instanceNumList;
    vector<float> pieceSizeList;

    for (auto & piece : piece_templates)
    {
        pieceSizeList.push_back(piece.pixels.size());
        instanceNumList.push_back(getCurrentTemplateInstanceNum(piece, shape));
    }

    normalizeVector(pieceSizeList);
    normalizeVector(instanceNumList);

    int i = 0;
    for (auto & piece : piece_templates)
    {
        if (isPieceinTemplates(piece, shape.blocked_templates))
        {
            possibList.push_back(0);
            ++i;
            continue;
        }

        float currPossib = round( 9.0 * (1.0 - pieceSizeList[i]) ) + round( 1.0 * instanceNumList[i]);
        // float currPossib = round( 10.0 * (1.0 - pieceSizeList[i]) ) + round(0.0 * instanceNumList[i]);

        if (currPossib < 1.0) 
            currPossib = 1.0;

        possibList.push_back(currPossib);

        ++i;
    }

    if (*std::max_element(possibList.begin(), possibList.end()) == 0)
    {
        isDisconnected = true;
    }
    else
    {
        isDisconnected = false;
    }

    Piece newPiece = piece_templates[GetRandomObjIndex(possibList, alpha)];

    return newPiece;
}

Piece selectFirstGrowingPieceInGuidingPieceClass(vector<Piece> & guiding_piece_instances, Shape & shape, float alpha)
{
    vector<float> possibList;

    for (auto & piece : guiding_piece_instances)
    {
        float count = 0;

        vector<Pixel> neighbourPixels = getNeighbouringPixels(piece, shape);

        for (auto & neighourPixel : neighbourPixels)
        {
            int currPixelID = shape.getPixelID(neighourPixel);

            // count += round(0.6 * shape.accessibilityList[currPixelID] + 0.4 * shape.blockabilityList[currPixelID]);
            count += round(shape.blockabilityList[currPixelID]);
        }

        possibList.push_back(count);
    }

    // cout << "possibList: " << endl;
    // for (int i = 0; i < possibList.size(); ++i)
    // {
    //     cout << possibList[i] << " " ;
    // }
    // cout << endl;

    Piece newPiece = guiding_piece_instances[GetRandomObjIndex(possibList, alpha)];

    return newPiece;
}

Piece selectEliminationPieceTemplate(vector<Piece> & piece_templates, Shape & shape, int pieceNumThreshold, bool & noPossiblePieceTemplate)
{
    vector<float> possibList;
    vector<float> instanceNumList;
    vector<float> pieceSizeList;
    vector<float> instanceNumList_normalized;

    for (auto & piece : piece_templates)
    {
        pieceSizeList.push_back(piece.pixels.size());
        instanceNumList.push_back(getCurrentTemplateInstanceNum(piece, shape));
    }

    cout << "instanceNumList: " << endl;
    for (int i = 0; i < instanceNumList.size(); ++i)
    {
        cout << instanceNumList[i] << " " ;
    }
    cout << endl;

    cout << "pieceNumThreshold: " << pieceNumThreshold << endl;

    instanceNumList_normalized = instanceNumList;

    normalizeVector(pieceSizeList);
    normalizeVector(instanceNumList_normalized);

    int i = 0;
    for (auto & piece : piece_templates)
    {
        if (instanceNumList[i] >= pieceNumThreshold)
        {
            possibList.push_back(0);
            ++i;
            continue;
        }
        else{
            float currPossib = round( 10.0 * (1.0 - instanceNumList_normalized[i]) );

            // float currPossib = round( 10.0 * pieceSizeList[i]) ;

            if (currPossib < 1.0) 
                currPossib = 1.0;

            possibList.push_back(currPossib);

            ++i;
        }
    }

    cout << "possibList: " << endl;
    for (int i = 0; i < possibList.size(); ++i)
    {
        cout << possibList[i] << " " ;
    }
    cout << endl;

    if (*std::max_element(possibList.begin(), possibList.end()) == 0)
    {
        noPossiblePieceTemplate = true;
    }
    else
    {
        noPossiblePieceTemplate = false;
    }

    Piece newPiece = piece_templates[GetRandomObjIndex(possibList, 3.0)];

    return newPiece;
}

void fillDisconnectedPixels(Shape & shape)
{
    vector<int> toAssignPieceIDs;
    vector<Pixel> alonePixels = getAlonePixels(shape, toAssignPieceIDs);

    // cout << "alone pixels number: " << alonePixels.size() << endl << endl;

    int i = 0;
    for (auto & id : toAssignPieceIDs)
    {
        shape.pieces[id].addPixel(alonePixels[i]);
        ++i;
    }
}

void updateTemplates(Shape & shape)
{
    shape.current_templates.clear();
    shape.blocked_templates.clear();

    for (int i = 0; i < shape.pieces.size(); ++i)
    {
        shape.pieces[i].isBlocked = isBlocked(shape.pieces[i], shape);

        if (!isPieceinTemplates(shape.pieces[i], shape.current_templates))
            shape.current_templates.push_back(shape.pieces[i]);
    }

    for (auto & tmp : shape.current_templates)
    {
        bool isBlockedTmp = true;
        vector<Piece> instances = getTemplateInstances(tmp, shape);

        for (auto & instance : instances)
        {
            if (!instance.isBlocked)
            {
                isBlockedTmp = false;
                break;
            }
        }

        if (isBlockedTmp)
            shape.blocked_templates.push_back(tmp);
    }
}

