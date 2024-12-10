#pragma once
#ifndef VIEWER_H
#define VIEWER_H
#include "Shape.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <vector>

using namespace std;

class Window {

public:
    int originX;
    int originY;
    int width;
    int height;
    string label;

    Shape solution;
    vector<Piece> pieces;
    vector<float> accessibilityList; // accessibility for unoccupied squares

    Window(int originX, int originY, int width, int height);
    void addPiece(const Piece &piece);
    void addShape(Shape &shape);
    void addLabel(const string &label);
    void clear();
};

class Viewer {

public:
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    std::vector<Window> windows;

    Viewer();
    void addWindow(Window window);
    void clearViewerMeshes();
    void updateWindows();
    void launch(int w, int h);

    void updateViewerBufferAndSave(vector<Piece> all_templates, Shape state, int pieceNum, string filepath, int stateID, vector<Eigen::RowVector3d> additionalColorList);
};



#endif