#pragma once
#ifndef IOUTILS_H
#define IOUTILS_H
#include "TriPatchSolver.h"
#include "TriTriSolver.h"
#include "PatchPatchSolver.h"
#include "BezierTri.h"
#include "BSplinePatch.h"
#include <fstream>
using namespace std;

BezierTri readTri(ifstream& file) {
    string header_token;
    file >> header_token;
    int degree = 3;

    size_t pos = header_token.find(':');
    if (pos != string::npos) {
        degree = stoi(header_token.substr(pos + 1));
    }

    int num_points = (degree + 1) * (degree + 2) / 2;

    Points pts;
    pts.reserve(num_points);
    double x, y, z;

    for (int i = 0; i < num_points; ++i) {
        if (file >> x >> y >> z) {
            pts.push_back(Vector3d(x, y, z));
        }
    }

    return BezierTri(pts, degree);
}

BSplinePatch readPatch(ifstream& file) {
    int du, dv;
    file >> du >> dv;

    string line;
    getline(file, line);

    while (getline(file, line)) {
        if (!line.empty() && line.find_first_not_of(" \r\n\t") != string::npos) break;
    }
    Knots uk;
    stringstream ss_u(line);
    double val;
    while (ss_u >> val) {
        uk.push_back(val);
    }

    while (getline(file, line)) {
        if (!line.empty() && line.find_first_not_of(" \r\n\t") != string::npos) break;
    }
    Knots vk;
    stringstream ss_v(line);
    while (ss_v >> val) {
        vk.push_back(val);
    }

    int nu = uk.size() - du - 1;
    int nv = vk.size() - dv - 1;

    Points pts(nu * nv);
    double x, y, z;

    for (int v = 0; v < nv; ++v) {
        for (int u = 0; u < nu; ++u) {
            if (file >> x >> y >> z) {
                pts[u * nv + v] = Vector3d(x, y, z);
            }
        }
    }

    return BSplinePatch(pts, uk, vk, du, dv);
}

std::pair<BezierTri, BezierTri> readTriTri(std::ifstream& file) {
    int degree;
    file >> degree;

    int num_points = (degree + 1) * (degree + 2) / 2;

    Points pts1, pts2;
    pts1.reserve(num_points);
    pts2.reserve(num_points);

    double x, y, z;

    for (int i = 0; i < num_points; ++i) {
        if (file >> x >> y >> z) {
            pts1.push_back(Vector3d(x, y, z));
        }
    }

    for (int i = 0; i < num_points; ++i) {
        if (file >> x >> y >> z) {
            pts2.push_back(Vector3d(x, y, z));
        }
    }

    return std::make_pair(BezierTri(pts1, degree), BezierTri(pts2, degree));
}

std::pair<BSplinePatch, BSplinePatch> readPatchPatch(std::ifstream& file) {
    auto read_single_patch = [&file]() -> BSplinePatch {
        int du, dv;
        file >> du >> dv;

        std::string line;
        std::getline(file, line);

        // 读 U knots
        while (std::getline(file, line)) {
            if (!line.empty() && line.find_first_not_of(" \r\n\t") != std::string::npos) break;
        }
        Knots uk;
        std::stringstream ss_u(line);
        double val;
        while (ss_u >> val) uk.push_back(val);

        // 读 V knots
        while (std::getline(file, line)) {
            if (!line.empty() && line.find_first_not_of(" \r\n\t") != std::string::npos) break;
        }
        Knots vk;
        std::stringstream ss_v(line);
        while (ss_v >> val) vk.push_back(val);

        // 动态计算控制点网格尺寸
        int nu = uk.size() - du - 1;
        int nv = vk.size() - dv - 1;

        Points pts(nu * nv);
        double x, y, z;

        // 列优先读入，行优先存储
        for (int v = 0; v < nv; ++v) {
            for (int u = 0; u < nu; ++u) {
                if (file >> x >> y >> z) {
                    pts[u * nv + v] = Vector3d(x, y, z);
                }
            }
        }

        return BSplinePatch(pts, uk, vk, du, dv);
        };
    BSplinePatch patch1 = read_single_patch();
    BSplinePatch patch2 = read_single_patch();

    return std::make_pair(patch1, patch2);
}

std::pair<BezierTri, BSplinePatch> readTriPatch(std::ifstream& file) {
    std::string line;

    int tri_degree = 3;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            std::stringstream ss(line);
            ss >> tri_degree;
            break;
        }
    }

    int num_tri_pts = (tri_degree + 1) * (tri_degree + 2) / 2;
    Points tri_pts;
    tri_pts.reserve(num_tri_pts);

    while (tri_pts.size() < num_tri_pts && std::getline(file, line)) {
        for (char& c : line) {
            if (c == ',') c = ' ';
        }

        std::stringstream ss(line);
        double x, y, z;
        if (ss >> x >> y >> z) {
            tri_pts.push_back(Vector3d(x, y, z));
        }
    }

    Points patch_pts;
    while (std::getline(file, line)) {
        for (char& c : line) {
            if (c == ',') c = ' ';
        }

        std::stringstream ss(line);
        double x, y, z;
        if (ss >> x >> y >> z) {
            patch_pts.push_back(Vector3d(x, y, z));
        }
    }

    int patch_degU = 3, patch_degV = 3;
    Knots uk = { 0, 0, 0, 0, 1, 1, 1, 1 };
    Knots vk = { 0, 0, 0, 0, 1, 1, 1, 1 };

    return std::make_pair(
        BezierTri(tri_pts, tri_degree),
        BSplinePatch(patch_pts, uk, vk, patch_degU, patch_degV)
    );
}

#endif // !IOUTILS_H