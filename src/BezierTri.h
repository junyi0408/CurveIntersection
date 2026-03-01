#pragma once
#ifndef BEZIERTRI_H
#define BEZIERTRI_H

#include <vector>
#include "Common.h"

class BezierTri
{
public:
	BezierTri() = default;
    BezierTri(const Points& p, const int d) : ctrl_points(p), degree(d) {};
    /* Index Map & Conversion */
    Vector3d map_param_coord_to_bary_coord(const Vector2d& uv) const;

    // map (i, j) in lower triangle matrix to the index in linear array.
    inline Index map_mat_idx_to_vec_idx(const Index i, const Index j) const { return (i * (i + 1)) / 2 + j; }

    // map (i, j) in lower triangle matrix to the polar form (i, j, k). n is degree.
    inline Vector3i map_mat_idx_to_polar(const Index i, const Index j, const uint32_t n) const { return Vector3i(n - i, j, i - j); }

    // map polar form to index in lower triangle matrix. n is degree.
    inline Vector2i map_polar_to_mat_idx(const Vector3i& ijk, const uint32_t n)  const { return Vector2i(n - ijk[0], ijk[1]); }

    Vector3d deCasteljau(const Vector3d& uvw) const;

    static const double fac_[];
    static double facT(int n)
    {
        // NOTE: in most cases, subscript n won't overflow. 
        return fac_[n];
    }
    double Berstein_polynomial(const Vector3i& ijk, const Vector3d& uvw) const;
    Vector3d eval_by_Berstein(const Vector3d& uvw) const;

    /* Subdivision */
    void subdivision_1_2(const Index local_edge_idx, const Vector2d& uv,
        std::array<Points, 2>& sub_ctrl_points, Vector3d& sub_split_ctrl_pnt) const;

    void subdivision_1_3(const Vector3d& uvw,
        std::array<Points, 3>& sub_ctrl_points, Vector3d& sub_center_ctrl_pnt) const;

    void subdivision_1_4(std::array<Points, 4>& sub_ctrl_points) const;

    Vector3d sub_control_point(int i, int j, uint32_t k0, uint32_t k1, uint32_t k2) const;


    /* Gradient & Jacobian */

    Points gradient_vectors(const int local_edge_idx) const;

    bool selfInt(int depth = 0) const;
    int getBestSplitDir() const;
public:
    Points ctrl_points;
    uint32_t degree;
};

#endif // !BEZIERTRI_H



