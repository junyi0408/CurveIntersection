#include "BezierTri.h"
#include "TriTriSolver.h"
Vector3d BezierTri::map_param_coord_to_bary_coord(const Vector2d& uv) const
{
    Vector3d uvw;
    uvw[0] = uv[0];
    uvw[1] = uv[1];
    uvw[2] = 1. - uvw[0] - uvw[1];
    return uvw;
}

Vector3d BezierTri::deCasteljau(const Vector3d& uvw) const
{
    Points intermediate_points = ctrl_points;

    for (int n = (int)degree; n > 0; n--)
    {
        for (int i = 0; i < n; i++)    // i-th row
        {
            int s = map_mat_idx_to_vec_idx(i, 0);
            int t = map_mat_idx_to_vec_idx(i + 1, 0);

            for (int j = 0; j <= i; j++, s++, t++)       // j-th column
            {
                intermediate_points[s] =
                    intermediate_points[t + 1] * uvw[0] +
                    intermediate_points[s] * uvw[1] +
                    intermediate_points[t] * uvw[2];
            }
        }
    }
    return intermediate_points[0];
}

const double BezierTri::fac_[] =
//0    1    2    3    4     5      6      7       8       9        10
{ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320., 362880., 3628800. };

/// @brief Calculate the Berstein polynomial, with given degree, polar index and position.
/// @note the relation between polar index and parameter is (u^j * v^i * w^k).
double BezierTri::Berstein_polynomial(const Vector3i& ijk, const Vector3d& uvw) const
{
    if (ijk[0] < 0 || ijk[1] < 0 || ijk[2] < 0)
        return 0.0;

    return facT(degree) * pow(uvw[0], ijk[1]) * pow(uvw[1], ijk[0]) * pow(uvw[2], ijk[2])
        / (facT(ijk[0]) * facT(ijk[1]) * facT(ijk[2]));
}

Vector3d BezierTri::eval_by_Berstein(const Vector3d& uvw) const
{
    Vector3d result(0, 0, 0);
    for (Index i = 0; i <= (Index)degree; i++)
    {
        for (Index j = 0; j <= i; j++)
        {
            const Vector3d& cp = ctrl_points[map_mat_idx_to_vec_idx(i, j)];
            Vector3i ijk = map_mat_idx_to_polar(i, j, degree);
            result += cp * Berstein_polynomial(ijk, uvw);
        }
    }
    return result;
}

void BezierTri::subdivision_1_3(const Vector3d& uvw,
    std::array<Points, 3>& sub_ctrl_points, Vector3d& sub_center_ctrl_pnt) const
{
    Points intermediate_points = ctrl_points;
    for (size_t i = 0; i < 3; i++)
        sub_ctrl_points[i].resize(ctrl_points.size());

    for (uint8_t n = degree; n > 0; n--)
    {
        // put intermediate control points into subdivided triangles.
        Index diff = degree - n;
        for (Index i = 0; i <= n; i++)
        {
            Index src = map_mat_idx_to_vec_idx(i, 0);
            Index dst = map_mat_idx_to_vec_idx(i + diff, diff);
            sub_ctrl_points[0][dst] = intermediate_points[src];
        }
        for (Index j = 0; j <= n; j++)
        {
            Index src_dst = map_mat_idx_to_vec_idx(n, j);
            sub_ctrl_points[1][src_dst] = intermediate_points[src_dst];
        }
        for (Index i = 0; i <= n; i++)
        {
            Index src = map_mat_idx_to_vec_idx(i, i);
            Index dst = map_mat_idx_to_vec_idx(i + diff, i);
            sub_ctrl_points[2][dst] = intermediate_points[src];
        }
        // calculate intermediate control points
        for (Index i = 0; i < n; i++)    // i-th row
        {
            Index s = map_mat_idx_to_vec_idx(i, 0);
            Index t = map_mat_idx_to_vec_idx(i + 1, 0);

            for (Index j = 0; j <= i; j++, s++, t++)       // j-th column
            {
                intermediate_points[s] =
                    intermediate_points[t + 1] * uvw[0] +
                    intermediate_points[s] * uvw[1] +
                    intermediate_points[t] * uvw[2];
            }
        }
    }
    // put the last one control points into three triangles.
    sub_ctrl_points[0][map_mat_idx_to_vec_idx(degree, degree)] = intermediate_points[0];
    sub_ctrl_points[1][0] = intermediate_points[0];
    sub_ctrl_points[2][map_mat_idx_to_vec_idx(degree, 0)] = intermediate_points[0];
    sub_center_ctrl_pnt = intermediate_points[0];
}

/// @brief Apply subdivision_1_2 operation on a Bezier triangle (one Bezier triangle is subdivided to two Bezier triangles).
/// @param [in] local_edge_idx The edge of Bezier triangle to split.
/// @param [in] uv The barycentric coordinate of subdivided center on the edge.
/// @param [out] sub_ctrl_points Three groups of control points of three sub-triangles.
/// @param [out] sub_center_ctrl_pnt Subdivide point.
void BezierTri::subdivision_1_2(const Index local_edge_idx, const Vector2d& uv,
    std::array<Points, 2>& sub_ctrl_points, Vector3d& sub_split_ctrl_pnt) const
{
    std::array<Points, 3> _sub_ctrl_points;
    Vector3d uvw;
    switch (local_edge_idx)
    {
    case 0:
    {
        uvw = { 0, uv[0], uv[1] };
        subdivision_1_3(uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
        sub_ctrl_points[0] = _sub_ctrl_points[2];
        sub_ctrl_points[1] = _sub_ctrl_points[1];
    }break;
    case 1:
    {
        uvw = { uv[1], 0, uv[0] };
        subdivision_1_3(uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
        sub_ctrl_points[0] = _sub_ctrl_points[0];
        sub_ctrl_points[1] = _sub_ctrl_points[2];
    }break;
    case 2:
    {
        uvw = { uv[0], uv[1], 0 };
        subdivision_1_3(uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
        sub_ctrl_points[0] = _sub_ctrl_points[1];
        sub_ctrl_points[1] = _sub_ctrl_points[0];
    }break;
    }
}

Vector3d BezierTri::sub_control_point(int i, int j, uint32_t k0, uint32_t k1, uint32_t k2) const
{
    if (k0 > 0)
    {
        return 0.5 * (sub_control_point(i, j, k0 - 1, k1, k2) +
            sub_control_point(i + 1, j, k0 - 1, k1, k2));
    }
    if (k1 > 0)
    {
        return 0.5 * (sub_control_point(i + 1, j, 0, k1 - 1, k2) +
            sub_control_point(i + 1, j + 1, 0, k1 - 1, k2));
    }
    if (k2 > 0)
    {
        return 0.5 * (sub_control_point(i + 1, j + 1, 0, 0, k2 - 1) +
            sub_control_point(i, j, 0, 0, k2 - 1));
    }
    return ctrl_points[map_mat_idx_to_vec_idx(i, j)];
}

void BezierTri::subdivision_1_4(std::array<Points, 4>& sub_ctrl_points) const
{
    for (size_t i = 0; i < 4; i++)
    {
        sub_ctrl_points[i].clear();
        sub_ctrl_points[i].reserve(ctrl_points.size());
    }

    for (uint32_t i = 0; i <= degree; i++)
    {
        Vector2i ij = map_polar_to_mat_idx(Vector3i(degree - i, 0, 0), degree - i);
        for (uint32_t n = 0; n <= i; n++)
            sub_ctrl_points[0].push_back(sub_control_point(ij[0], ij[1], i - n, 0, n));

        ij = map_polar_to_mat_idx(Vector3i(0, 0, degree - i), degree - i);
        for (uint32_t n = 0; n <= i; n++)
            sub_ctrl_points[1].push_back(sub_control_point(ij[0], ij[1], n, i - n, 0));

        ij = map_polar_to_mat_idx(Vector3i(0, degree - i, 0), degree - i);
        for (uint32_t n = 0; n <= i; n++)
            sub_ctrl_points[2].push_back(sub_control_point(ij[0], ij[1], 0, n, i - n));
    }

    for (uint32_t i = 0; i <= degree; i++)
    {
        for (uint32_t j = 0; j <= i; j++)
        {
            // k0,k1,k2 are slightly different from polar index.
            sub_ctrl_points[3].push_back(sub_control_point(0, 0, degree - i, i - j, j));
        }
    }
}

Points BezierTri::gradient_vectors(const int local_edge_idx) const
{
    Points vectors; vectors.reserve((degree * (degree + 1)) / 2);
    switch (local_edge_idx)
    {
    case 0:
    {
        for (int i = 1; i <= (int)degree; i++)
        {
            int s = map_mat_idx_to_vec_idx(i, 0);
            for (int j = 0; j < i; j++, s++)
                vectors.push_back(ctrl_points[s + 1] - ctrl_points[s]);
        }
    }break;
    case 1:
    {
        for (int i = 1; i <= (int)degree; i++)
        {
            int s = map_mat_idx_to_vec_idx(i, 0);
            int t = map_mat_idx_to_vec_idx(i - 1, 0);
            for (int j = 0; j < i; j++, s++, t++)
                vectors.push_back(ctrl_points[t] - ctrl_points[s]);
        }
    }break;
    case 2:
    {
        for (int i = 1; i <= (int)degree; i++)
        {
            int s = map_mat_idx_to_vec_idx(i, 0);
            for (int j = 0; j < i; j++, s++)
                vectors.push_back(ctrl_points[s] - ctrl_points[s + 1]);
        }
    }break;
    default:
        break;
    }
    return vectors;
}


int BezierTri::getBestSplitDir() const
{
    double s1 = 0, s2 = 0, s3 = 0;
    for (int i = 1; i <= degree; i++)
    {
		Vector3d temp = ctrl_points[map_mat_idx_to_vec_idx(i, 0)] - ctrl_points[map_mat_idx_to_vec_idx(i - 1, 0)];
		s1 += temp.norm();
    }
    for (int i = 1; i <= degree; i++)
    {
        Vector3d temp = ctrl_points[map_mat_idx_to_vec_idx(degree, i)] - ctrl_points[map_mat_idx_to_vec_idx(degree, i - 1)];
        s2 += temp.norm();
    }
    for (int i = 1; i <= degree; i++)
    {
        Vector3d temp = ctrl_points[map_mat_idx_to_vec_idx(i, i)] - ctrl_points[map_mat_idx_to_vec_idx(i - 1, i - 1)];
        s3 += temp.norm();
    }
	if (s1 > s2 && s1 > s3) return 0;
	if (s2 > s1 && s2 > s3) return 1;
	return 2;
}

bool BezierTri::selfInt(int depth) const
{
    Points gradU, gradV, gradNegV;
    gradU = gradient_vectors(0);
    gradV = gradient_vectors(1);
    gradNegV.resize(gradV.size());
    std::transform(gradV.begin(), gradV.end(), gradNegV.begin(), [](const Vector3d& v) {return -v; });

    bool cone = !Cone_Utils::ConeIntersect(gradU, gradV) && !Cone_Utils::ConeIntersect(gradU, gradNegV);

    if (cone)
        return false;
    /*AABB box(ctrl_points);
    if (box.scale < 1.0e-6) {
        if (stats) stats->aabb_pruned++;
        return false;
    }*/
    if (depth >= 5)
        return true;

    int dir = getBestSplitDir();
	std::array < Points, 2> sub_ctrl_points;
	Vector3d sub_split_ctrl_pnt;
	subdivision_1_2(dir, Vector2d(0.5, 0.5), sub_ctrl_points, sub_split_ctrl_pnt);
    BezierTri b1(sub_ctrl_points[0], degree), b2(sub_ctrl_points[1], degree);

    int index = 0;
    if (dir == 0) index = 2;
    if (dir == 2) index = 1;
    bool intersect = TriTriSolver::checkInt(b1, b2);
    if (intersect)
        return true;
    return b1.selfInt(depth + 1) || b2.selfInt(depth + 1);
}