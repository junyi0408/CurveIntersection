#include "TriPatchSolver.h"

int TriPatchSolver::getTriSplitDir(const BezierTri& tri)
{
	double sum0 = 0.0, sum1 = 0.0;
	for (int i = 1; i <= tri.degree; i++)
	{
		Vector3d temp = tri.ctrl_points[tri.map_mat_idx_to_vec_idx(i, 0)] - tri.ctrl_points[tri.map_mat_idx_to_vec_idx(i - 1, 0)];
		sum0 += temp.norm();
	}
	for (int j = 1; j <= tri.degree; j++)
	{
		Vector3d temp = tri.ctrl_points[tri.map_mat_idx_to_vec_idx(tri.degree, j)] - tri.ctrl_points[tri.map_mat_idx_to_vec_idx(tri.degree, j - 1)];
		sum1 += temp.norm();
	}
	return sum0 > sum1 ? 0 : 1;
}

bool TriPatchSolver::checkInt(const BezierTri& tri, const BSplinePatch& patch)
{
	Points left, edge, right;
	left = tri.gradient_vectors(2);
	patch.gradient_vectors(1, right);
	edge.reserve(patch.NUPoles() - 1);
	for (size_t i = 1; i < patch.NUPoles(); i++)
		edge.push_back(patch.CntPoint(i, 0) - patch.CntPoint(i - 1, 0));
	Points L, R;
	L.reserve(left.size() + edge.size());
	R.reserve(right.size() + edge.size());
	L.insert(L.end(), left.begin(), left.end());
	L.insert(L.end(), edge.begin(), edge.end());
	R.insert(R.end(), right.begin(), right.end());
	R.insert(R.end(), edge.begin(), edge.end());

	bool conetest = !Cone_Utils::ConeIntersect(L, right) && !Cone_Utils::ConeIntersect(R, left);

	if (conetest)
		return false;
	/*AABB box(p1.ctrl_points);
	if (box.scale < 1.0e-6)
		return false;*/
	std::array<Points, 2> sub_ctrl_points;
	Vector3d temp;
	int local_edge_index = getTriSplitDir(tri);
	tri.subdivision_1_2(local_edge_index, Vector2d(0.5, 0.5), sub_ctrl_points, temp);
	BezierTri t1, t2;
	if (local_edge_index == 0)
	{
		t1 = BezierTri(sub_ctrl_points[0], tri.degree);
		t2 = BezierTri(sub_ctrl_points[1], tri.degree);
	}
	else
	{
		t1 = BezierTri(sub_ctrl_points[1], tri.degree);
		t2 = BezierTri(sub_ctrl_points[0], tri.degree);
	}
	BSplinePatch p1, p2;
	patch.sub1to2(1, 0.1, p1, p2);

	if (Intersection::checkInt_no_contact(tri.ctrl_points, p2.ctrl_points))
		return true;
	if (local_edge_index == 0)
	{
		if (Intersection::checkInt_common_vertex(t2.ctrl_points, patch.ctrl_points, patch.CntPoint(patch.udegree, 0)))
			return true;
	}
	else
	{
		if (Intersection::checkInt_common_vertex(t2.ctrl_points, patch.ctrl_points, patch.CntPoint(0, 0)))
			return true;
	}
	return checkInt(t1, p1);
}