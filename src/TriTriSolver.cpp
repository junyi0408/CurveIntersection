#include "TriTriSolver.h"

bool TriTriSolver::checkIntCommonVertex(const BezierTri& tri1, const BezierTri& tri2, const Vector3d& ver)
{
	Points points(tri1.ctrl_points);
	Vector3d temp(0, 0, 0);
	for(const auto& p : points)
		temp += p - ver;
	for (auto& p : points) {
		if((p - ver).norm() < 1e-9)
			p += 0.001 * temp.normalized();
	}
	BezierTri tempTri(points, tri1.degree);
	return checkIntNoContact(tempTri, tri2);
}

bool TriTriSolver::checkIntNoContact(const BezierTri& tri1, const BezierTri& tri2, int depth)
{
	Intersection::AABB ab1(tri1.ctrl_points), ab2(tri2.ctrl_points);
	if (ab1.IsOut(ab2))
		return false;
	if (!Intersection::checkInt_no_contact(tri1.ctrl_points, tri2.ctrl_points))
		return false;
	if (depth >= 5)
		return true;

	double sizeA = ab1.scale;
	double sizeB = ab2.scale;
	if (sizeA > sizeB) {
		// 细分 A 为 3 个子三角形
		std::array<Points, 4> subA;
		tri1.subdivision_1_4(subA);

		// 只要 B 和 A 的任意一个子块相交，就说明整体相交
		for (int i = 0; i < 4; ++i) {
			BezierTri subTriA(subA[i], tri1.degree);
			if (checkIntNoContact(subTriA, tri2, depth + 1)) {
				return true;
			}
		}
	}
	else {
		// 细分 B 为 3 个子三角形
		std::array<Points, 4> subB;
		tri2.subdivision_1_4(subB);

		for (int i = 0; i < 4; ++i) {
			BezierTri subTriB(subB[i], tri2.degree);
			if (checkIntNoContact(tri1, subTriB, depth + 1)) {
				return true;
			}
		}
	}

	// 所有子块都不交
	return false;
}

int TriTriSolver::getTriSplitDir(const BezierTri& tri, bool isleft)
{
	if (isleft) 
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
	else
	{
		double sum0 = 0.0, sum1 = 0.0;
		for (int i = 1; i <= tri.degree; i++)
		{
			Vector3d temp = tri.ctrl_points[tri.map_mat_idx_to_vec_idx(i, i)] - tri.ctrl_points[tri.map_mat_idx_to_vec_idx(i - 1, i - 1)];
			sum0 += temp.norm();
		}
		for (int j = 1; j <= tri.degree; j++)
		{
			Vector3d temp = tri.ctrl_points[tri.map_mat_idx_to_vec_idx(tri.degree, j)] - tri.ctrl_points[tri.map_mat_idx_to_vec_idx(tri.degree, j - 1)];
			sum1 += temp.norm();
		}
		return sum0 > sum1 ? 2 : 1;
	}
}

bool TriTriSolver::checkInt(const BezierTri& tri1, const BezierTri& tri2, int depth)
{
	Points left, edge, right;
	left = tri1.gradient_vectors(2);
	right = tri2.gradient_vectors(0);
	edge.reserve(tri1.degree);
	for (size_t i = 1; i <= tri2.degree; i++)
		edge.push_back(tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(i, 0)] - tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(i - 1, 0)]);
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
	/*Intersection::AABB box(tri1.ctrl_points);
	if (box.scale < 1.0e-6)
		return true;*/
	if (depth >= 5)
		return true;
	std::array<Points, 2> sub_ctrl_points;
	Vector3d temp;
	int local_edge_index1 = getTriSplitDir(tri1, true);
	int local_edge_index2 = getTriSplitDir(tri2, false);
	tri1.subdivision_1_2(local_edge_index1, Vector2d(0.5, 0.5), sub_ctrl_points, temp);
	BezierTri t1, t2, t3, t4;
	if (local_edge_index1 == 0)
	{
		t1 = BezierTri(sub_ctrl_points[1], tri1.degree);
		t2 = BezierTri(sub_ctrl_points[0], tri1.degree);
	}
	else
	{
		t1 = BezierTri(sub_ctrl_points[0], tri1.degree);
		t2 = BezierTri(sub_ctrl_points[1], tri1.degree);
	}
	tri2.subdivision_1_2(local_edge_index2, Vector2d(0.5, 0.5), sub_ctrl_points, temp);
	if (local_edge_index2 == 2)
	{
		t3 = BezierTri(sub_ctrl_points[1], tri2.degree);
		t4 = BezierTri(sub_ctrl_points[0], tri2.degree);
	}
	else
	{
		t3 = BezierTri(sub_ctrl_points[0], tri2.degree);
		t4 = BezierTri(sub_ctrl_points[1], tri2.degree);
	}
	
	if (local_edge_index1 == 0) {
		if (checkIntCommonVertex(t1,tri2, tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(tri2.degree, 0)]))
			return true;
	}
	else {
		if (checkIntCommonVertex(t1, tri2, tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(0, 0)]))
			return true;
	}
	if (local_edge_index2 == 2) {
		if (checkIntCommonVertex(t4, tri1, tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(tri2.degree, 0)]))
			return true;
	}
	else {
		if (checkIntCommonVertex(t4, tri1, tri2.ctrl_points[tri2.map_mat_idx_to_vec_idx(0, 0)]))
			return true;
	}
	return checkInt(t2, t3, depth + 1);
}