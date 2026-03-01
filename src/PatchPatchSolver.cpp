#include "PatchPatchSolver.h"
bool PatchPatchSolver::checkIntNoContact(const BSplinePatch& p1, const BSplinePatch& p2, int depth)
{
	Intersection::AABB ab1(p1.ctrl_points), ab2(p2.ctrl_points);
	if (ab1.IsOut(ab2))
		return false;
	if (!Intersection::checkInt_no_contact(p1.ctrl_points, p2.ctrl_points))
		return false;
	if (depth >= 15)
		return true;

	double sizeA = ab1.scale;
	double sizeB = ab2.scale;
	if (sizeA > sizeB) {
		// 细分 A 为 3 个子三角形
		std::array<BSplinePatch, 4> subA;
		p1.sub1to4(subA);

		// 只要 B 和 A 的任意一个子块相交，就说明整体相交
		for (int i = 0; i < 4; ++i) {
			if (checkIntNoContact(subA[i], p2, depth + 1)) {
				return true;
			}
		}
	}
	else {
		// 细分 B 为 3 个子三角形
		std::array<BSplinePatch, 4> subB;
		p2.sub1to4(subB);

		for (int i = 0; i < 4; ++i) {
			if (checkIntNoContact(p1, subB[i], depth + 1)) {
				return true;
			}
		}
	}

	// 所有子块都不交
	return false;
}

bool PatchPatchSolver::checkInt(const BSplinePatch& p1, const BSplinePatch& p2, int depth)
{
	Points left, edge, right;
	p1.gradient_vectors(1, left);
	p2.gradient_vectors(1, right);
	std::transform(left.begin(), left.end(), left.begin(), std::negate<>());
	edge.reserve(p1.NUPoles() - 1);
	for (size_t i = 1; i < p1.NUPoles(); i++)
		edge.push_back(p2.CntPoint(i, 0) - p2.CntPoint(i - 1, 0));
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
	if (depth >= 5)
		return true;

	BSplinePatch b1, b2, b3, b4;
	p1.sub1to2(1, 0.9, b1, b2);
	p2.sub1to2(1, 0.1, b3, b4);

	bool flag = checkIntNoContact(b1, p2) || checkIntNoContact(b4, p1);
	if (flag)
		return true;
	return checkInt(b2, b3, depth + 1);
}