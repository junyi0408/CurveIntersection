#include "BSplinePatch.h"
#include "PatchPatchSolver.h"
void BSplinePatch::gradient_vectors(size_t dir, Points& vecs) const
{
	vecs.clear();
	if (dir) {
		vecs.reserve(NUPoles() * (NVPoles() - 1));
		for (size_t i = 0; i < NUPoles(); i++)
		{
			for (size_t j = 1; j < NVPoles(); j++)
			{
				Vector3d temp = CntPoint(i, j) - CntPoint(i, j - 1);
				vecs.push_back(temp);
			}
		}
	}
	else
	{
		vecs.reserve((NUPoles() - 1) * NVPoles());
		for (size_t i = 1; i < NUPoles(); i++)
		{
			for (size_t j = 0; j < NVPoles(); j++)
			{
				Vector3d temp = CntPoint(i, j) - CntPoint(i - 1, j);
				vecs.push_back(temp);
			}
		}
	}
}

void BSplinePatch::sub1to2(size_t dir, double ratio, BSplinePatch& p1, BSplinePatch& p2, bool evalKnots) const
{
	Points Lu, Ru;
	Knots lu, ru;
	double t = ratio;
	if (dir % 2 == 1)   // v
	{
		if (evalKnots) t = (1.0 - ratio) * vknots[0] + ratio * vknots.back();
		auto it = std::upper_bound(vknots.begin(), vknots.end(), t);
		size_t k = (size_t)std::distance(vknots.begin(), it) - 1;
		if (abs(t - vknots[k]) < 1.0e-9) {
			t = vknots[k];
		}
		size_t K = abs(t - vknots[k]) < 1.0e-9 ? k : k + 1;
		Lu.resize(NUPoles() * K);
		Ru.resize(NUPoles() * (NVPoles() - k + vdegree));
		Points col(NVPoles()), Lc, Rc;
		for (uint32_t i = 0; i < NUPoles(); ++i) {
			for (uint32_t j = 0; j < NVPoles(); ++j) {
				col[j] = CntPoint(i, j);
			}
			BSpline_Utils::subBSplineCurve(col, vknots, vdegree, k, t, Lc, Rc);
			for (uint32_t j = 0; j < K; ++j) {
				Lu[i * K + j] = Lc[j];
			}
			for (uint32_t j = 0; j < (NVPoles() - k + vdegree); ++j) {
				Ru[i * (NVPoles() - k + vdegree) + j] = Rc[j];
			}
		}
		BSpline_Utils::subBSplineKnots(vknots, vdegree, k, t, lu, ru);
		p1 = BSplinePatch(Lu, uknots, lu, udegree, vdegree);
		p2 = BSplinePatch(Ru, uknots, ru, udegree, vdegree);
	}
	else   // u
	{
		if (evalKnots) t = (1.0 - ratio) * uknots[0] + ratio * uknots.back();
		auto it = std::upper_bound(uknots.begin(), uknots.end(), t);
		size_t k = (size_t)std::distance(uknots.begin(), it) - 1;
		if (abs(t - uknots[k]) < 1.0e-9) {
			t = uknots[k];
		}
		size_t K = abs(t - uknots[k]) < 1.0e-9 ? k : k + 1;
		Lu.resize(K * NVPoles());
		Ru.resize((NUPoles() - k + udegree) * NVPoles());
		Points col(NUPoles()), Lc, Rc;
		for (size_t j = 0; j < NVPoles(); j++) {
			for (size_t i = 0; i < NUPoles(); i++) {
				col[i] = CntPoint(i, j);
			}
			BSpline_Utils::subBSplineCurve(col, uknots, udegree, k, t, Lc, Rc);
			for (size_t i = 0; i < K; ++i) {
				Lu[i * NVPoles() + j] = Lc[i];
			}
			for (size_t i = 0; i < (NUPoles() - k + udegree); ++i) {
				Ru[i * NVPoles() + j] = Rc[i];
			}
		}
		BSpline_Utils::subBSplineKnots(uknots, udegree, k, t, lu, ru);
		p1 = BSplinePatch(Lu, lu, vknots, udegree, vdegree);
		p2 = BSplinePatch(Ru, ru, vknots, udegree, vdegree);
	}
}

void BSplinePatch::sub1to4(std::array<BSplinePatch, 4>& sub_patches) const
{
	BSplinePatch p_u_left, p_u_right;

	this->sub1to2(0, 0.5, p_u_left, p_u_right, true);
	p_u_left.sub1to2(1, 0.5, sub_patches[0], sub_patches[1], true);
	p_u_right.sub1to2(1, 0.5, sub_patches[2], sub_patches[3], true);
}


int BSplinePatch::getBestSplitDir() const
{
	double vlen = -1;
	for (size_t i = 0; i < NUPoles(); i++) {
		double sum = 0.0;
		for (size_t j = 1; j < NVPoles(); j++) {
			sum += (CntPoint(i, j) - CntPoint(i, j - 1)).norm();
		}
		vlen = std::max(vlen, sum);
	}
	double ulen = -1;
	for (size_t j = 0; j < NVPoles(); j++) {
		double sum = 0.0;
		for (size_t i = 1; i < NUPoles(); i++) {
			sum += (CntPoint(i, j) - CntPoint(i - 1, j)).norm();
		}
		ulen = std::max(ulen, sum);
	}
	return ulen > vlen ? 0 : 1;
}


bool BSplinePatch::selfInt(int depth) const
{
	Points gradU, gradV, gradNegV;
	gradient_vectors(0, gradU);
	gradient_vectors(1, gradV);
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
	BSplinePatch b1, b2, b3, b4;

	int dir = getBestSplitDir();
	sub1to2(dir, 0.45, b1, b2);
	sub1to2(dir, 0.55, b3, b4);

	bool intersect = PatchPatchSolver::checkIntNoContact(b1, b4);
	if (intersect)
		return true;
	return b2.selfInt(depth + 1) || b3.selfInt(depth + 1);
}