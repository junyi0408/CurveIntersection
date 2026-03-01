#pragma once
#ifndef BSPLINEPATCH_H
#define BSPLINEPATCH_H
#include "Common.h"

class BSplinePatch
{
public:
	BSplinePatch() = default;
	~BSplinePatch() = default;
	BSplinePatch(const Points& cp, const Knots& uk, const Knots& vk, int ud, int vd) :
		ctrl_points(cp), uknots(uk), vknots(vk), udegree(ud), vdegree(vd) {};

	inline size_t NUPoles() const { return uknots.size() - udegree - 1; };
	inline size_t NVPoles() const { return vknots.size() - vdegree - 1; };
	inline size_t map_mat_to_idx(size_t i, size_t j) const { return i * NVPoles() + j; };
	inline size_t uspan() const { return uknots.size() - 2 * udegree - 1; };
	inline size_t vspan() const { return vknots.size() - 2 * vdegree - 1; };
	const Vector3d& CntPoint(size_t i, size_t j) const { return ctrl_points[map_mat_to_idx(i, j)]; };
	/*const LightBezier& GetCachedBezier(const Region& r) const {
		size_t vspan = vknots.size() - 2 * vdegree - 1;
		size_t id = (r.ustart - udegree) * vspan + r.vstart - vdegree;
		return bezier_cache[id];
	}
	void getCntByRegion(const Region& r, Points& points) const;
	void getBezierByRegion(const Region& r, LightBezier& bezier) const;
	void BuildBezierCache(const std::vector<uint8_t>& mask) const;
	void BuildBezierCacheGlobal() const;
	int getBestSplitDir() const;

	bool selfIntersect(size_t depth, IntersectionStats* stats) const;
	bool selfIntersect() const { return selfIntersect(0, nullptr); };*/

	// dir = 0 for u,dir = 1 for v
	void gradient_vectors(size_t dir, Points& vecs) const;
	void sub1to2(size_t dir, double ratio, BSplinePatch& p1, BSplinePatch& p2, bool evalKnots = true) const;
	void sub1to4(std::array<BSplinePatch, 4>& sub_patches) const;

	bool selfInt(int depth = 0) const;
	int getBestSplitDir() const;

public:
	Points ctrl_points;
	Knots uknots, vknots;
	size_t udegree, vdegree;

	/*mutable std::vector<LightBezier> bezier_cache;
	mutable bool is_cache_built = false;*/
};
#endif // !BSPLINEPATCH_H



