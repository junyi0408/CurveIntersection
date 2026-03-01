#pragma once
#ifndef COMMON_H
#define COMMON_H
#include <Eigen\Dense>
#include <openGJK\openGJK.h>
#include <iostream>

using namespace Eigen;
typedef std::vector<Vector3d> Points;
typedef std::vector<double> Knots;
typedef std::vector<int> Multis;
namespace BSpline_Utils
{
	inline void subBSplineCurve(const Points& cps, const Knots& knots, size_t p, size_t k, double t, Points& left, Points& right)
	{
		std::vector<Eigen::Vector3d> tempCPs;
		tempCPs.reserve(p + 1);
		for (int i = 0; i <= p; ++i) {
			tempCPs.push_back(cps[k - p + i]);
		}

		std::vector<Eigen::Vector3d> leftChain; leftChain.reserve(p + 1);
		std::vector<Eigen::Vector3d> rightChain; rightChain.reserve(p + 1);

		leftChain.push_back(tempCPs.front());
		rightChain.push_back(tempCPs.back());

		size_t s = (abs(t - knots[k]) < 1.0e-9) ? (size_t)1 : 0;
		size_t h = p - s;
		for (int r = 1; r <= h; ++r) {
			for (size_t i = 0; i <= p - r - s; ++i) {
				int original_idx = k - p + i + r;
				double u_start = knots[original_idx];
				double u_end = knots[k + 1 + i];
				double alpha = (t - u_start) / (u_end - u_start);
				tempCPs[i] = (1.0 - alpha) * tempCPs[i] + alpha * tempCPs[i + 1];
			}
			leftChain.push_back(tempCPs[0]);           // 新一层的最左点
			rightChain.push_back(tempCPs[p - r]);      // 新一层的最右点
		}

		left.assign(cps.begin(), cps.begin() + k - p);
		left.insert(left.end(), leftChain.begin(), leftChain.end());
		if (h < p) left.insert(left.end(), tempCPs.begin() + 1, tempCPs.begin() + (p - h + 1));

		std::reverse(rightChain.begin(), rightChain.end());
		right = rightChain;
		if (h < p) right.insert(right.begin(), tempCPs.begin(), tempCPs.begin() + (p - h));
		right.insert(right.end(), cps.begin() + k + 1, cps.end());
	}

	inline void subBSplineKnots(const Knots& knots, size_t p, size_t k, double t, Knots& left, Knots& right)
	{
		left.reserve(k + p + 2);
		right.reserve(knots.size() - k + 2 * p + 1);

		left.assign(knots.begin(), knots.begin() + k + 1);
		int existing_mult = 0;
		while (existing_mult < left.size() &&
			std::abs(left[left.size() - 1 - existing_mult] - t) < 1e-9) {
			existing_mult++;
		}
		int needed = p + 1 - existing_mult;
		for (int i = 0; i < needed; ++i) left.push_back(t);

		for (int i = 0; i <= p; ++i) right.push_back(t);
		right.insert(right.end(), knots.begin() + k + 1, knots.end());
	}
}

namespace Intersection
{
	class AABB
	{
	public:
		AABB() = default;
		~AABB() = default;
		AABB(const Points& points) {
			double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX, xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;
			for (const Vector3d& p : points) {
				if (p.x() < xmin) xmin = p.x();
				if (p.x() > xmax) xmax = p.x();
				if (p.y() < ymin) ymin = p.y();
				if (p.y() > ymax) ymax = p.y();
				if (p.z() < zmin) zmin = p.z();
				if (p.z() > zmax) zmax = p.z();
			}
			lb = Vector3d(xmin, ymin, zmin);
			ru = Vector3d(xmax, ymax, zmax);
			scale = sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) + (zmax - zmin) * (zmax - zmin));
		}
		bool IsOut(const AABB& other) const {
			if (ru(0) < other.lb(0) || lb(0) > other.ru(0)) return true;
			if (ru(1) < other.lb(1) || lb(1) > other.ru(1)) return true;
			if (ru(2) < other.lb(2) || lb(2) > other.ru(2)) return true;
			return false;
		};
	public:
		Vector3d lb, ru;
		double scale = 0;
	};


	struct GJKBodyWrapper {
		std::vector<gkFloat*> ptrs; // 存储二级指针
		// 注意：我们不需要存储数据的拷贝，只需存储指向外部数据的指针
		// 但为了安全，如果外部数据是临时的，这里最好有一份拷贝。
		// 在本例中，我们专门构建一个 Minkowski 点云，所以需要自己持有数据。
		std::vector<Vector3d> data_storage;
		gkPolytope polytope;

		// 构造函数：接受一组点，构建 openGJK 结构
		GJKBodyWrapper(const std::vector<Vector3d>& points) {
			data_storage = points; // 拷贝数据，保证内存安全
			size_t n = data_storage.size();
			ptrs.reserve(n);

			for (size_t i = 0; i < n; ++i) {
				// 强转 double* 为 gkFloat*
				ptrs.push_back(reinterpret_cast<gkFloat*>(data_storage[i].data()));
			}

			polytope.numpoints = static_cast<int>(n);
			polytope.coord = ptrs.data(); // 获取 vector 内部数组的指针
		}
	};
	inline bool checkInt_no_contact(const Points& ch1, const Points& ch2) {
		/*for (const auto& v : ch1) {
			std::cout << v.transpose() << std::endl;
		}*/
		GJKBodyWrapper body1(ch1), body2(ch2);
		gkSimplex s;
		s.nvrtx = 0;
		gkFloat dist = compute_minimum_distance(body1.polytope, body2.polytope, &s);
		return dist < 1e-9;
	}

	inline bool checkInt_common_vertex(const Points& ch1, const Points& ch2, Vector3d p) {
		/*Points cloud;
		cloud.reserve(ch1.size() + ch2.size());
		for (const auto& v : ch1) {
			Vector3d temp = v - p;
			if (temp.norm() > 1e-12)
				cloud.push_back(temp.normalized());
		}
		for (const auto& v : ch2) {
			Vector3d temp = v - p;
			if (temp.norm() > 1e-12)
				cloud.push_back(-temp.normalized());
		}*/
		double shrink_factor = 0.9999;
		Points cloud1(ch1), cloud2(ch2);
		Vector3d centroid(0, 0, 0);
		for (const auto& p : cloud1) {
			centroid += p;
		}
		centroid /= static_cast<double>(cloud1.size());
		for (auto& p : cloud1) {
			p = centroid + (p - centroid) * shrink_factor;
		}

		centroid.setZero();
		for (const auto& p : cloud2) {
			centroid += p;
		}
		centroid /= static_cast<double>(cloud2.size());
		for (auto& p : cloud2) {
			p = centroid + (p - centroid) * shrink_factor;
		}
		return checkInt_no_contact(cloud1, cloud2);
	}
}

namespace Cone_Utils
{
	struct GJKBodyWrapper {
		std::vector<gkFloat*> ptrs; // 存储二级指针
		// 注意：我们不需要存储数据的拷贝，只需存储指向外部数据的指针
		// 但为了安全，如果外部数据是临时的，这里最好有一份拷贝。
		// 在本例中，我们专门构建一个 Minkowski 点云，所以需要自己持有数据。
		std::vector<Vector3d> data_storage;
		gkPolytope polytope;

		// 构造函数：接受一组点，构建 openGJK 结构
		GJKBodyWrapper(const std::vector<Vector3d>& points) {
			data_storage = points; // 拷贝数据，保证内存安全
			size_t n = data_storage.size();
			ptrs.reserve(n);

			for (size_t i = 0; i < n; ++i) {
				// 强转 double* 为 gkFloat*
				ptrs.push_back(reinterpret_cast<gkFloat*>(data_storage[i].data()));
			}

			polytope.numpoints = static_cast<int>(n);
			polytope.coord = ptrs.data(); // 获取 vector 内部数组的指针
		}
	};
	inline bool ConeIntersect(const Points& cone1, const Points& cone2) {
		//BoundingCone bc1 = BoundingCone::Build(cone1);
		//BoundingCone bc2 = BoundingCone::Build(cone2);
		//if (SimpleConeSeparated(bc1, bc2)) {
		//	return false; // 快速剔除成功！
		//}

		Points cloud;
		cloud.reserve(cone1.size() + cone2.size());
		for (const auto& v : cone1) {
			if (v.norm() > 1e-12)
				cloud.push_back(v.normalized());
		}
		for (const auto& v : cone2) {
			if (v.norm() > 1e-12)
				cloud.push_back(-v.normalized());
		}
		if (cloud.empty()) std::cerr << "ConeIntersect!!!" << std::endl;
		GJKBodyWrapper body_diff(cloud);
		std::vector<Vector3d> origin_point = { Vector3d(0, 0, 0) };
		GJKBodyWrapper body_origin(origin_point);
		gkSimplex s;
		s.nvrtx = 0;
		gkFloat dist = compute_minimum_distance(body_diff.polytope, body_origin.polytope, &s);
		return dist < 1e-9;
	}
}

#endif // !COMMON_H
