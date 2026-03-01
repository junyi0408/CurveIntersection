#pragma once
#ifndef TRITRISOLVER_H
#define TRITRISOLVER_H
#include "BezierTri.h"
#include "Common.h"
class TriTriSolver
{
public:
	static bool checkInt(const BezierTri& tri1, const BezierTri& tri2, int depth = 0);
	static int getTriSplitDir(const BezierTri& tri, bool isleft);

	static bool checkIntNoContact(const BezierTri& tri1, const BezierTri& tri2, int depth = 0);
	static bool checkIntCommonVertex(const BezierTri& tri1, const BezierTri& tri2, const Vector3d& ver);
};
#endif // !TRITRISOLVER_H



