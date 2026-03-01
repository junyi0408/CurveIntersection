#pragma once
#ifndef TRIPATCHSOLVER_H
#define TRIPATCHSOLVER_H
#include "BezierTri.h"
#include "BSplinePatch.h"

class TriPatchSolver
{
public:
	static bool checkInt(const BezierTri& tri, const BSplinePatch& patch);
	static int getTriSplitDir(const BezierTri& tri);
};

#endif // !TRIPATCHSOLVER_H