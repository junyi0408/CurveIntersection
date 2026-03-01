#pragma once
#ifndef PATCHPATCHSOLVER_H
#define PATCHPATCHSOLVER_H
#include "BSplinePatch.h"
#include "Common.h"

class PatchPatchSolver
{
public:
	static bool checkInt(const BSplinePatch& p1, const BSplinePatch& p2, int depth = 0);

	static bool checkIntNoContact(const BSplinePatch& p1, const BSplinePatch& p2, int depth = 0);
};

#endif // !PATCHPATCHSOLVER_H



