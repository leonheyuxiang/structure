#ifndef TIN_TRIANGLE_H
#define TIN_TRIANGLE_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// macros for accessing triangle corners

#define TIN_TRIANGLE(i) (i/3)
#define TIN_CORNER(i) (i%3)
#define TIN_INDEX(t,c) (t*3+c)
#define TIN_NEXT(c) ((c+1)%3)
#define TIN_PREV(c) ((c+2)%3)

#define POINT_BRIO_BUFFER 10000
#define MAX_LOCATE_STEPS 10000

#define stkPOP() (dfs_stack[--dfs_stack_cur])
#define stkEMPTY() (dfs_stack_cur == 0)
#define stkPUSH(value) { \
  if (dfs_stack_cur == dfs_stack_max) {dfs_stack_max += 500; dfs_stack = (int*)realloc(dfs_stack, sizeof(int)*dfs_stack_max);} \
  dfs_stack[dfs_stack_cur++] = value; \
}
#define stkDESTROY() { if (dfs_stack) free(dfs_stack); dfs_stack = 0; dfs_stack_cur = 0; dfs_stack_max = 0; }

class  __declspec(dllexport) TINtriangle
{
public:
	float* V[3];
	int N[3];
	int next;

public:
	bool IsValid() const
	{
		if (this->V[0])
			return true;
		else
			return false;
	}
};

class  __declspec(dllexport) TINMesh
{
public:
	TINMesh() {};
	TINMesh(int npoints, float* buffer);
	~TINMesh();

	void Initialize(int num_of_points, bool convex_hull = false);
	void AddPoint(float* p);
	void Derive();
	void Destory();

	TINtriangle* TINlocate_brute_query(float* p);
	TINtriangle* TINlocate_query(float* p);

	int Size();
	TINtriangle* GetTriangle(int t);

	bool IsValidTriangle(const TINtriangle* t);
	bool IsValidTriangle(int index);

	int TriangleIndex(TINtriangle* t_locate);
	int* PointIndices(TINtriangle* t_locate);

public:
	TINtriangle* triangle_buffer = 0;
	int triangle_buffer_size = 0;
	int triangle_buffer_alloc = 0;
	int triangle_next;
	int triangle_newest;
	bool initialized = false;
	bool buffer_full = false;
	bool convex_hull_only = false;
	float* points[POINT_BRIO_BUFFER];
	int pointer;
	int* dfs_stack = 0;
	int dfs_stack_cur = 0;
	int dfs_stack_max = 0;

public:
	float* point_buffer = 0;

protected:
	double incirclef(const float *pa, const float* pb, const float* pc, const float* pd);
	double orient2df(const float *pa, const float* pb, const float* pc);
	double orient2df(const float *pa, const double* dpb, const float* pc);
	bool Init(float* v0, float* v1, float* v2);
	int InSegment(const float* a, const float* b, const float* p);
	bool IsOnSegment(const float* a, const float* b, const float* p);
	bool IsInsideCircle(const TINtriangle* t, const float* p);
	TINtriangle* BrutePointLocation(float* p);
	TINtriangle* SpecialPointLocation(const float* p, TINtriangle* t);
	TINtriangle* PointLocation(float* p);
	void Update(float* p, TINtriangle* t);
	bool InsertPoint(float* p);
};

#endif