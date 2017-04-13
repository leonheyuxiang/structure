#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "predicates.h"
#include "TINTriangle.h"

extern void exactinit();
extern double incircle(const double *pa, const double* pb, const double* pc, const double* pd);
extern double orient2d(const double *pa, const double* pb, const double* pc);

TINMesh::TINMesh(int npoints, float* buffer)
{
	point_buffer = (float*)malloc(sizeof(float)* 3 * npoints);
	Initialize(npoints);

	for (int i = 0; i < npoints; ++i)
	{
		point_buffer[3 * i + 0] = buffer[3 * i + 0];
		point_buffer[3 * i + 1] = buffer[3 * i + 1];
		point_buffer[3 * i + 2] = buffer[3 * i + 2];
		// add the point to the triangulation
		AddPoint(&(point_buffer[3 * i]));
	}
	Derive();
}

TINMesh::~TINMesh()
{
	Destory();
	if (point_buffer != NULL)
		free(point_buffer);
}

double TINMesh::incirclef(const float *pa, const float* pb, const float* pc, const float* pd)
{
	double dpa[2];
	double dpb[2];
	double dpc[2];
	double dpd[2];
	dpa[0] = pa[0]; dpa[1] = pa[1];
	dpb[0] = pb[0]; dpb[1] = pb[1];
	dpc[0] = pc[0]; dpc[1] = pc[1];
	dpd[0] = pd[0]; dpd[1] = pd[1];
	return incircle(dpa, dpb, dpc, dpd);
};

double TINMesh::orient2df(const float *pa, const float* pb, const float* pc)
{
	double dpa[2];
	double dpb[2];
	double dpc[2];
	dpa[0] = pa[0]; dpa[1] = pa[1];
	dpb[0] = pb[0]; dpb[1] = pb[1];
	dpc[0] = pc[0]; dpc[1] = pc[1];
	return orient2d(dpa, dpb, dpc);
};

double TINMesh::orient2df(const float *pa, const double* dpb, const float* pc)
{
	double dpa[2];
	double dpc[2];
	dpa[0] = pa[0]; dpa[1] = pa[1];
	dpc[0] = pc[0]; dpc[1] = pc[1];
	return orient2d(dpa, dpb, dpc);
};

bool TINMesh::Init(float* v0, float* v1, float* v2)
{
	int i, j, k;
	double orient = orient2df(v0, v1, v2);
	// create the finite triangle
	triangle_buffer[0].next = -1; // in use
	if (orient > 0)
	{
		triangle_buffer[0].V[0] = v0;
		triangle_buffer[0].V[1] = v1;
		triangle_buffer[0].V[2] = v2;
	}
	else if (orient < 0)
	{
		triangle_buffer[0].V[0] = v1;
		triangle_buffer[0].V[1] = v0;
		triangle_buffer[0].V[2] = v2;
	}
	else
	{
		return false;
	}
	// create the other three infinite triangles
	triangle_buffer[1].next = -1; // in use
	triangle_buffer[1].V[0] = 0;
	triangle_buffer[1].V[1] = triangle_buffer[0].V[2];
	triangle_buffer[1].V[2] = triangle_buffer[0].V[1];

	triangle_buffer[2].next = -1; // in use
	triangle_buffer[2].V[0] = 0;
	triangle_buffer[2].V[1] = triangle_buffer[0].V[0];
	triangle_buffer[2].V[2] = triangle_buffer[0].V[2];

	triangle_buffer[3].next = -1; // in use
	triangle_buffer[3].V[0] = 0;
	triangle_buffer[3].V[1] = triangle_buffer[0].V[1];
	triangle_buffer[3].V[2] = triangle_buffer[0].V[0];
	// set up neighbor pointers
	for (i = 0; i<4; i++)
	{
		for (j = i + 1; j<4; j++)
		{
			int c1 = -1;
			int c2 = -1;
			for (k = 0; k<3; k++)
			{
				if (triangle_buffer[i].V[k] != triangle_buffer[j].V[0] && triangle_buffer[i].V[k] != triangle_buffer[j].V[1] && triangle_buffer[i].V[k] != triangle_buffer[j].V[2])
				{
					assert(c1 == -1);
					c1 = k;
				}
			}
			assert(c1 != -1);
			for (k = 0; k<3; k++)
			{
				if (triangle_buffer[j].V[k] != triangle_buffer[i].V[0] && triangle_buffer[j].V[k] != triangle_buffer[i].V[1] && triangle_buffer[j].V[k] != triangle_buffer[i].V[2])
				{
					assert(c2 == -1);
					c2 = k;
				}
			}
			assert(c2 != -1);
			triangle_buffer[i].N[c1] = TIN_INDEX(j, c2);
			triangle_buffer[j].N[c2] = TIN_INDEX(i, c1);
		}
	}
	triangle_newest = 0;
	triangle_next = 4;
	triangle_buffer_size = 4;
	return true;
}

int TINMesh::InSegment(const float* a, const float* b, const float* p)
{
	if (a[0] < b[0])
	{
		if ((p[0] < a[0]) || (p[0] > b[0]))
		{
			return 0;
		}
		else if (p[0] == a[0])
		{
			return 1;
		}
		else if (p[0] == b[0])
		{
			return 2;
		}
		else
		{
			return 3;
		}
	}
	else if (a[0] > b[0])
	{
		if ((p[0] > a[0]) || (p[0] < b[0]))
		{
			return 0;
		}
		else if (p[0] == a[0])
		{
			return 1;
		}
		else if (p[0] == b[0])
		{
			return 2;
		}
		else
		{
			return 3;
		}
	}
	else if (a[1] < b[1])
	{
		if ((p[1] < a[1]) || (p[1] > b[1]))
		{
			return 0;
		}
		else if (p[1] == a[1])
		{
			return 1;
		}
		else if (p[1] == b[1])
		{
			return 2;
		}
		else
		{
			return 3;
		}
	}
	else if (a[1] > b[1])
	{
		if ((p[1] > a[1]) || (p[1] < b[1]))
		{
			return 0;
		}
		else if (p[1] == a[1])
		{
			return 1;
		}
		else if (p[1] == b[1])
		{
			return 2;
		}
		else
		{
			return 3;
		}
	}
	else
	{
		fprintf(stderr, "FATAL ERROR: identical endpoints in InSegment() test\n");
		exit(1);
		return -1;
	}
}

bool TINMesh::IsOnSegment(const float* a, const float* b, const float* p)
{
	if (a[0] < b[0])
	{
		return ((a[0] <= p[0]) && (p[0] <= b[0]));
	}
	else if (a[0] > b[0])
	{
		return ((a[0] >= p[0]) && (p[0] >= b[0]));
	}
	else if (a[1] < b[1])
	{
		return ((a[1] <= p[1]) && (p[1] <= b[1]));
	}
	else if (a[1] > b[1])
	{
		return ((a[1] >= p[1]) && (p[1] >= b[1]));
	}
	else
	{
		fprintf(stderr, "FATAL ERROR: identical endpoints in IsOnSegment() test\n");
		exit(1);
		return false;
	}
}

bool TINMesh::IsInsideCircle(const TINtriangle* t, const float* p)
{
	if (t->V[0] == 0)  // infinite circle
	{
		double d = -orient2df(t->V[1], t->V[2], p);
		return (d < 0 ? true : (d > 0 ? false : IsOnSegment(t->V[1], t->V[2], p)));
	}
	else // finite circle
	{
		double d = -incirclef(t->V[0], t->V[1], t->V[2], p);
		return (d <= 0);
	}
}

bool TINMesh::IsValidTriangle(const TINtriangle* t)
{
	if (t->V[0])
		return true;
	else
		return false;
}

bool TINMesh::IsValidTriangle(int index)
{
	TINtriangle *t = GetTriangle(index);
	bool valid = IsValidTriangle(t);
	return valid;
}

int TINMesh::TriangleIndex(TINtriangle* t_locate)
{
	return t_locate - triangle_buffer;
}

int* TINMesh::PointIndices(TINtriangle* t_locate)
{
	int *indices = (int *)malloc(3 * sizeof(int));
	for (int i = 0; i < 3; ++i)
	{
		indices[i] = (t_locate->V[i] - point_buffer) / 3;
	}
	return indices;
}

TINtriangle* TINMesh::TINlocate_brute_query(float* p)
{
	int i;
	TINtriangle* t;
	double d[3];

	for (i = 0; i < triangle_buffer_alloc; i++)
	{
		t = triangle_buffer + i;
		if (t->next >= 0) continue;

		if (t->V[0] == 0) // infinite triangle?
		{
			if ((d[0] = orient2df(p, t->V[1], t->V[2])) >= 0)
			{
				if (d[0] == 0)
				{
					if ((d[0] = InSegment(t->V[1], t->V[2], p)) > 0)
					{
						if (d[0] == 3)
						{
							return t;
						}
						else
						{
							return t; // DUPLICATE_POINT
						}
					}
				}
				else
				{
					return t;
				}
			}
		}
		else
		{
			if (((d[0] = orient2df(p, t->V[1], t->V[2])) >= 0) &&
				((d[1] = orient2df(p, t->V[2], t->V[0])) >= 0) &&
				((d[2] = orient2df(p, t->V[0], t->V[1])) >= 0))
			{
				if (d[0] == 0)
				{
					if (d[1] == 0)
					{
						assert((p[0] == t->V[2][0]) && (p[1] == t->V[2][1]));
						return t; // DUPLICATE_POINT
					}
					else if (d[2] == 0)
					{
						assert((p[0] == t->V[1][0]) && (p[1] == t->V[1][1]));
						return t; // DUPLICATE_POINT
					}
				}
				else if (d[1] == 0)
				{
					if (d[2] == 0)
					{
						assert((p[0] == t->V[0][0]) && (p[1] == t->V[0][1]));
						return 0; // DUPLICATE_POINT
					}
				}
				return t;
			}
		}
	}
	// we should never get here
	fprintf(stderr, "FATAL ERROR: fail to brute_locate a point\n");
	exit(1);
	return 0;
}

TINtriangle* TINMesh::BrutePointLocation(float* p)
{
	int i;
	TINtriangle* t;
	double d[3];

	for (i = 0; i < triangle_buffer_alloc; i++)
	{
		t = triangle_buffer + i;
		if (t->next >= 0) continue;

		if (t->V[0] == 0) // infinite triangle?
		{
			if ((d[0] = orient2df(p, t->V[1], t->V[2])) >= 0)
			{
				if (d[0] == 0)
				{
					if ((d[0] = InSegment(t->V[1], t->V[2], p)) > 0)
					{
						if (d[0] == 3)
						{
							return t;
						}
						else
						{
							return 0; // DUPLICATE_POINT
						}
					}
				}
				else
				{
					return t;
				}
			}
		}
		else
		{
			if (((d[0] = orient2df(p, t->V[1], t->V[2])) >= 0) &&
				((d[1] = orient2df(p, t->V[2], t->V[0])) >= 0) &&
				((d[2] = orient2df(p, t->V[0], t->V[1])) >= 0))
			{
				if (d[0] == 0)
				{
					if (d[1] == 0)
					{
						assert((p[0] == t->V[2][0]) && (p[1] == t->V[2][1]));
						return 0; // DUPLICATE_POINT
					}
					else if (d[2] == 0)
					{
						assert((p[0] == t->V[1][0]) && (p[1] == t->V[1][1]));
						return 0; // DUPLICATE_POINT
					}
				}
				else if (d[1] == 0)
				{
					if (d[2] == 0)
					{
						assert((p[0] == t->V[0][0]) && (p[1] == t->V[0][1]));
						return 0; // DUPLICATE_POINT
					}
				}
				return t;
			}
		}
	}
	// we should never get here
	fprintf(stderr, "FATAL ERROR: fail to brute_locate a point\n");
	exit(1);
	return 0;
}

TINtriangle* TINMesh::SpecialPointLocation(const float* p, TINtriangle* t)
{
	do
	{
		assert(t->V[0] == 0);
		if (t->V[2][0] > t->V[1][0])
		{
			if (p[0] > t->V[2][0]) // is p beyond t->V[2]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[1]);
			}
			else if (p[0] < t->V[1][0]) // is p beyond t->V[1]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[2]);
			}
			else if (p[0] == t->V[2][0]) // is p a duplicate of t->V[2]
			{
				assert(p[1] == t->V[2][1]);
				return 0; // DUPLICATE_POINT
			}
			else if (p[0] == t->V[1][0]) // is p a duplicate of t->V[1]
			{
				assert(p[1] == t->V[1][1]);
				return 0; // DUPLICATE_POINT
			}
			else // found it
			{
				break;
			}
		}
		else if (t->V[2][0] < t->V[1][0])
		{
			if (p[0] < t->V[2][0]) // is p beyond t->V[2]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[1]);
			}
			else if (p[0] > t->V[1][0]) // is p beyond t->V[1]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[2]);
			}
			else if (p[0] == t->V[2][0]) // is p a duplicate of t->V[2]
			{
				assert(p[1] == t->V[2][1]);
				return 0; // DUPLICATE_POINT
			}
			else if (p[0] == t->V[1][0]) // is p a duplicate of t->V[1]
			{
				assert(p[1] == t->V[1][1]);
				return 0; // DUPLICATE_POINT
			}
			else // found it
			{
				break;
			}
		}
		else if (t->V[2][1] > t->V[1][1])
		{
			if (p[1] > t->V[2][1]) // is p beyond t->V[2]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[1]);
			}
			else if (p[1] < t->V[1][1]) // is p beyond t->V[1]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[2]);
			}
			else if (p[1] == t->V[2][1]) // is p a duplicate of t->V[2]
			{
				assert(p[0] == t->V[2][0]);
				return 0; // DUPLICATE_POINT
			}
			else if (p[1] == t->V[1][1]) // is p a duplicate of t->V[1]
			{
				assert(p[0] == t->V[1][0]);
				return 0; // DUPLICATE_POINT
			}
			else // found it
			{
				break;
			}
		}
		else if (t->V[2][1] < t->V[1][1])
		{
			if (p[1] < t->V[2][1]) // is p beyond t->V[2]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[1]);
			}
			else if (p[1] > t->V[1][1]) // is p beyond t->V[1]
			{
				t = triangle_buffer + TIN_TRIANGLE(t->N[2]);
			}
			else if (p[1] == t->V[2][1]) // is p a duplicate of t->V[2]
			{
				assert(p[0] == t->V[2][0]);
				return 0; // DUPLICATE_POINT
			}
			else if (p[1] == t->V[1][1])
			{
				assert(p[0] == t->V[1][0]); // is p a duplicate of t->V[1]
				return 0; // DUPLICATE_POINT
			}
			else // found it
			{
				break;
			}
		}
		else
		{
			fprintf(stderr, "FATAL ERROR: degenerate finite edge of infinite triangle\n");
			exit(1);
		}
	} while (orient2df(p, t->V[2], t->V[1]) == 0);

	assert(IsInsideCircle(t, p));
	return t;
}

TINtriangle* TINMesh::TINlocate_query(float* p)
{
	int steps = 0;
	double d;
	int c1;
	int c2;
	int ci;
	TINtriangle* t = triangle_buffer + triangle_newest;

	// are point p and triangle t on the same side of the oriented edge (v1,v2) 

	if ((d = orient2df(p, t->V[1], t->V[2])) >= 0)
	{
		// yes ... use this oriented edge as the starting point
		ci = 0;
	}
	else
	{
		// no ... try reversing the orientation of the edge using the adjacent triangle
		c1 = t->N[0];
		ci = TIN_CORNER(c1);
		t = triangle_buffer + TIN_TRIANGLE(c1);
	}

	c1 = TIN_NEXT(ci);
	c2 = TIN_PREV(ci);

	// early out?

	if (d == 0)
	{
		if ((steps = InSegment(t->V[1], t->V[2], p)))
		{
			if (steps == 3)
			{
				return t;
			}
			else
			{
				return t; // DUPLICATE_POINT
			}
		}
	}

	// compute q, the other end of the straight line we walk.
	double q[2];
	q[0] = ((double)(t->V[c2][0]) + (double)(t->V[c1][0])) / (double)(2.0);
	q[1] = ((double)(t->V[c2][1]) + (double)(t->V[c1][1])) / (double)(2.0);

	// walk a maximum of MAX_LOCATE_STEPS steps and maintain the edge (t, ci)

	while ((steps++) < MAX_LOCATE_STEPS)
	{
		assert(t->next == -1);

		if (t->V[0] == 0) // is the third vertex the infinite vertex?
		{
			if (d == 0)
			{
				// special handling because p is *on* the line defined by the finite edge of the infinite triangle
				// fprintf(stderr,"WARNING ... using special_locate\n");
				return SpecialPointLocation(p, t);
			}
			else
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
		else if ((d = orient2df(p, q, t->V[ci])) > 0) // does the ray from q to p pass vertex t->V[ci] on the left
		{
			// if p lies past the corresponding triangle edge we continue walking
			c1 = TIN_NEXT(ci);
			if ((d = orient2df(p, t->V[c1], t->V[ci])) > 0)
			{
				c2 = t->N[TIN_PREV(ci)];
				ci = TIN_CORNER(c2);
				t = triangle_buffer + TIN_TRIANGLE(c2);
			}
			else if (d == 0) // if p lies on the line through the corresponding triangle edge we check
			{
				if ((c2 = InSegment(t->V[c1], t->V[ci], p)))
				{
					if (c2 == 3)
					{
						assert(IsInsideCircle(t, p));
						return t;
					}
					else
					{
						return t; // DUPLICATE_POINT
					}
				}
				else
				{
					c2 = t->N[TIN_PREV(ci)];
					ci = TIN_CORNER(c2);
					t = triangle_buffer + TIN_TRIANGLE(c2);
				}
			}
			else // if p lies before the corresponding triangle edge we are done
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
		else // the ray from q to p pass passes t->V[ci] on the right or directly through it
		{
			// if the ray passes directly through t->V[ci] then p and t->V[ci] may be identical
			if (d == 0)
			{
				if (p[0] == t->V[ci][0] && p[1] == t->V[ci][1]) // p could be a duplicate of V[ci] 
				{
					return t; // DUPLICATE_POINT
				}
			}
			// if p lies past the corresponding triangle edge we continue walking
			c2 = TIN_PREV(ci);
			if ((d = orient2df(p, t->V[ci], t->V[c2])) > 0)
			{
				c1 = t->N[TIN_NEXT(ci)];
				ci = TIN_CORNER(c1);
				t = triangle_buffer + TIN_TRIANGLE(c1);
			}
			else if (d == 0) // if p lies on the line through the corresponding triangle edge we check
			{
				if ((c1 = InSegment(t->V[ci], t->V[c2], p)))
				{
					if (c1 == 3)
					{
						assert(IsInsideCircle(t, p));
						return t;
					}
					else
					{
						return t; // DUPLICATE_POINT
					}
				}
				else
				{
					c1 = t->N[TIN_NEXT(ci)];
					ci = TIN_CORNER(c1);
					t = triangle_buffer + TIN_TRIANGLE(c1);
				}
			}
			else // if p lies before the corresponding triangle edge we are done
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
	}
	fprintf(stderr, "WARNING ... using brute_locate\n");
	return BrutePointLocation(p);
}

TINtriangle* TINMesh::PointLocation(float* p)
{
	int steps = 0;
	double d;
	int c1;
	int c2;
	int ci;
	TINtriangle* t = triangle_buffer + triangle_newest;

	// are point p and triangle t on the same side of the oriented edge (v1,v2) 

	if ((d = orient2df(p, t->V[1], t->V[2])) >= 0)
	{
		// yes ... use this oriented edge as the starting point
		ci = 0;
	}
	else
	{
		// no ... try reversing the orientation of the edge using the adjacent triangle
		c1 = t->N[0];
		ci = TIN_CORNER(c1);
		t = triangle_buffer + TIN_TRIANGLE(c1);
	}

	c1 = TIN_NEXT(ci);
	c2 = TIN_PREV(ci);

	// early out?

	if (d == 0)
	{
		if ((steps = InSegment(t->V[1], t->V[2], p)))
		{
			if (steps == 3)
			{
				return t;
			}
			else
			{
				return 0; // DUPLICATE_POINT
			}
		}
	}

	// compute q, the other end of the straight line we walk.
	double q[2];
	q[0] = ((double)(t->V[c2][0]) + (double)(t->V[c1][0])) / (double)(2.0);
	q[1] = ((double)(t->V[c2][1]) + (double)(t->V[c1][1])) / (double)(2.0);

	// walk a maximum of MAX_LOCATE_STEPS steps and maintain the edge (t, ci)

	while ((steps++) < MAX_LOCATE_STEPS)
	{
		assert(t->next == -1);

		if (t->V[0] == 0) // is the third vertex the infinite vertex?
		{
			if (d == 0)
			{
				// special handling because p is *on* the line defined by the finite edge of the infinite triangle
				// fprintf(stderr,"WARNING ... using special_locate\n");
				return SpecialPointLocation(p, t);
			}
			else
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
		else if ((d = orient2df(p, q, t->V[ci])) > 0) // does the ray from q to p pass vertex t->V[ci] on the left
		{
			// if p lies past the corresponding triangle edge we continue walking
			c1 = TIN_NEXT(ci);
			if ((d = orient2df(p, t->V[c1], t->V[ci])) > 0)
			{
				c2 = t->N[TIN_PREV(ci)];
				ci = TIN_CORNER(c2);
				t = triangle_buffer + TIN_TRIANGLE(c2);
			}
			else if (d == 0) // if p lies on the line through the corresponding triangle edge we check
			{
				if ((c2 = InSegment(t->V[c1], t->V[ci], p)))
				{
					if (c2 == 3)
					{
						assert(IsInsideCircle(t, p));
						return t;
					}
					else
					{
						return 0; // DUPLICATE_POINT
					}
				}
				else
				{
					c2 = t->N[TIN_PREV(ci)];
					ci = TIN_CORNER(c2);
					t = triangle_buffer + TIN_TRIANGLE(c2);
				}
			}
			else // if p lies before the corresponding triangle edge we are done
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
		else // the ray from q to p pass passes t->V[ci] on the right or directly through it
		{
			// if the ray passes directly through t->V[ci] then p and t->V[ci] may be identical
			if (d == 0)
			{
				if (p[0] == t->V[ci][0] && p[1] == t->V[ci][1]) // p could be a duplicate of V[ci] 
				{
					return 0; // DUPLICATE_POINT
				}
			}
			// if p lies past the corresponding triangle edge we continue walking
			c2 = TIN_PREV(ci);
			if ((d = orient2df(p, t->V[ci], t->V[c2])) > 0)
			{
				c1 = t->N[TIN_NEXT(ci)];
				ci = TIN_CORNER(c1);
				t = triangle_buffer + TIN_TRIANGLE(c1);
			}
			else if (d == 0) // if p lies on the line through the corresponding triangle edge we check
			{
				if ((c1 = InSegment(t->V[ci], t->V[c2], p)))
				{
					if (c1 == 3)
					{
						assert(IsInsideCircle(t, p));
						return t;
					}
					else
					{
						return 0; // DUPLICATE_POINT
					}
				}
				else
				{
					c1 = t->N[TIN_NEXT(ci)];
					ci = TIN_CORNER(c1);
					t = triangle_buffer + TIN_TRIANGLE(c1);
				}
			}
			else // if p lies before the corresponding triangle edge we are done
			{
				assert(IsInsideCircle(t, p));
				return t;
			}
		}
	}
	fprintf(stderr, "WARNING ... using brute_locate\n");
	return BrutePointLocation(p);
}

void TINMesh::Update(float* p, TINtriangle* t)
{
	// the circumcircle of triangle tri_idx must enclose the new point p.
	int c;
	int ci;

	assert(stkEMPTY());

	// add edges of seed triangle in CCW order to the stack
	for (ci = 0; ci < 3; ci++)
	{
		stkPUSH(t->N[ci]);
	}

	// delete the triangle
	t->next = triangle_next;
	triangle_next = (int)(t - triangle_buffer);

	int prev_idx = -1;
	int prev_ci = -1;
	int last_idx = -1;
	int last_ci = -1;

	while (!stkEMPTY())
	{
		// get a corner from the horizon search stack
		c = stkPOP();

		ci = TIN_CORNER(c);
		t = triangle_buffer + TIN_TRIANGLE(c);

		assert(t->next < 0);

		if (IsInsideCircle(t, p)) // if the circumcircle of t encloses p we need to expand the horizon
		{
			// add the two unchecked edges of this triangle in CCW order to the stack
			for (c = 0; c < 2; c++)
			{
				ci = TIN_NEXT(ci);
				stkPUSH(t->N[ci]);
			}
			// delete the triangle
			t->next = triangle_next;
			triangle_next = (int)(t - triangle_buffer);
		}
		else // we found a horizon edge. create the adjacent new triangle.
		{
			int tn_idx = triangle_next;
			TINtriangle* tn = triangle_buffer + tn_idx;
			triangle_newest = tn_idx;
			triangle_next = tn->next;
			tn->next = -1;
			if (triangle_buffer_size < triangle_next) triangle_buffer_size = triangle_next;

			// one of the vertices on the horizon could be the infinite vertex ... make sure it is stored in tn->V[0] 
			if (ci == 1)
			{
				// assign vertices with correct orientation
				tn->V[0] = t->V[0];
				tn->V[1] = t->V[2];
				tn->V[2] = p;
				// link the triangles across the horizon edge
				tn->N[2] = c;
				t->N[1] = TIN_INDEX(tn_idx, 2);
				// link the new triangles around the new vertex
				if (prev_idx != -1) // does a previous triangle exist?
				{
					tn->N[0] = TIN_INDEX(prev_idx, prev_ci); // tn->N[0] should point to previous
					(triangle_buffer + prev_idx)->N[prev_ci] = TIN_INDEX(tn_idx, 0);
				}
				else // no ... this pointer will be set at the end
				{
					last_idx = tn_idx;
					last_ci = 0;
				}
				prev_idx = tn_idx;
				prev_ci = 1;
			}
			else if (ci == 2)
			{
				// assign vertices with correct orientation
				tn->V[0] = t->V[0];
				tn->V[1] = p;
				tn->V[2] = t->V[1];
				// link the triangles across the horizon edge
				tn->N[1] = c;
				t->N[2] = TIN_INDEX(tn_idx, 1);
				// link the new triangles around the new vertex
				if (prev_idx != -1) // does a previous triangle exist?
				{
					tn->N[2] = TIN_INDEX(prev_idx, prev_ci); // tn->N[2] should point to previous
					(triangle_buffer + prev_idx)->N[prev_ci] = TIN_INDEX(tn_idx, 2);
				}
				else // no ... this pointer will be set at the end
				{
					last_idx = tn_idx;
					last_ci = 2;
				}
				prev_idx = tn_idx;
				prev_ci = 0;
			}
			else
			{
				// assign vertices with correct orientation
				tn->V[0] = p;
				tn->V[1] = t->V[2];
				tn->V[2] = t->V[1];
				// link the triangles across the horizon edge
				tn->N[0] = c;
				t->N[0] = TIN_INDEX(tn_idx, 0);
				// link the new triangles around the new vertex
				if (prev_idx != -1) // does a previous triangle exist?
				{
					tn->N[1] = TIN_INDEX(prev_idx, prev_ci); // tn->N[1] should point to previous
					(triangle_buffer + prev_idx)->N[prev_ci] = TIN_INDEX(tn_idx, 1);
				}
				else // no ... this pointer will be set at the end
				{
					last_idx = tn_idx;
					last_ci = 1;
				}
				prev_idx = tn_idx;
				prev_ci = 2;
			}
		}
	}
	// link the last two triangles
	assert((prev_idx != -1) && (last_idx != -1));
	(triangle_buffer + last_idx)->N[last_ci] = TIN_INDEX(prev_idx, prev_ci);
	(triangle_buffer + prev_idx)->N[prev_ci] = TIN_INDEX(last_idx, last_ci);
}

bool TINMesh::InsertPoint(float* p)
{
	// find a triangle whose circumcircle contains p
	TINtriangle* t = PointLocation(p);
	// discard duplicate points 
	if (t == 0)
	{
		return false;
	}
	// do we only care about the convex hull
	if (convex_hull_only)
	{
		// then we only insert points that fall into an infinite triangle
		if (t->V[0])
		{
			return false;
		}
	}
	// update triangulation with Bowyer-Watson
	Update(p, t);
	return true;
}

int TINMesh::Size()
{
	return triangle_buffer_size;
}

TINtriangle* TINMesh::GetTriangle(int t)
{
	return triangle_buffer + t;
}

void TINMesh::Initialize(int num, bool convex_hull)
{
	// initialize jonathan's predicate
	exactinit();
	// initialize triangle buffer
	if (triangle_buffer_alloc < num * 2)
	{
		triangle_buffer_alloc = num * 2;
		if (triangle_buffer) free(triangle_buffer);
		triangle_buffer = (TINtriangle*)malloc(sizeof(TINtriangle)*triangle_buffer_alloc);
		if (triangle_buffer == 0)
		{
			fprintf(stderr, "ERROR (TINInitial): failed malloc for %d TINtriangles\n", triangle_buffer_alloc);
			exit(1);
		}
	}
	for (int i = 0; i < triangle_buffer_alloc; i++)
	{
		triangle_buffer[i].next = i + 1; // points to next unused triangle
	}
	triangle_buffer_size = 0;
	pointer = 0;
	buffer_full = false;
	initialized = false;
	convex_hull_only = convex_hull;
}

void TINMesh::Derive()
{
	int i;
	if (buffer_full)
	{
		for (i = 0; i < POINT_BRIO_BUFFER; i++) InsertPoint(points[i]);
	}
	else
	{
		for (i = 0; i < pointer; i++) InsertPoint(points[i]);
	}
}

void TINMesh::Destory()
{
	if (triangle_buffer) free(triangle_buffer);
	triangle_buffer = NULL;
	stkDESTROY();
}