#ifndef __VECTOR_3F_H__
#define __VECTOR_3F_H__

#include <algorithm>
#include <vector>
#include <math.h>
using namespace std;
#define E 0.0001
#define VECTOR3F_PRINTF(v3) ((v3).x), ((v3).y), ((v3).z)

struct vector3f{
	double x, y, z;
	vector3f(){ x = y = z = 0; }
	vector3f(double p){ x = y = z = p; }
	vector3f(double nx, double ny, double nz){ x = nx; y = ny; z = nz; }
	friend vector3f operator +(vector3f p, vector3f q)
	{
		vector3f r = vector3f(p.x + q.x, p.y + q.y, p.z + q.z);
		return r;
	}
	friend vector3f operator *(double p, vector3f q)
	{
		vector3f r = vector3f(p*q.x, p*q.y, p*q.z);
		return r;
	}
	friend vector3f operator *(vector3f q, double p)
	{
		vector3f r = vector3f(p*q.x, p*q.y, p*q.z);
		return r;
	}
	friend vector3f operator /(vector3f q, double p)
	{
		vector3f r = vector3f(q.x / p, q.y / p, q.z / p);
		return r;
	}
	friend double operator *(vector3f p, vector3f q)
	{
		double r = (p.x*q.x + p.y*q.y + p.z*q.z);
		return r;
	}
	friend vector3f operator -(vector3f p, vector3f q)
	{
		vector3f r = vector3f(p.x - q.x, p.y - q.y, p.z - q.z);
		return r;
	}
	friend vector3f operator %(vector3f p, vector3f q)
	{
		vector3f r = vector3f(p.y*q.z - p.z*q.y, p.z*q.x - p.x*q.z, p.x*q.y - p.y*q.x);
		return r;
	}
	friend vector3f operator ^(vector3f p, vector3f q)
	{
		vector3f r = vector3f(p.x*q.x, p.y*q.y, p.z*q.z);
		return r;
	}
	void operator += (const vector3f &p)			{ x += p.x;  y += p.y; z += p.z; }
	void operator -= (const vector3f &p)			{ x -= p.x;  y -= p.y; z -= p.z; }
	void operator *= (double c)						{ x *= c;  y *= c; z *= c; }
	void operator /= (double c)						{ x /= c;  y /= c; z /= c; }
	void operator ^= (const vector3f &p)			{ x *= p.x;  y *= p.y; z *= p.z; }

	double dist()
	{
		return sqrt(x*x + y*y + z*z);
	}
	double dist(vector3f p)
	{
		vector3f r = *this - p;
		return r.dist();
	}
	void normalize()
	{
		double w = dist();
		if (w < E) return;
		x /= w;
		y /= w;
		z /= w;
	}
};

struct matrix33f{
	double a[3][3];
	matrix33f(){ for (int i = 0; i<3; i++)for (int j = 0; j<3; j++)a[i][j] = 0; }

	friend matrix33f operator *(matrix33f p, matrix33f q)
	{
		int i, j, k;
		matrix33f r;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				for (k = 0; k<3; k++)
					r.a[i][j] += p.a[i][k] * q.a[k][j];
		return r;
	}
	friend matrix33f operator *(matrix33f p, double q)
	{
		int i, j;
		matrix33f r;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				r.a[i][j] = p.a[i][j] * q;
		return r;
	}
	friend matrix33f operator /(matrix33f p, double q)
	{
		int i, j;
		matrix33f r;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				r.a[i][j] = p.a[i][j] / q;
		return r;
	}

	friend matrix33f operator +(matrix33f p, matrix33f q)
	{
		int i, j;
		matrix33f r;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				r.a[i][j] = p.a[i][j] + q.a[i][j];
		return r;
	}
	friend matrix33f operator -(matrix33f p, matrix33f q)
	{
		int i, j;
		matrix33f r;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				r.a[i][j] = p.a[i][j] - q.a[i][j];
		return r;
	}
	friend vector3f operator *(matrix33f p, vector3f q)
	{
		vector3f r;
		r.x += p.a[0][0] * q.x; r.x += p.a[0][1] * q.y; r.x += p.a[0][2] * q.z;
		r.y += p.a[1][0] * q.x; r.y += p.a[1][1] * q.y; r.y += p.a[1][2] * q.z;
		r.z += p.a[2][0] * q.x; r.z += p.a[2][1] * q.y; r.z += p.a[2][2] * q.z;
		return r;
	}

	void operator +=(matrix33f p) { for (int i = 0; i<3; i++)for (int j = 0; j<3; j++) a[i][j] += p.a[i][j]; }
	void operator -=(matrix33f p) { for (int i = 0; i<3; i++)for (int j = 0; j<3; j++) a[i][j] -= p.a[i][j]; }
	void operator /=(double p){ for (int i = 0; i<3; i++)for (int j = 0; j<3; j++)a[i][j] /= p; }

	double invert(matrix33f &r)
	{
		int i, j;
		for (i = 0; i<3; i++)
			for (j = 0; j<3; j++)
				r.a[j][i] = a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3] - a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3];
		double det = 0;
		for (i = 0; i<3; i++)
			det += r.a[i][0] * a[0][i];
		r = r / det;
		return det;
	}
};

struct vector2f{
	double x, y;
	vector2f(){ x = y = 0; }
	vector2f(double p){ x = y = p; }
	vector2f(double nx, double ny){ x = nx; y = ny; }
	vector2f(int nx, int ny){ x = nx; y = ny; }
	friend vector2f operator +(vector2f p, vector2f q)
	{
		vector2f r = vector2f(p.x + q.x, p.y + q.y);
		return r;
	}
	friend vector2f operator *(double p, vector2f q)
	{
		vector2f r = vector2f(p*q.x, p*q.y);
		return r;
	}
	friend vector2f operator *(vector2f q, double p)
	{
		vector2f r = vector2f(p*q.x, p*q.y);
		return r;
	}
	friend vector2f operator /(vector2f q, double p)
	{
		vector2f r = vector2f(q.x / p, q.y / p);
		return r;
	}
	friend double operator *(vector2f p, vector2f q)
	{
		double r = (p.x*q.x + p.y*q.y);
		return r;
	}
	friend vector2f operator -(vector2f p, vector2f q)
	{
		vector2f r = vector2f(p.x - q.x, p.y - q.y);
		return r;
	}
	void operator += (const vector2f &p)			{ x += p.x;  y += p.y; }
	void operator -= (const vector2f &p)			{ x -= p.x;  y -= p.y; }
	void operator *= (double c)						{ x *= c;  y *= c; }
	void operator /= (double c)						{ x /= c;  y /= c; }

	double dist()
	{
		return sqrt(x*x + y*y);
	}
	double dist(vector2f p)
	{
		vector2f r = *this - p;
		return r.dist();
	}
	void normalize()
	{
		double w = dist();
		if (w < E) return;
		x /= w;
		y /= w;
	}
}; 
#endif /* __VECTOR_3F_H__ */