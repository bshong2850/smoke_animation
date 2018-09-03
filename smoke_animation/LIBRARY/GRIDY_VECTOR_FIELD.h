/*
 *	GRIDY_VECTOR_FIELD
 *
 *	현재는 3차원만 제공한다.
 *
*/

#ifndef _GRIDY_VECTOR_FIELD_H_
#define _GRIDY_VECTOR_FIELD_H_





#include <omp.h>
#include <math.h>
//#include "GRIDY_COMMON.h"



class GRIDY_VECTOR_FIELD
{
private:
	int mWidth, mHeight, mDepth;
	int mSize;


public:
	double *u, *v, *w;

	GRIDY_VECTOR_FIELD(int width, int height, int depth);
	GRIDY_VECTOR_FIELD() {};
	~GRIDY_VECTOR_FIELD();

	void Clear_data(double default_value = 0.0f);
	//void Clear_data_(double default_value = 0.0f);
	void Init(int width, int height, int depth);

	double Get_norm2(short iu, short iv, short iw);
	double Sum_data(void);

	int inline Get_size(void) { return mSize; }
	int inline Get_width_size(void) { return mWidth; }
	int inline Get_height_size(void) { return mHeight; }
	int inline Get_depth_size(void) { return mDepth; }
};





GRIDY_VECTOR_FIELD::GRIDY_VECTOR_FIELD(int width, int height, int depth) 
{
	Init(width, height, depth);
}



GRIDY_VECTOR_FIELD::~GRIDY_VECTOR_FIELD() 
{
	if(u) delete[] u;		//	free array
	if(v) delete[] v;
	if(w) delete[] w;
}



void GRIDY_VECTOR_FIELD::Init(int width, int height, int depth)
{
	//	init scalar field
	mSize = (width+2)*(height+2)*(depth+2);

	mWidth = width;
	mHeight = height;
	mDepth = depth;

	//	allocate array
	u = new double[mSize];
	v = new double[mSize];
	w = new double[mSize];

	Clear_data();
}



void GRIDY_VECTOR_FIELD::Clear_data(double default_value)
{
	#pragma omp parallel for
	for(int i=0; i<mSize; i++)
	{
		u[i] = v[i] = w[i] = default_value;
	}
}



double GRIDY_VECTOR_FIELD::Get_norm2(short iu, short iv, short iw)
{ 
	//double DoDx = (NORM2(vortex.u[IX3(i+1,j,k)],vortex.v[IX3(i,j,k)],vortex.w[IX3(i,j,k)]) - NORM2(vortex.u[IX3(i-1,j,k)],vortex.v[IX3(i,j,k)],vortex.w[IX3(i,j,k)])) * 0.5f;
	double vu, vv, vw, sqv;
	
	vu = u[(iu) + (mWidth+2)*(iv) + (mWidth+2)*(mHeight+2)*(iw)];
	vv = v[(iu) + (mWidth+2)*(iv) + (mWidth+2)*(mHeight+2)*(iw)];
	vw = w[(iu) + (mWidth+2)*(iv) + (mWidth+2)*(mHeight+2)*(iw)];
	sqv = sqrt(vu*vu + vv*vv + vw*vw);

	return sqv;
}




#endif // _GRIDY_VECTOR_FIELD_H_