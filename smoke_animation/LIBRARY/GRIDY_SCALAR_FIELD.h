/*
 *	GRIDY_SCALAR_FIELD
 *
 *	현재로서는 3차원 스칼라 필드만 제공한다.
 *
*/

#ifndef _GRIDY_SCALAR_FIELD_H_
#define _GRIDY_SCALAR_FIELD_H_





//#include "GRIDY_GRID.h"
#include <omp.h>



class GRIDY_SCALAR_FIELD
{
protected:
	int mWidth, mHeight, mDepth;
	int mSize;


public:
	double *value;

	GRIDY_SCALAR_FIELD(int width, int height, int depth);
	GRIDY_SCALAR_FIELD() {
		mWidth = mHeight = mDepth = -1;
	};
	~GRIDY_SCALAR_FIELD();

	void Clear_data(double default_value = 0.0f);
	void Init(int width, int height, int depth);

	int inline Get_size(void) { return mSize; }

	int inline Get_width_size(){ return mWidth; }
	int inline Get_height_size(){ return mHeight; }
	int inline Get_depth_size(){ return mDepth; }
};





GRIDY_SCALAR_FIELD::GRIDY_SCALAR_FIELD(int width, int height, int depth) 
{
	Init(width, height, depth);
}



GRIDY_SCALAR_FIELD::~GRIDY_SCALAR_FIELD() 
{
	if(value) delete[] value;		//	free array
}



void GRIDY_SCALAR_FIELD::Init(int width, int height, int depth)
{
	//	init scalar field
	mSize = (width+2)*(height+2)*(depth+2);

	mWidth = width;
	mHeight = height;
	mDepth = depth;

	//	allocate array
	value = new double[mSize];

	Clear_data();
}



void GRIDY_SCALAR_FIELD::Clear_data(double default_value)
{
	#pragma omp parallel for
	for(int i=0; i<mSize; i++)
	{
		value[i] = default_value;
	}
}


#endif // _GRIDY_SCALAR_FIELD_H_