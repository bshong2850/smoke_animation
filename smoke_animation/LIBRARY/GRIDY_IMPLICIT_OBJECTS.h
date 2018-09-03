#ifndef _GRIDY_IMPLICIT_OBJECTS_H_
#define _GRIDY_IMPLICIT_OBJECTS_H_

#include <math.h>
#include "GRIDY_COMMON.h"
#include "GRIDY_ADT.h"
#include "GRIDY_SCALAR_FIELD.h"
//#include "GRIDY_LOG.h"





GRIDY_COMMON mUTIL;





//	implicit 3D (or 2D) sphere
class GRIDY_IMPLICIT_OBJECT_SPHERE
{
public:
	double center_x;
	double center_y;
	double center_z;
	double radius;


public:
	void Init(double mCenter_x, double mCenter_y, double mRadius)	// 2D circle
	{
		center_x = mCenter_x;
		center_y = mCenter_y;
		center_z = -1.0f;
		radius = mRadius;
	}

	void Init(double mCenter_x, double mCenter_y, double mCenter_z, double mRadius)	// 3D sphere
	{
		center_x = mCenter_x;
		center_y = mCenter_y;
		center_z = mCenter_z;
		radius = mRadius;
	}

	int Levelset(int x, int y)	// for 2D
	{
		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2));

		if (dist <= radius) return -1;
		else if (abs(dist - radius) <= 1.0f) return 0;
		else return 1;
	}

	int Levelset(int x, int y, int z)	// for 3D
	{
		int r = 1;

		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

		if (dist < radius) r = -1;
		if (abs(dist - radius) <= 1.74f) r = 0;

		return r;
	}

	double Get_distance(int x, int y, int z)
	{
		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

		return abs(dist - radius);
	}


	double Get_distance(double x, double y, double z)
	{
		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

		return abs(dist - radius);
	}



	void Get_normal(double *normal, int x, int y, int z)
	{
		double norm;

		if (center_x < 0 || center_y < 0 || center_z < 0)
		{
			normal[0] = normal[1] = normal[2] = 0.0f;
		}
		else
		{
			norm = sqrt(((double)x - center_x)*((double)x - center_x) + ((double)y - center_y)*((double)y - center_y) + ((double)z - center_z)*((double)z - center_z));

			normal[0] = (x - center_x) / norm;
			normal[1] = (y - center_y) / norm;
			normal[2] = (z - center_z) / norm;
		}
	}



	void Get_normal(double *normal, double x, double y, double z)
	{
		double norm;

		if (center_x < 0 || center_y < 0 || center_z < 0)
		{
			normal[0] = normal[1] = normal[2] = 0.0f;
		}
		else
		{
			norm = sqrt(((double)x - center_x)*((double)x - center_x) + ((double)y - center_y)*((double)y - center_y) + ((double)z - center_z)*((double)z - center_z));

			normal[0] = (x - center_x) / norm;
			normal[1] = (y - center_y) / norm;
			normal[2] = (z - center_z) / norm;

			//printf("!!%f %f %f\n", normal[0], normal[1], normal[2]);
		}
	}
};










// Flame shape control에서 사용하는 class
class GRIDY_IMPLICIT_OBJECT_SPHERES
{
public:
	int mNo;			// 총 파티클 갯수.
	double mMax_radius;	// 파티클이 가질 수 있는 최대 반지름

	double *mCx;			// 현재 파티클 center x
	double *mCy;			// 현재 파티클 center y
	double *mCz;			// 현재 파티클 center z
	double *mNx;			// 현재 파티클 normal x
	double *mNy;			// 현재 파티클 normal y
	double *mNz;			// 현재 파티클 normal z
	double *mAngle;			// 현재 파티클 회전각도
	double *mRadius;		// 현재 파티클 반지름
	double *mTemperature;	// 현재 파티클 온도

	bool mOption_Vibrate; // 진동 모드인가?

	bool *isActivate;	// 현재 동작하고 있나?
	bool *isGrowing;	// 현재 자라고 있나?
	bool *isVibrate;	// 진동 모드인가?



public:
	GRIDY_IMPLICIT_OBJECT_SPHERES()
	{
		mNo = -1;
		mMax_radius = 0;
		mOption_Vibrate = false;
	}


	~GRIDY_IMPLICIT_OBJECT_SPHERES()
	{
		Free();
	}


	GRIDY_IMPLICIT_OBJECT_SPHERES(int n, double max_radius)
	{
		mNo = -1;
		mMax_radius = 0;
		mOption_Vibrate = false;

		Init(n, max_radius);
	}


	void Init(int n, double max_radius = 5.0f)
	{
		if (n < 0)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR.\n";
			exit(0);
		}
		else
		{
			mNo = n;

			mCx = new double[mNo];
			mCy = new double[mNo];
			mCz = new double[mNo];
			mNx = new double[mNo];
			mNy = new double[mNo];
			mNz = new double[mNo];
			mAngle = new double[mNo];
			mRadius = new double[mNo];
			mTemperature = new double[mNo];
			isActivate = new bool[mNo];
			isGrowing = new bool[mNo];
			isVibrate = new bool[mNo];
		}

		mMax_radius = max_radius;
	}





	void Free(void)
	{
		if (mCx) delete[] mCx;
		if (mCy) delete[] mCy;
		if (mCz) delete[] mCz;
		if (mRadius) delete[] mRadius;
		if (mNx) delete[] mNx;
		if (mNy) delete[] mNy;
		if (mNz) delete[] mNz;
		if (mAngle) delete[] mAngle;
		if (mTemperature) delete[] mTemperature;
		if (isActivate) delete[] isActivate;
		if (isGrowing) delete[] isGrowing;
	}



	void Set(int n, double center_x, double center_y, double center_z, double radius)
	{
		if (n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Set).\n";
			cout << "n is " << n << ", and mNo is " << mNo << endl;
			exit(0);
		}
		else
		{
			mCx[n] = center_x;
			mCy[n] = center_y;
			mCz[n] = center_z;
			mRadius[n] = radius;
			isActivate[n] = isGrowing[n] = isVibrate[n] = false;
		}
	}



	int Levelset(int n, int x, int y, int z)
	{
		if (n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Levelset).\n";
			exit(0);
		}
		else
		{
			int levelset = 1;

			double dist = sqrt((mCx[n] - x)*(mCx[n] - x) + (mCy[n] - y)*(mCy[n] - y) + (mCz[n] - z)*(mCz[n] - z));

			if (dist < mRadius[n]) levelset = -1;
			else if (abs(dist - mRadius[n]) < 1.74f) levelset = 0;

			return levelset;
		}
	}



	double Get_distance(int n, int x, int y, int z)
	{
		if (n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Get_distance).\n";
			exit(0);
		}
		else
		{
			double dist = sqrt(pow((mCx[n] - x), 2) + pow((mCy[n] - y), 2) + pow((mCz[n] - z), 2));

			return abs(dist - mRadius[n]);
		}
	}



	void Get_normal(int n, int x, int y, int z, double *normal)
	{
		if (n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Get_normal).\n";
			exit(0);
		}
		else
		{
			double dx = (double)x - mCx[n];
			double dy = (double)y - mCy[n];
			double dz = (double)z - mCz[n];
			double norm = sqrt(dx*dx + dy*dy + dz*dz);

			if (norm != 0.0f)
			{
				normal[0] = dx / norm;
				normal[1] = dy / norm;
				normal[2] = dz / norm;
			}
			else
			{
				//				printf("NORM ERROR : %f %f %f - %d %d %d\n", mCx[n], mCy[n], mCz[n], x, y, z);
				normal[0] = normal[1] = normal[2] = 0.0f;
			}
		}
	}




	void Growing_radius(int n, double coef_grow = 1.0f)
	{
		if (n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Set).\n";
			exit(0);
		}
		else if (mRadius[n] < 1.0f || mCx[n] < 0.0f || mCy[n] < 0.0f || mCz[n] < 0.0f)
		{
			isActivate[n] = false;
		}
		else
		{
			if (isGrowing[n] == true && isActivate[n] == true)
			{
				if (isVibrate[n] == true && mOption_Vibrate == true)
				{
					mRadius[n] -= coef_grow;
					if (mRadius[n] < mMax_radius / 2.0f)
					{
						mRadius[n] = mMax_radius / 2.0f;
						isVibrate[n] = false;
					}
				}
				else
				{
					mRadius[n] += coef_grow;
					if (mRadius[n] > mMax_radius)
					{
						mRadius[n] = mMax_radius;
						isVibrate[n] = true;
					}
				}
			}
		}
	}

};










class GRIDY_IMPLICIT_OBJECT_MOVING_SPHERES : GRIDY_IMPLICIT_OBJECT_SPHERES
{
public:
	struct str_GRIDY_MOVING_SPHERES {
		double x, y, z;
		double vx, vy, vz;
		double ax, ay, az;
		double r;
	};

	double mDef_r;
	double mDef_vx, mDef_vy, mDef_vz;


	//str_GRIDY_MOVING_SPHERES *mSpheres;
	//mAmount = 0;
	list <str_GRIDY_MOVING_SPHERES> mSpheres;



public:
	GRIDY_IMPLICIT_OBJECT_MOVING_SPHERES(void)	{ };
	~GRIDY_IMPLICIT_OBJECT_MOVING_SPHERES(void)	{ };


	void Insert(double r, double x, double y, double z, double vx, double vy, double vz, double ax, double ay, double az)
	{
		str_GRIDY_MOVING_SPHERES t;

		t.r = r;
		t.x = x;		t.y = y;		t.z = z;
		t.vx = vx;		t.vy = vy;		t.vz = vz;
		t.ax = ax;		t.ay = ay;		t.az = az;

		mSpheres.push_back(t);
	}

	/*
	void Init(double def_r, double def_vx, double def_vy, double def_vz)
	{
	mDef_r = def_r;
	mDef_vx = def_r;
	mDef_vy = def_r;
	mDef_vz = def_r;
	}
	*/
	/*
	void Init(double mCenter_x, double mCenter_y, double mRadius)
	{
	center_x = mCenter_x;
	center_y = mCenter_y;
	radius = mRadius;
	}

	void Init(double mCenter_x, double mCenter_y, double mCenter_z, double mRadius)
	{
	center_x = mCenter_x;
	center_y = mCenter_y;
	center_z = mCenter_z;
	radius = mRadius;
	}


	int Levelset(int x, int y, int z)
	{
	int r = 1;

	double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

	if (dist < radius) r = -1;
	if (abs(dist - radius) <= 1.74f) r = 0;
	//if (abs(dist - radius) <= 0.5f) r = 0;

	return r;
	}

	double Get_distance(int x, int y, int z)
	{
	double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

	return abs(dist - radius);
	}


	void Get_normal(double *normal, int x, int y, int z)
	{
	double norm;

	if (center_x < 0 || center_y < 0 || center_z < 0)
	{
	normal[0] = normal[1] = normal[2] = 0.0f;
	}
	else
	{
	norm = sqrt(((double)x - center_x)*((double)x - center_x) + ((double)y - center_y)*((double)y - center_y) + ((double)z - center_z)*((double)z - center_z));

	normal[0] = (x - center_x) / norm;
	normal[1] = (y - center_y) / norm;
	normal[2] = (z - center_z) / norm;
	}
	}
	*/
};




#endif