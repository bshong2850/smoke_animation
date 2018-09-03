#ifndef _GRIDY_IMPLICIT_OBJECT_H_
#define _GRIDY_IMPLICIT_OBJECT_H_



#pragma warning(disable:4996)


#include <math.h>
#include "GRIDY_COMMON.h"
#include "GRIDY_ADT.h"
#include "GRIDY_SCALAR_FIELD.h"

//GRIDY_COMMON mUTIL;
//GRIDY_LOG mLOG;



class GRIDY_IMPLICIT_OBJECT
{
private:

public:
	GRIDY_IMPLICIT_OBJECT(void)	{};
	~GRIDY_IMPLICIT_OBJECT(void) {};

	


protected :
	int mWidth, mHeight, mDepth;

	inline double Phi_on_grid(double vertex_pos, double start_AABB, double half_AABB_size, int grid_edge_size)
	{
		return (double)((abs(start_AABB-vertex_pos)*grid_edge_size)/(half_AABB_size*2.0f));
	}

	void Fill_density_between_two_vertices(double x0, double y0, double z0, double x1, double y1, double z1, double *arr, double mDefValue = 0.1f);
	//void 



public:
	
	void Put_obj_to_grid(double *arr, char* full_file_name, int width, int height, int depth, double scale = 0.5f, double mDefValue = 0.1f, int margin_x = 0, int margin_y = 0, int margin_z = 0)
	{
		mWidth = width;
		mHeight = height;
		mDepth = depth;
		//int mCount = 0;
	
		char mId;
		char buf_to_drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;
		//int* face_idx[3];

		double mAA[3], mBB[3];			// boundary of obj
		//double scale = 0.5f;

		errno_t err;

		FILE *fp;

		struct str_vertex {
			double x; double y; double z;
		};

		struct str_face_idx {
			int idx0; int idx1; int idx2;
		};

		std::vector<str_vertex> mVtx;
		std::vector<str_face_idx> mFace;

		//	init arr mAA, mBB
		for(int i=0; i<3; i++)
		{
			mAA[i] = 9999999.0f;
			mBB[i] = -9999999.0f;
		}
			
		if((err = fopen_s(&fp, full_file_name, "r"))!= 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}
		
		// set mAA, mBB
		while(!feof(fp))
		{
			//fscanf(fp, "%c %lf %lf %lf", &mId, &mVertex_x, &mVertex_y, &mVertex_z);
			fscanf_s(fp, "%c %lf %lf %lf", &mId, 1, &mVertex_x, &mVertex_y, &mVertex_z);
			//fscanf_s(fp, )
			
			if(mId == '#')	// have to drop next strings
			{
				fgets(buf_to_drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				if (abs(mVertex_x) < 99999999.9f && abs(mVertex_y) < 99999999.9f && abs(mVertex_z) < 99999999.9f) 
				{
					// make mAA get the smallest position
					if(mAA[0] > mVertex_x) mAA[0] = mVertex_x;
					if(mAA[1] > mVertex_y) mAA[1] = mVertex_y;
					if(mAA[2] > mVertex_z) mAA[2] = mVertex_z;

					// make mBB get the largest position
					if(mBB[0] < mVertex_x) mBB[0] = mVertex_x;
					if(mBB[1] < mVertex_y) mBB[1] = mVertex_y;
					if(mBB[2] < mVertex_z) mBB[2] = mVertex_z;

					// vector container
					str_vertex temp_vtx = { mVertex_x, mVertex_y, mVertex_z };
					mVtx.push_back(temp_vtx);
				}
			}
			else if (mId == 'f')
			{
				//str_vertex temp_face = { mVertex_x, mVertex_y, mVertex_z };
				//mFace.push_back([);
				mFace.push_back({ (int)mVertex_x, (int)mVertex_y, (int)mVertex_z });
				//mFace.data;
			}
		}

		// reload the obj file
		rewind(fp);
		fclose(fp);

		// set half_AABB_size
		double edge_size_width = (double)(mBB[0]-mAA[0]);
		double edge_size_height = (double)(mBB[1]-mAA[1]);
		double edge_size_depth = (double)(mBB[2]-mAA[2]);

		double max_edge_size = MAX(MAX(edge_size_width, edge_size_height), edge_size_depth);
		
		int max_edge_res = 0;
		if (max_edge_size <= edge_size_width) max_edge_res = mWidth;
		else if (max_edge_size <= edge_size_height) max_edge_res = mHeight;
		else max_edge_res = mDepth;

		double delta = max_edge_size / max_edge_res;

		//double max_edge = MAX(edge_size_width, edge_size_height)

		// 이게 잘못 되었다! 아닌가?
		//double delta_width = (double)abs(mBB[0]-mAA[0])/mWidth;
		//double delta_height = (double)abs(mBB[1]-mAA[1])/mHeight;
		//double delta_depth = (double)abs(mBB[2]-mAA[2])/mDepth;


		//printf("%lf %lf %lf \t %lf %lf %lf \n", edge_size_width, edge_size_height, edge_size_depth, delta_width, delta_height, delta_depth);

		double margin[3];
		margin[0] = (double)width * (1.0f - scale) * 0.5f;
		margin[1] = (double)height * (1.0f - scale) * 0.5f;
		margin[2] = (double)depth * (1.0f - scale) * 0.5f;

		//printf("margin : %f %f %f\n", margin[0], margin[1], margin[2]);
		/*
		for (vector<str_vertex>::size_type i = 0; i < mVtx.size(); i++)
		{
			printf("#%d %lf %lf %lf\n", i, mVtx[i].x, mVtx[i].y, mVtx[i].z);
		}

		for (vector<str_vertex>::size_type i = 0; i < mFace.size(); i++)
		{
			printf("#%d %d %d %d\n", i, mFace[i].idx0, mFace[i].idx1, mFace[i].idx2);
		}
		*/
		/*
		//	put DefValue on grid
		while(!feof(fp))
		{
			fscanf(fp, "%c %lf %lf %lf", &mId, &mVertex_x, &mVertex_y, &mVertex_z);

			if(mId == '#')	// have to drop next strings
			{
				fgets(buf_to_drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				//printf("%c %lf %lf %lf ", mId, mVertex_x, mVertex_y, mVertex_z);

				mX = (int)(abs((mVertex_x - mAA[0])/delta_width*scale) + margin[0]);
				mY = (int)(abs((mVertex_y - mAA[1])/delta_height*scale) + margin[1]);
				mZ = (int)(abs((mVertex_z - mAA[2])/delta_depth*scale) + margin[2]);
				
				//printf("%f ", scale);
				//printf("%lf %lf %lf ", delta_width, delta_height, delta_depth);
				//printf("%lf %lf %lf ", mAA[0], mAA[1], mAA[2]);
				//printf("%lf %lf %lf ", mBB[0], mBB[1], mBB[2]);

				if ((mX > 0 && mX <= mWidth) && (mY > 0 && mY <= mHeight) && (mZ > 0 && mZ <= mDepth))
					arr[(mX)+((mWidth+2)*mY+((mWidth+2)*(mHeight+2)*mZ))] += mDefValue;
				//printf("%d %d %d %d %d %d\n", mX, mY, mZ, mWidth, mHeight, mDepth);
			}

		}
	*/

		if (mVtx.size() < 3) 
		{
			mUTIL.Terminate_system("Mesh model is too small.");
		}

		//for (vector<str_vertex>::size_type i = 0; i < mVtx.size()-2; i++)
		for (vector<str_vertex>::size_type i = 0; i < mFace.size(); i++)
		{
			int idx0 = mFace[i].idx0 - 1;
			int idx1 = mFace[i].idx1 - 1;
			int idx2 = mFace[i].idx2 - 1;

			//printf("%d %d %d\t", idx0, idx1, idx2);

			int x0 = (int)(abs((mVtx[idx0].x - mAA[0]) / delta * scale) + margin[0]);
			int y0 = (int)(abs((mVtx[idx0].y - mAA[1]) / delta * scale) + margin[1]);
			int z0 = (int)(abs((mVtx[idx0].z - mAA[2]) / delta * scale) + margin[2]);

			int x1 = (int)(abs((mVtx[idx1].x - mAA[0]) / delta * scale) + margin[0]);
			int y1 = (int)(abs((mVtx[idx1].y - mAA[1]) / delta * scale) + margin[1]);
			int z1 = (int)(abs((mVtx[idx1].z - mAA[2]) / delta * scale) + margin[2]);

			int x2 = (int)(abs((mVtx[idx2].x - mAA[0]) / delta * scale) + margin[0]);
			int y2 = (int)(abs((mVtx[idx2].y - mAA[1]) / delta * scale) + margin[1]);
			int z2 = (int)(abs((mVtx[idx2].z - mAA[2]) / delta * scale) + margin[2]);


			// add margin
			x0 += margin_x;		x1 += margin_x;		x2 += margin_x;
			y0 += margin_y;		y1 += margin_y;		y2 += margin_y;
			z0 += margin_z;		z1 += margin_z;		z2 += margin_z;


			// 가장 긴 변을 찾는다
			double dist1 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0);
			double dist2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
			double dist3 = (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) + (z0 - z2) * (z0 - z2);

			double max_dist_edge = MAX(MAX(dist1, dist2), dist3);

			if (max_dist_edge <= 1.0f)	
			{
				arr[IX3(x0, y0, z0)] += mDefValue;
				
				//arr[IX3((int)(x0 + delta_x * (double)i), (int)(y0 + delta_y * (double)i), (int)(z0 + delta_z * (double)i))] = mDefValue;
				continue;
			}

			vector3f current, target1, target2, delta;
			
			if (dist1 >= max_dist_edge)
			{
				current.x = x2;		current.y = y2;		current.z = z2;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x1;		target2.y = y1;		target2.z = z1;
			}
			else if (dist2 >= max_dist_edge)
			{
				current.x = x0;		current.y = y0;		current.z = z0;
				target1.x = x1;		target1.y = y1;		target1.z = z1;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			}
			else
			{
				current.x = x1;		current.y = y1;		current.z = z1;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			}
			
			double max_dist_axis = MAX(MAX(abs(target1.x - target2.x), abs(target1.y - target2.y)), abs(target1.z - target2.z));
			max_dist_edge = sqrt(max_dist_edge);
			//delta = (target2 - target1) / max_dist_axis;
			delta.x = (target2.x - target1.x) / max_dist_axis;
			delta.y = (target2.y - target1.y) / max_dist_axis;
			delta.z = (target2.z - target1.z) / max_dist_axis;
			//printf("delta : %lf %lf %lf\n", delta.x, delta.y, delta.z);

			//while ((target1.x <= target2.x) && (target1.y <= target2.y) && (target1.z <= target2.z))
			//while ((abs(target1.x - target2.x) > abs(delta.x)) || (abs(target1.y - target2.y) > abs(delta.y)) || (abs(target1.z - target2.z) > abs(delta.z)))
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);
			//Fill_density_between_two_vertices(target1.x, target1.y, target1.z, target2.x, target2.y, target2.z, arr, 100.0f);
			
			for (int i = 0; i <= (int)max_dist_axis; i++)
			{
				//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, arr, 100.0f);
				//printf("go : %lf %lf %lf -> %lf %lf %lf\n", target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, target2.x, target2.y, target2.z);
				//target1 += delta;
			}
			//printf("(%lf %lf %lf)\n", delta.x, delta.y, delta.z);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);
			
			
			//Fill_density_between_two_vertices(x0, y0, z0, x1, y1, z1, arr, 10.0f);
			//Fill_density_between_two_vertices(x1, y1, z1, x2, y2, z2, arr, 10.0f);
			//Fill_density_between_two_vertices(x2, y2, z2, x0, y0, z0, arr, 10.0f);

			// 가장 긴 변을 찾는다

			// 1. 직선 구하기
			// 가장 긴 변을 찾는다
			// 긴 변을 카운터로 한 줄 그리기

			// 2. 각 점에 대해 나머지 한 줄에 대해서 한 줄 그리기

		}
	}


	void Get_AABB(double* AA, double* BB, int width, int height, int depth, double scale = 0.5f, double mDefValue = 0.1f, int margin_x = 0, int margin_y = 0, int margin_z = 0)
	{
		mWidth = width;
		mHeight = height;
		mDepth = depth;
		//int mCount = 0;

		char mId;
		char buf_to_drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;

		errno_t err;

		FILE *fp;




		//	init arr mAA, mBB
		for (int i = 0; i<3; i++)
		{
			AA[i] = 9999999.0f;
			BB[i] = -9999999.0f;
		}

		for (int i = 0; i < mMeshLoader.Data.verts.size(); ++i)
		{
			auto &v = mMeshLoader.Data.verts[i];

			if (AA[0] > v.x) AA[0] = v.x;
			if (AA[1] > v.y) AA[1] = v.y;
			if (AA[2] > v.z) AA[2] = v.z;

			if (BB[0] < v.x) BB[0] = v.x;
			if (BB[1] < v.y) BB[1] = v.y;
			if (BB[2] < v.z) BB[2] = v.z;
		}
	}

	void Get_AABB_from_obj(double* AA, double* BB, char* full_file_name, int width, int height, int depth, double scale = 0.5f, double mDefValue = 0.1f, int margin_x = 0, int margin_y = 0, int margin_z = 0)
	{
		mWidth = width;
		mHeight = height;
		mDepth = depth;
		//int mCount = 0;

		char mId;
		char buf_to_drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;

		errno_t err;

		FILE *fp;

		


		//	init arr mAA, mBB
		for (int i = 0; i<3; i++)
		{
			AA[i] = 9999999.0f;
			BB[i] = -9999999.0f;
		}

		if ((err = fopen_s(&fp, full_file_name, "r")) != 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}

		// set mAA, mBB
		while (!feof(fp))
		{
			fscanf_s(fp, "%c %lf %lf %lf", &mId, 1, &mVertex_x, &mVertex_y, &mVertex_z);

			if (mId == '#')	// have to drop next strings
			{
				fgets(buf_to_drop_string, 255, fp);
			}
			else if (mId == 'v')
			{
				if (abs(mVertex_x) < 99999999.9f && abs(mVertex_y) < 99999999.9f && abs(mVertex_z) < 99999999.9f)
				{
					// make mAA get the smallest position
					if (AA[0] > mVertex_x) AA[0] = mVertex_x;
					if (AA[1] > mVertex_y) AA[1] = mVertex_y;
					if (AA[2] > mVertex_z) AA[2] = mVertex_z;

					// make mBB get the largest position
					if (BB[0] < mVertex_x) BB[0] = mVertex_x;
					if (BB[1] < mVertex_y) BB[1] = mVertex_y;
					if (BB[2] < mVertex_z) BB[2] = mVertex_z;
				}
			}
		}

		// reload the obj file
		rewind(fp);
		fclose(fp);

	}





	void Put_obj_to_grid_with_AABB(double *arr, double *AA, double *BB, char* full_file_name, int width, int height, int depth, double scale = 0.5f, double mDefValue = 0.1f, int margin_x = 0, int margin_y = 0, int margin_z = 0)
	{
		mWidth = width;
		mHeight = height;
		mDepth = depth;
		//int mCount = 0;

		char mId;
		char buf_to_drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;

		struct str_vertex {
			double x; double y; double z;
		};
		std::vector<str_vertex> mVtx;

		struct str_face_idx {
			int idx0; int idx1; int idx2;
		};
		std::vector<str_face_idx> mFace;
		//int* face_idx[3];

		errno_t err;

		FILE *fp;



		double mAA[3], mBB[3];			// boundary of obj
		//	init arr mAA, mBB
		for (int i = 0; i < 3; i++)
		{
			mAA[i] = AA[i];
			mBB[i] = BB[i];
		}

		if ((err = fopen_s(&fp, full_file_name, "r")) != 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}

		// set mAA, mBB
		while (!feof(fp))
		{
			fscanf_s(fp, "%c %lf %lf %lf", &mId, 1, &mVertex_x, &mVertex_y, &mVertex_z);

			if (mId == '#')	// have to drop next strings
			{
				fgets(buf_to_drop_string, 255, fp);
			}
			else if (mId == 'v')
			{
				if (abs(mVertex_x) < 99999999.9f && abs(mVertex_y) < 99999999.9f && abs(mVertex_z) < 99999999.9f)
				{
					// vector container
					str_vertex temp_vtx = { mVertex_x, mVertex_y, mVertex_z };
					mVtx.push_back(temp_vtx);
				}
			}
			else if (mId == 'f')
			{
				//str_vertex temp_face = { mVertex_x, mVertex_y, mVertex_z };
				//mFace.push_back([);
				mFace.push_back({ (int)mVertex_x, (int)mVertex_y, (int)mVertex_z });
				//mFace.data;
			}
		}

		// reload the obj file
		rewind(fp);
		fclose(fp);




		// set half_AABB_size
		double edge_size_width = (double)(mBB[0] - mAA[0]);
		double edge_size_height = (double)(mBB[1] - mAA[1]);
		double edge_size_depth = (double)(mBB[2] - mAA[2]);

		double max_edge_size = MAX(MAX(edge_size_width, edge_size_height), edge_size_depth);

		int max_edge_res = 0;
		if (max_edge_size <= edge_size_width) max_edge_res = mWidth;
		else if (max_edge_size <= edge_size_height) max_edge_res = mHeight;
		else max_edge_res = mDepth;

		double delta = max_edge_size / max_edge_res;

		//double max_edge = MAX(edge_size_width, edge_size_height)

		// 이게 잘못 되었다! 아닌가?
		//double delta_width = (double)abs(mBB[0]-mAA[0])/mWidth;
		//double delta_height = (double)abs(mBB[1]-mAA[1])/mHeight;
		//double delta_depth = (double)abs(mBB[2]-mAA[2])/mDepth;


		//printf("%lf %lf %lf \t %lf %lf %lf \n", edge_size_width, edge_size_height, edge_size_depth, delta_width, delta_height, delta_depth);

		double margin[3];
		margin[0] = (double)mWidth * (1.0f - scale) * 0.5f;
		margin[1] = (double)mHeight * (1.0f - scale) * 0.5f;
		margin[2] = (double)mDepth * (1.0f - scale) * 0.5f;


		if (mVtx.size() < 3)
		{
			mUTIL.Terminate_system("Mesh model is too small.");
		}

		//for (vector<str_vertex>::size_type i = 0; i < mVtx.size()-2; i++)
		for (vector<str_vertex>::size_type i = 0; i < mFace.size(); i++)
		{
			int idx0 = mFace[i].idx0 - 1;
			int idx1 = mFace[i].idx1 - 1;
			int idx2 = mFace[i].idx2 - 1;

			//printf("%d %d %d\t", idx0, idx1, idx2);

			int x0 = (int)(abs((mVtx[idx0].x - mAA[0]) / delta * scale) + margin[0]);
			int y0 = (int)(abs((mVtx[idx0].y - mAA[1]) / delta * scale) + margin[1]);
			int z0 = (int)(abs((mVtx[idx0].z - mAA[2]) / delta * scale) + margin[2]);

			int x1 = (int)(abs((mVtx[idx1].x - mAA[0]) / delta * scale) + margin[0]);
			int y1 = (int)(abs((mVtx[idx1].y - mAA[1]) / delta * scale) + margin[1]);
			int z1 = (int)(abs((mVtx[idx1].z - mAA[2]) / delta * scale) + margin[2]);

			int x2 = (int)(abs((mVtx[idx2].x - mAA[0]) / delta * scale) + margin[0]);
			int y2 = (int)(abs((mVtx[idx2].y - mAA[1]) / delta * scale) + margin[1]);
			int z2 = (int)(abs((mVtx[idx2].z - mAA[2]) / delta * scale) + margin[2]);

			// add margin
			x0 += margin_x;		x1 += margin_x;		x2 += margin_x;
			y0 += margin_y;		y1 += margin_y;		y2 += margin_y;
			z0 += margin_z;		z1 += margin_z;		z2 += margin_z;


			// 바운더리를 넘어갔나?
			if ((x0 < 1) || (x0 > width) || (x1 < 1) || (x1 > width) || (x2 < 1) || (x2 > width) ||
				(y0 < 1) || (y0 > height) || (y1 < 1) || (y1 > height) || (y2 < 1) || (y2 > height) ||
				(z0 < 1) || (z0 > depth) || (z1 < 1) || (z1 > depth) || (z2 < 1) || (z2 > depth)
				) continue;

			// 가장 긴 변을 찾는다
			double dist1 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0);
			double dist2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
			double dist3 = (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) + (z0 - z2) * (z0 - z2);

			double max_dist_edge = MAX(MAX(dist1, dist2), dist3);

			if (max_dist_edge <= 1.0f)
			{
				arr[IX3(x0, y0, z0)] += mDefValue;

				//arr[IX3((int)(x0 + delta_x * (double)i), (int)(y0 + delta_y * (double)i), (int)(z0 + delta_z * (double)i))] = mDefValue;
				continue;
			}

			vector3f current, target1, target2, delta;

			if (dist1 >= max_dist_edge)
			{
				current.x = x2;		current.y = y2;		current.z = z2;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x1;		target2.y = y1;		target2.z = z1;
			}
			else if (dist2 >= max_dist_edge)
			{
				current.x = x0;		current.y = y0;		current.z = z0;
				target1.x = x1;		target1.y = y1;		target1.z = z1;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			}
			else
			{
				current.x = x1;		current.y = y1;		current.z = z1;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			}

			double max_dist_axis = MAX(MAX(abs(target1.x - target2.x), abs(target1.y - target2.y)), abs(target1.z - target2.z));
			max_dist_edge = sqrt(max_dist_edge);
			//delta = (target2 - target1) / max_dist_axis;
			delta.x = (target2.x - target1.x) / max_dist_axis;
			delta.y = (target2.y - target1.y) / max_dist_axis;
			delta.z = (target2.z - target1.z) / max_dist_axis;
			//printf("delta : %lf %lf %lf\n", delta.x, delta.y, delta.z);

			//while ((target1.x <= target2.x) && (target1.y <= target2.y) && (target1.z <= target2.z))
			//while ((abs(target1.x - target2.x) > abs(delta.x)) || (abs(target1.y - target2.y) > abs(delta.y)) || (abs(target1.z - target2.z) > abs(delta.z)))
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);
			//Fill_density_between_two_vertices(target1.x, target1.y, target1.z, target2.x, target2.y, target2.z, arr, 100.0f);

			for (int i = 0; i <= (int)max_dist_axis; i++)
			{
				//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, arr, 100.0f);
				//printf("go : %lf %lf %lf -> %lf %lf %lf\n", target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, target2.x, target2.y, target2.z);
				//target1 += delta;
			}
			//printf("(%lf %lf %lf)\n", delta.x, delta.y, delta.z);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);


			//Fill_density_between_two_vertices(x0, y0, z0, x1, y1, z1, arr, 10.0f);
			//Fill_density_between_two_vertices(x1, y1, z1, x2, y2, z2, arr, 10.0f);
			//Fill_density_between_two_vertices(x2, y2, z2, x0, y0, z0, arr, 10.0f);

			// 가장 긴 변을 찾는다

			// 1. 직선 구하기
			// 가장 긴 변을 찾는다
			// 긴 변을 카운터로 한 줄 그리기

			// 2. 각 점에 대해 나머지 한 줄에 대해서 한 줄 그리기

		}
	}


	void Put_to_grid_with_AABB(double *arr, double *nx, double *ny, double *nz, double *AA, double *BB, int width, int height, int depth, double scale = 0.5f, double mDefValue = 0.1f, int margin_x = 0, int margin_y = 0, int margin_z = 0)
	{
		mWidth = width;
		mHeight = height;
		mDepth = depth;
		//int mCount = 0;

		char mId;
		char buf_to_drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;

		struct str_vertex {
			double x; double y; double z;
		};
		std::vector<str_vertex> mVtx;

		struct str_face_idx {
			int idx0; int idx1; int idx2;
		};
		std::vector<str_face_idx> mFace;
		//int* face_idx[3];

		errno_t err;

		FILE *fp;



		double mAA[3], mBB[3];			// boundary of obj
		//	init arr mAA, mBB
		for (int i = 0; i < 3; i++) {
			mAA[i] = AA[i];
			mBB[i] = BB[i];
		}

		//if ((err = fopen_s(&fp, full_file_name, "r")) != 0)
		//{
		//	printf("Cannot open %s\n", full_file_name);
		//	exit(0);
		//}

		//// set mAA, mBB
		//while (!feof(fp))
		//{
		//	fscanf_s(fp, "%c %lf %lf %lf", &mId, 1, &mVertex_x, &mVertex_y, &mVertex_z);

		//	if (mId == '#')	// have to drop next strings
		//	{
		//		fgets(buf_to_drop_string, 255, fp);
		//	}
		//	else if (mId == 'v')
		//	{
		//		if (abs(mVertex_x) < 99999999.9f && abs(mVertex_y) < 99999999.9f && abs(mVertex_z) < 99999999.9f)
		//		{
		//			// vector container
		//			str_vertex temp_vtx = { mVertex_x, mVertex_y, mVertex_z };
		//			mVtx.push_back(temp_vtx);
		//		}
		//	}
		//	else if (mId == 'f')
		//	{
		//		//str_vertex temp_face = { mVertex_x, mVertex_y, mVertex_z };
		//		//mFace.push_back([);
		//		mFace.push_back({ (int)mVertex_x, (int)mVertex_y, (int)mVertex_z });
		//		//mFace.data;
		//	}
		//}

		//// reload the obj file
		//rewind(fp);
		//fclose(fp);


		for (int i = 0; i < mMeshLoader.Data.verts.size(); ++i) {
			auto &v = mMeshLoader.Data.verts[i];

			mVtx.push_back({ v.x, v.y, v.z });
		}

		for (int i = 0; i < mMeshLoader.Data.inds.size(); i += 3) {
			mFace.push_back({
				mMeshLoader.Data.inds[i + 0],
				mMeshLoader.Data.inds[i + 1],
				mMeshLoader.Data.inds[i + 2],
			});
		}

		// set half_AABB_size
		double edge_size_width = (double)(mBB[0] - mAA[0]);
		double edge_size_height = (double)(mBB[1] - mAA[1]);
		double edge_size_depth = (double)(mBB[2] - mAA[2]);

		double max_edge_size = MAX(MAX(edge_size_width, edge_size_height), edge_size_depth);

		int max_edge_res = 0;
		if (max_edge_size <= edge_size_width) max_edge_res = mWidth;
		else if (max_edge_size <= edge_size_height) max_edge_res = mHeight;
		else max_edge_res = mDepth;

		double delta = max_edge_size / max_edge_res;

		//double max_edge = MAX(edge_size_width, edge_size_height)

		// 이게 잘못 되었다! 아닌가?
		//double delta_width = (double)abs(mBB[0]-mAA[0])/mWidth;
		//double delta_height = (double)abs(mBB[1]-mAA[1])/mHeight;
		//double delta_depth = (double)abs(mBB[2]-mAA[2])/mDepth;


		//printf("%lf %lf %lf \t %lf %lf %lf \n", edge_size_width, edge_size_height, edge_size_depth, delta_width, delta_height, delta_depth);

		double margin[3];
		margin[0] = (double)mWidth * (1.0f - scale) * 0.5f;
		margin[1] = (double)mHeight * (1.0f - scale) * 0.5f;
		margin[2] = (double)mDepth * (1.0f - scale) * 0.5f;


		if (mVtx.size() < 3) {
			mUTIL.Terminate_system("Mesh model is too small.");
		}

		//for (vector<str_vertex>::size_type i = 0; i < mVtx.size()-2; i++)
		for (vector<str_vertex>::size_type i = 0; i < mFace.size(); i++) {
			int idx0 = mFace[i].idx0;
			int idx1 = mFace[i].idx1;
			int idx2 = mFace[i].idx2;

			//printf("%d %d %d\t", idx0, idx1, idx2);

			int x0 = (int)(abs((mVtx[idx0].x - mAA[0]) / delta * scale) + margin[0]);
			int y0 = (int)(abs((mVtx[idx0].y - mAA[1]) / delta * scale) + margin[1]);
			int z0 = (int)(abs((mVtx[idx0].z - mAA[2]) / delta * scale) + margin[2]);

			int x1 = (int)(abs((mVtx[idx1].x - mAA[0]) / delta * scale) + margin[0]);
			int y1 = (int)(abs((mVtx[idx1].y - mAA[1]) / delta * scale) + margin[1]);
			int z1 = (int)(abs((mVtx[idx1].z - mAA[2]) / delta * scale) + margin[2]);

			int x2 = (int)(abs((mVtx[idx2].x - mAA[0]) / delta * scale) + margin[0]);
			int y2 = (int)(abs((mVtx[idx2].y - mAA[1]) / delta * scale) + margin[1]);
			int z2 = (int)(abs((mVtx[idx2].z - mAA[2]) / delta * scale) + margin[2]);

			// add margin
			x0 += margin_x;		x1 += margin_x;		x2 += margin_x;
			y0 += margin_y;		y1 += margin_y;		y2 += margin_y;
			z0 += margin_z;		z1 += margin_z;		z2 += margin_z;


			// 바운더리를 넘어갔나?
			if ((x0 < 1) || (x0 > width) || (x1 < 1) || (x1 > width) || (x2 < 1) || (x2 > width) ||
				(y0 < 1) || (y0 > height) || (y1 < 1) || (y1 > height) || (y2 < 1) || (y2 > height) ||
				(z0 < 1) || (z0 > depth) || (z1 < 1) || (z1 > depth) || (z2 < 1) || (z2 > depth)
				) continue;

			// 가장 긴 변을 찾는다
			double dist1 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) + (z1 - z0) * (z1 - z0);
			double dist2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
			double dist3 = (x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) + (z0 - z2) * (z0 - z2);

			double max_dist_edge = MAX(MAX(dist1, dist2), dist3);

			glm::vec3 v0(mVtx[idx0].x, mVtx[idx0].y, mVtx[idx0].z);
			glm::vec3 v1(mVtx[idx1].x, mVtx[idx1].y, mVtx[idx1].z);
			glm::vec3 v2(mVtx[idx2].x, mVtx[idx2].y, mVtx[idx2].z);

			glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

			if (max_dist_edge <= 1.0f) {
				arr[IX3(x0, y0, z0)] += mDefValue;

				nx[IX3(x0, y0, z0)] += normal.x;
				ny[IX3(x0, y0, z0)] += normal.y;
				nz[IX3(x0, y0, z0)] += normal.z;

				//arr[IX3((int)(x0 + delta_x * (double)i), (int)(y0 + delta_y * (double)i), (int)(z0 + delta_z * (double)i))] = mDefValue;
				continue;
			}

			vector3f current, target1, target2, delta;

			if (dist1 >= max_dist_edge) {
				current.x = x2;		current.y = y2;		current.z = z2;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x1;		target2.y = y1;		target2.z = z1;
			} else if (dist2 >= max_dist_edge) {
				current.x = x0;		current.y = y0;		current.z = z0;
				target1.x = x1;		target1.y = y1;		target1.z = z1;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			} else {
				current.x = x1;		current.y = y1;		current.z = z1;
				target1.x = x0;		target1.y = y0;		target1.z = z0;
				target2.x = x2;		target2.y = y2;		target2.z = z2;
			}

			double max_dist_axis = MAX(MAX(abs(target1.x - target2.x), abs(target1.y - target2.y)), abs(target1.z - target2.z));
			max_dist_edge = sqrt(max_dist_edge);
			//delta = (target2 - target1) / max_dist_axis;
			delta.x = (target2.x - target1.x) / max_dist_axis;
			delta.y = (target2.y - target1.y) / max_dist_axis;
			delta.z = (target2.z - target1.z) / max_dist_axis;
			//printf("delta : %lf %lf %lf\n", delta.x, delta.y, delta.z);

			//while ((target1.x <= target2.x) && (target1.y <= target2.y) && (target1.z <= target2.z))
			//while ((abs(target1.x - target2.x) > abs(delta.x)) || (abs(target1.y - target2.y) > abs(delta.y)) || (abs(target1.z - target2.z) > abs(delta.z)))
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);
			//Fill_density_between_two_vertices(target1.x, target1.y, target1.z, target2.x, target2.y, target2.z, arr, 100.0f);

			for (int i = 0; i <= (int)max_dist_axis; i++) {
				//Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x, target1.y, target1.z, arr, 100.0f);
				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, arr, mDefValue);

				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, nx, normal.x);
				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, ny, normal.y);
				Fill_density_between_two_vertices(current.x, current.y, current.z, target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, nz, normal.z);
				//printf("go : %lf %lf %lf -> %lf %lf %lf\n", target1.x + delta.x * (double)i, target1.y + delta.y * (double)i, target1.z + delta.z * (double)i, target2.x, target2.y, target2.z);
				//target1 += delta;
			}
			//printf("(%lf %lf %lf)\n", delta.x, delta.y, delta.z);
			//Fill_density_between_two_vertices(current.x, current.y, current.z, target2.x, target2.y, target2.z, arr, 100.0f);


			//Fill_density_between_two_vertices(x0, y0, z0, x1, y1, z1, arr, 10.0f);
			//Fill_density_between_two_vertices(x1, y1, z1, x2, y2, z2, arr, 10.0f);
			//Fill_density_between_two_vertices(x2, y2, z2, x0, y0, z0, arr, 10.0f);

			// 가장 긴 변을 찾는다

			// 1. 직선 구하기
			// 가장 긴 변을 찾는다
			// 긴 변을 카운터로 한 줄 그리기

			// 2. 각 점에 대해 나머지 한 줄에 대해서 한 줄 그리기

			//printf("%d %d %d = %d\n", width, height, depth, (width + 2)*(height + 2)*(depth + 2));
		}

#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++) {
			for (int j = 1; j <= mHeight; j++) {
				for (int k = 1; k <= mDepth; k++) {
					if (arr[IX3(i, j, k)] > 0.0f) {
						float x = nx[IX3(i, j, k)];
						float y = ny[IX3(i, j, k)];
						float z = nz[IX3(i, j, k)];

						glm::fvec3 n(x, y, z);
						n = glm::normalize(n);

						nx[IX3(i, j, k)] = n.x;
						ny[IX3(i, j, k)] = n.y;
						nz[IX3(i, j, k)] = n.z;
					}
				}
			}
		}
	}


	/*
	void Put_obj_to_grid(double *arr, char* full_file_name, int mWidth, int mHeight, int mDepth, int AAx, int AAy, int AAz, int BBx, int BBy, int BBz, double mDefValue = 0.1)
	{
		int mCount = 0;
	
		char mId;
		char drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		int mX, mY, mZ;

		double mAA[3], mBB[3];			// boundary of obj
		//	double mDefValue = 0.5;		// default adding value for vertex on grid

		errno_t err;

		FILE *fp;



		//	init arr mAA, mBB
		for(int i=0; i<3; i++)
		{
			mAA[i] = 999999.0f;
			mBB[i] = -999999.0f;
		}


	
		if((err = fopen_s(&fp, full_file_name, "r"))!= 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}



		// set mAA, mBB
		while(!feof(fp))
		{
			fscanf_s(fp, "%c %f %f %f", &mId, &mVertex_x, &mVertex_y, &mVertex_z);

			if(mId == '#')	// have to drop next strings
			{
				fgets(drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				// make mAA get the smallest position
				if(mAA[0] > mVertex_x) mAA[0] = mVertex_x;
				if(mAA[1] > mVertex_y) mAA[1] = mVertex_y;
				if(mAA[2] > mVertex_z) mAA[2] = mVertex_z;

				// make mBB get the largest position
				if(mBB[0] < mVertex_x) mBB[0] = mVertex_x;
				if(mBB[1] < mVertex_y) mBB[1] = mVertex_y;
				if(mBB[2] < mVertex_z) mBB[2] = mVertex_z;
			}
		}

		// reload the obj file
		rewind(fp);		 

		// set half_AABB_size

		double half_AABB_size_width = (double)(mBB[0]-mAA[0])/2.0f;
		double half_AABB_size_height = (double)(mBB[1]-mAA[1])/2.0f;
		double half_AABB_size_depth = (double)(mBB[2]-mAA[2])/2.0f;

	
		//	half_AABB_size_width = half_AABB_size_height = half_AABB_size_depth = 2.0f;
		//	half_AABB_size_height = half_AABB_size_depth = half_AABB_size_width;
	

		//	put DefValue on grid
		while(!feof(fp))
		{
			fscanf_s(fp, "%c %f %f %f", &mId, &mVertex_x, &mVertex_y, &mVertex_z);

			if(mId == '#')	// have to drop next strings
			{
				fgets(drop_string, 255, fp);	
			}
			else if(mId == 'v')
			{
				int scale = 1;

				mX = Position_on_grid(mVertex_x, mAA[0], half_AABB_size_width, mWidth/scale);
				mY = Position_on_grid(mVertex_y, mAA[1], half_AABB_size_height, mHeight/scale);
				mZ = Position_on_grid(mVertex_z, mAA[2], half_AABB_size_depth, mDepth/scale);


				//	printf("#%d %d(%f) %d(%f) %d(%f)\n", mCount, mX, mVertex_x, mY, mVertex_y, mZ, mVertex_z);
				arr[(mX)+((mWidth+2)*mY+((mWidth+2)*(mHeight+2)*mZ))] += mDefValue;
				
				mCount++;
			}

		}
	
		printf("\t%d data were read.\n", mCount);

		fclose(fp);
	}
	*/


	void Fill_density(double *dens, double *temperature, int mWidth, int mHeight, int mDepth)
	{
		int toggle = 0;

		for(int i=1; i<=mWidth; i++)
			for(int j=1; j<=mHeight; j++)
				for(int k=1; k<=mDepth; k++)
				{
					toggle = 0;
					if(dens[IX3(i,j,k)] > 0) toggle = !toggle;

					if(toggle) 
					{
						dens[IX3(i,j,k)] = 0.8f;
						temperature[IX3(i,j,k)] = 0.5f;
					}
					//else dens[IX3(i,j,k)] = 0;
				}
	}





	void Fill_density(double *dens, int mWidth, int mHeight, int mDepth)
	{
		int toggle = 0;

		for(int i=1; i<=mWidth; i++)
			for(int j=1; j<=mHeight; j++)
				for(int k=1; k<=mDepth; k++)
				{
					toggle = 0;
					if(dens[IX3(i,j,k)] > 0) toggle = !toggle;

					if(toggle) 
					{
						dens[IX3(i,j,k)] = 1.0f;
					}
					else dens[IX3(i,j,k)] = 0;
				}
	}





	void Put_phi_from_obj_to_grid(double *arr, char* full_file_name, int mWidth, int mHeight, int mDepth)
	{
		char mId;
		char drop_string[256];
		double mVertex_x, mVertex_y, mVertex_z;
		double mX, mY, mZ;

		double mAA[3], mBB[3];			// boundary of obj

		errno_t err;

		FILE *fp;



		//	init arr mAA, mBB
		for(int i=0; i<3; i++)
		{
			mAA[i] = 999999.0f;
			mBB[i] = -999999.0f;
		}


	
		if((err = fopen_s(&fp, full_file_name, "r"))!= 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}



		// set mAA, mBB
		while(!feof(fp))
		{
			fscanf_s(fp, "%c %lf %lf %lf", &mId, &mVertex_x, &mVertex_y, &mVertex_z);

			if(mId == '#')	// have to drop next strings
			{
				fgets(drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				// make mAA get the smallest position
				if(mAA[0] > mVertex_x) mAA[0] = mVertex_x;
				if(mAA[1] > mVertex_y) mAA[1] = mVertex_y;
				if(mAA[2] > mVertex_z) mAA[2] = mVertex_z;

				// make mBB get the largest position
				if(mBB[0] < mVertex_x) mBB[0] = mVertex_x;
				if(mBB[1] < mVertex_y) mBB[1] = mVertex_y;
				if(mBB[2] < mVertex_z) mBB[2] = mVertex_z;
			}
		}

		// reload the obj file
		rewind(fp);

		// set half_AABB_size
		double half_AABB_size_width = (double)(mBB[0]-mAA[0])/2.0f;
		double half_AABB_size_height = (double)(mBB[1]-mAA[1])/2.0f;
		double half_AABB_size_depth = (double)(mBB[2]-mAA[2])/2.0f;

		int min_pos[3] = {-1, -1, -1};
		double min_val = 999999;

		// init arr
		for(int i=0; i<=mWidth+1; i++)
			for(int j=0; j<=mHeight+1; j++)
				for(int k=0; k<=mDepth+1; k++)
					arr[IX3(i,j,k)] = 999999;



		//	put DefValue on grid
		while(!feof(fp))
		{
			fscanf_s(fp, "%c %lf %lf %lf", &mId, &mVertex_x, &mVertex_y, &mVertex_z);

			if(mId == '#')	// have to drop next strings
			{
				fgets(drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				int scale = 2;

				mX = Phi_on_grid(mVertex_x, mAA[0], half_AABB_size_width, mWidth/scale);
				mY = Phi_on_grid(mVertex_y, mAA[1], half_AABB_size_height, mHeight/scale);
				mZ = Phi_on_grid(mVertex_z, mAA[2], half_AABB_size_depth, mDepth/scale);

				mX += mWidth/(scale*2);
				//mY += height/(scale*2);
				mY += mHeight/(scale*6);
				mZ += mDepth/(scale*2);

				double diff_x = mX - (int)mX;
				double diff_y = mY - (int)mY;
				double diff_z = mZ - (int)mZ;

				double s = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);

				if(s < min_val)
				{
					min_val = s;
					min_pos[0] = (int)mX, min_pos[1] = (int)mY, min_pos[2] = (int)mZ;
				}
				
				if(arr[IX3((int)mX,(int)mY,(int)mZ)] > s) 
					arr[IX3((int)mX,(int)mY,(int)mZ)] = s;
				
				//	printf("#%d %d(%f) %d(%f) %d(%f)\n", mCount, mX, mVertex_x, mY, mVertex_y, mZ, mVertex_z);
				//arr[(mX)+((width+2)*mY+((width+2)*(height+2)*mZ))] += mDefValue;
			}
		}

		double h = 1.0f;						// length of edge
		
		printf("\t\tMake a unsiged distance field\n");
		Make_unsigned_distance_field(arr, mWidth, mHeight, mDepth, min_pos[0], min_pos[1], min_pos[2], h);
		
		printf("\t\tAdd sign to a unsigned distance field\n");
		Add_sign_to_unsigned_distance_field(arr, mWidth, mHeight, mDepth, h);
		

		


/*
			// check neighbor
			if(mPos[0] > 1)
			{
				if(flag[IX3(mPos[0]-1, mPos[1], mPos[2])] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0]-1, mPos[1], mPos[2])]) > h ) 
					{
						flag[IX3(mPos[0]-1, mPos[1], mPos[2])] = false;
						stack.push(mPos[0]-1, mPos[1], mPos[2]);
					}
					else 
					{
						//flag[IX3(mPos[0]-1, mPos[1], mPos[2])] = true;
						mDist_min = MIN(arr[IX3(mPos[0]-1, mPos[1], mPos[2])]+h, mDist_min);
					}
				}
				else stack.push(mPos[0]-1, mPos[1], mPos[2]);
			}
			if(mPos[0] < width)
			{
				if(flag[IX3(mPos[0]+1, mPos[1], mPos[2])] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0]+1, mPos[1], mPos[2])]) > h )
					{
						flag[IX3(mPos[0]+1, mPos[1], mPos[2])] = false;
						stack.push(mPos[0]+1, mPos[1], mPos[2]);
					}
					else 
					{
						//flag[IX3(mPos[0]+1, mPos[1], mPos[2])] = true;
						mDist_min = MIN(arr[IX3(mPos[0]+1, mPos[1], mPos[2])]+h, mDist_min);
					}
				}
				else stack.push(mPos[0]+1, mPos[1], mPos[2]);
			}
			if(mPos[1] > 1)
			{
				if(flag[IX3(mPos[0], mPos[1]-1, mPos[2])] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0], mPos[1]-1, mPos[2])]) > h )
					{
						flag[IX3(mPos[0], mPos[1]-1, mPos[2])] = false;
						stack.push(mPos[0], mPos[1]-1, mPos[2]);
					}
					else 
					{
						//flag[IX3(mPos[0], mPos[1]-1, mPos[2])] = true;
						mDist_min = MIN(arr[IX3(mPos[0], mPos[1]-1, mPos[2])]+h, mDist_min);
					}
				}
				else stack.push(mPos[0], mPos[1]-1, mPos[2]);
			}
			if(mPos[1] < height)
			{
				if(flag[IX3(mPos[0], mPos[1]+1, mPos[2])] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0], mPos[1]+1, mPos[2])]) > h ) 
					{
						flag[IX3(mPos[0], mPos[1]+1, mPos[2])] = false;
						stack.push(mPos[0], mPos[1]+1, mPos[2]);
					}
					else 
					{
						//flag[IX3(mPos[0], mPos[1]+1, mPos[2])] = true;
						mDist_min = MIN(arr[IX3(mPos[0], mPos[1]+1, mPos[2])]+h, mDist_min);
					}
				}
				else stack.push(mPos[0], mPos[1]+1, mPos[2]);
			}			
			if(mPos[2] > 1)
			{
				if(flag[IX3(mPos[0], mPos[1], mPos[2]-1)] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0], mPos[1], mPos[2]-1)]) > h ) 
					{
						flag[IX3(mPos[0], mPos[1], mPos[2]-1)] = false;
						stack.push(mPos[0], mPos[1], mPos[2]-1);
					}
					else 
					{
						//flag[IX3(mPos[0], mPos[1], mPos[2]-1)] = true;
						mDist_min = MIN(arr[IX3(mPos[0], mPos[1], mPos[2]-1)]+h, mDist_min);
					}
				}
				else stack.push(mPos[0], mPos[1], mPos[2]-1);
			}
			if(mPos[2] < depth)
			{
				if(flag[IX3(mPos[0], mPos[1], mPos[2]+1)] == true)
				{
					if( abs(mDist-arr[IX3(mPos[0], mPos[1], mPos[2]+1)]) > h ) 
					{
						flag[IX3(mPos[0], mPos[1], mPos[2]+1)] = false;
						stack.push(mPos[0], mPos[1], mPos[2]+1);
					}
					else 
					{
						//flag[IX3(mPos[0], mPos[1], mPos[2]+1)] = true;
						mDist_min = MIN(arr[IX3(mPos[0], mPos[1], mPos[2]+1)]+h, mDist_min);
					}
				}
				else stack.push(mPos[0], mPos[1], mPos[2]+1);
			}		
		

			arr[IX3(mPos[0], mPos[1], mPos[2])] = mDist_min;
			flag[IX3(mPos[0], mPos[1], mPos[2])] = true;
		}
*/

		/*
		double mDist = 999999;
		double mSmallerDist = 999999;

		while(stack.getSize() > 0)
		{
			mPos[2] = stack.pop();	mPos[1] = stack.pop();	mPos[0] = stack.pop();	
			mDist = arr[IX3(mPos[0], mPos[1], mPos[2])];
			printf("POP : %d %d %d(%f), %d left.\n", mPos[0], mPos[1], mPos[2], mDist, stack.getSize());
			
			// check neighbor
			if(mPos[0] > 1)
			{
				if(arr[IX3(mPos[0]-1, mPos[1], mPos[2])] > mDist+h) 
				{
					arr[IX3(mPos[0]-1, mPos[1], mPos[2])] = mDist+h; 
					stack.push(mPos[0]-1, mPos[1], mPos[2]);
				}
				else if(arr[IX3(mPos[0]-1, mPos[1], mPos[2])] < mDist+h) 
				{
					mSmallerDist = arr[IX3(mPos[0]-1, mPos[1], mPos[2])]-h;
					stack.push(mPos[0]-1, mPos[1], mPos[2]);
					//stack.push(mPos[0], mPos[1], mPos[2]);
				}
			}
			if(mPos[0] < width)
			{
				if(arr[IX3(mPos[0]+1, mPos[1], mPos[2])] > mDist+h)
				{
					arr[IX3(mPos[0]+1, mPos[1], mPos[2])] = mDist+h; 
					stack.push(mPos[0]+1, mPos[1], mPos[2]);
				}
				else if(arr[IX3(mPos[0]+1, mPos[1], mPos[2])] < mDist+h)
				{
					mSmallerDist = arr[IX3(mPos[0]+1, mPos[1], mPos[2])]-h;
					stack.push(mPos[0]+1, mPos[1], mPos[2]);
					//stack.push(mPos[0], mPos[1], mPos[2]);
				}
			}
			if(mPos[1] > 1)
			{
				if(arr[IX3(mPos[0], mPos[1]-1, mPos[2])] > mDist+h) 
				{
					arr[IX3(mPos[0], mPos[1]-1, mPos[2])] = mDist+h; 
					stack.push(mPos[0], mPos[1]-1, mPos[2]);
				}
				else if(arr[IX3(mPos[0], mPos[1]-1, mPos[2])] < mDist+h)
				{
					mSmallerDist = arr[IX3(mPos[0], mPos[1]-1, mPos[2])]-h;
					//stack.push(mPos[0], mPos[1], mPos[2]);
					stack.push(mPos[0], mPos[1]-1, mPos[2]);
				}
			}
			if(mPos[1] < height)
			{
				if(arr[IX3(mPos[0], mPos[1]+1, mPos[2])] > mDist+h)
				{
					arr[IX3(mPos[0], mPos[1]+1, mPos[2])] = mDist+h; 
					stack.push(mPos[0], mPos[1]+1, mPos[2]);
				}
				else if(arr[IX3(mPos[0], mPos[1]+1, mPos[2])] < mDist+h)
				{
					mSmallerDist = arr[IX3(mPos[0], mPos[1]+1, mPos[2])]-h;
					//stack.push(mPos[0], mPos[1], mPos[2]);
					stack.push(mPos[0], mPos[1]+1, mPos[2]);
				}
			}
			if(mPos[2] > 1)
			{
				if(arr[IX3(mPos[0], mPos[1], mPos[2])-1] > mDist+h) 
				{
					arr[IX3(mPos[0], mPos[1], mPos[2])-1] = mDist+h; 
					stack.push(mPos[0], mPos[1], mPos[2]-1);
				}
				else if(arr[IX3(mPos[0], mPos[1], mPos[2])-1] < mDist+h)
				{
					mSmallerDist = arr[IX3(mPos[0], mPos[1], mPos[2]-1)]-h;
					//stack.push(mPos[0], mPos[1], mPos[2]);
					stack.push(mPos[0], mPos[1], mPos[2]-1);
				}
			}
			if(mPos[2] < depth)
			{
				if(arr[IX3(mPos[0], mPos[1], mPos[2]+1)] > mDist+h)
				{
					arr[IX3(mPos[0], mPos[1], mPos[2]+1)] = mDist+h; 
					stack.push(mPos[0], mPos[1], mPos[2]+1);
				}
				else if(arr[IX3(mPos[0], mPos[1], mPos[2]+1)] < mDist+h)
				{
					mSmallerDist = arr[IX3(mPos[0], mPos[1], mPos[2]+1)]-h;
					//stack.push(mPos[0], mPos[1], mPos[2]);
					stack.push(mPos[0], mPos[1], mPos[2]+1);
				}
			}
			arr[IX3(mPos[0], mPos[1], mPos[2])] = mSmallerDist;
		}
		*/
		/*
		for(int i=1; i<=width; i++)
			for(int j=1; j<=height; j++)
			{
				for(int k=1; k<=depth; k++)
				{
					temp_min = 999999;

					if(arr[IX3(i,j,k)] < temp_min) 
					{
						temp_min = arr[IX3(i,j,k)];
					}
					else
					{
						arr[IX3(i,j,k)] = temp_min;
						temp_min += h;
					}
				}
				for(int k=depth; k>=1; k--)
				{
					temp_min = 999999;
		
					if(arr[IX3(i,j,k)] < temp_min) 
					{
						temp_min = arr[IX3(i,j,k)];
					}
					else
					{
						arr[IX3(i,j,k)] = temp_min;
						temp_min += h;
					}
				}
			}
		*/


		fclose(fp);
	}





	void Make_unsigned_distance_field(double *arr, int mWidth, int mHeight, int mDepth, int start_x, int start_y, int start_z, double h)
	{
		bool *flag = (bool*)malloc(sizeof(bool)*(mWidth+2)*(mHeight+2)*(mDepth+2));		

		// init arr
		for(int i=0; i<=mWidth+1; i++)
			for(int j=0; j<=mHeight+1; j++)
				for(int k=0; k<=mDepth+1; k++)
				{
					if(arr[IX3(i,j,k)] == 999999) flag[IX3(i,j,k)] = false;
					else flag[IX3(i,j,k)] = true;
				}

		GRIDY_LIST_STACK3<int> stack;
		
		flag[IX3(start_x, start_y, start_z)] = true;
		stack.push(start_x-1, start_y, start_z);
		stack.push(start_x+1, start_y, start_z);
		stack.push(start_x, start_y-1, start_z);
		stack.push(start_x, start_y+1, start_z);
		stack.push(start_x, start_y, start_z-1);
		stack.push(start_x, start_y, start_z+1);
		
		int mPos[3];						// temp arr for stack


		double mDist = 999999;

		while(stack.getSize() > 0)
		{
			mPos[2] = stack.pop();	mPos[1] = stack.pop();	mPos[0] = stack.pop();	
			mDist = arr[IX3(mPos[0], mPos[1], mPos[2])];
			double mDist_min = 999999;
			// printf("POP : %d %d %d(%f), %d left.\n", mPos[0], mPos[1], mPos[2], mDist, stack.getSize());

			// check neighbor
			if(mPos[0] >= 1)
			{
				if(flag[IX3(mPos[0]-1, mPos[1], mPos[2])] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0]-1, mPos[1], mPos[2])]+h, mDist_min);
				}
				else stack.push(mPos[0]-1, mPos[1], mPos[2]);
			}
			if(mPos[0] <= mWidth)
			{
				if(flag[IX3(mPos[0]+1, mPos[1], mPos[2])] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0]+1, mPos[1], mPos[2])]+h, mDist_min);
				}
				else stack.push(mPos[0]+1, mPos[1], mPos[2]);
			}
			if(mPos[1] >= 1)
			{
				if(flag[IX3(mPos[0], mPos[1]-1, mPos[2])] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0], mPos[1]-1, mPos[2])]+h, mDist_min);
				}
				else stack.push(mPos[0], mPos[1]-1, mPos[2]);
			}
			if(mPos[1] <= mHeight)
			{
				if(flag[IX3(mPos[0], mPos[1]+1, mPos[2])] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0], mPos[1]+1, mPos[2])]+h, mDist_min);
				}
				else stack.push(mPos[0], mPos[1]+1, mPos[2]);
			}			
			if(mPos[2] >= 1)
			{
				if(flag[IX3(mPos[0], mPos[1], mPos[2]-1)] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0], mPos[1], mPos[2]-1)]+h, mDist_min);
				}
				else stack.push(mPos[0], mPos[1], mPos[2]-1);
			}
			if(mPos[2] <= mDepth)
			{
				if(flag[IX3(mPos[0], mPos[1], mPos[2]+1)] == true)
				{
					mDist_min = MIN(arr[IX3(mPos[0], mPos[1], mPos[2]+1)]+h, mDist_min);
				}
				else stack.push(mPos[0], mPos[1], mPos[2]+1);
			}		
		

			arr[IX3(mPos[0], mPos[1], mPos[2])] = mDist_min;
			flag[IX3(mPos[0], mPos[1], mPos[2])] = true;
		}


		free(flag);
	}



	void Add_sign_to_unsigned_distance_field(double *arr, int mWidth, int mHeight, int mDepth, double h)
	{
		int sign = -1;
		double last_dist = 0.0f;
		//double max_dist = pow(h, 1.0/3.0);

		for(int i=0; i<=mWidth+1; i++)
			for(int j=0; j<=mHeight+1; j++)
			{
				sign = -1;
				for(int k=0; k<=mDepth+1; k++)
				{
					
					if( abs(arr[IX3(i,j,k)]) < 1.73)
					{
						if(abs(abs(last_dist)-abs(arr[IX3(i,j,k)])) != h)
						{
							sign *= -1;		// boundary found
						}
					}
					arr[IX3(i,j,k)] *= sign;
					last_dist = arr[IX3(i,j,k)];
				}
			}
	}




	void Make_SDF_from_obj(double *phi, char* full_file_name, int mWidth, int mHeight, int mDepth)
	{
		char mId;
		char mBuf_to_drop_string[256];
		double mVertex_position[3];
		
		double mBounding_box_AA[3], mBounding_box_BB[3];			// boundary of obj
		double mEdge_length[3];
		double mDelta_grid[3];


		FILE *fp;
		errno_t err;



		//	init arr mAA, mBB
		for(int i=0; i<3; i++)
		{
			mBounding_box_AA[i] = 9999999.0f;
			mBounding_box_BB[i] = -9999999.0f;
		}
			
		if((err = fopen_s(&fp, full_file_name, "r"))!= 0)
		{
			printf("Cannot open %s\n", full_file_name);
			exit(0);
		}
		
		//	set bounding box
		while(!feof(fp))
		{
			fscanf(fp, "%c %f %f %f", &mId, &mVertex_position[0], &mVertex_position[1], &mVertex_position[2]);
			
			if(mId == '#')	// have to drop next strings
			{
				fgets(mBuf_to_drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				for(int i=0; i<3; i++)
				{
					if(mBounding_box_AA[i] > mVertex_position[i]) mBounding_box_AA[i] = mVertex_position[i];	// make mAA get the smallest value
					if(mBounding_box_BB[i] < mVertex_position[i]) mBounding_box_BB[i] = mVertex_position[i];	// make mBB get the largest value
				}
			}
		}



		//	reload the obj file
		rewind(fp);

		//	get the length of each edge
		for(int i=0; i<3; i++)
			mEdge_length[i] = (double)abs(mBounding_box_BB[i]-mBounding_box_AA[i]);

		mDelta_grid[0] = (double)mEdge_length[0]/mWidth;
		mDelta_grid[1] = (double)mEdge_length[1]/mHeight;
		mDelta_grid[2] = (double)mEdge_length[2]/mDepth;
	

		//	put DefValue on grid
		while(!feof(fp))
		{
			fscanf(fp, "%c %f %f %f", &mId, &mVertex_position[0], &mVertex_position[1], &mVertex_position[2]);

			if(mId == '#')	// have to drop next strings
			{
				fgets(mBuf_to_drop_string, 255, fp);
			}
			else if(mId == 'v')
			{
				double temp_distance_vertex[3], temp_phi, temp_phi0;

				for(int i=0; i<3; i++)
					temp_distance_vertex[i] = abs((mVertex_position[i] - mBounding_box_AA[i])/mDelta_grid[i]);

				temp_phi = sqrt(temp_distance_vertex[0]*temp_distance_vertex[0] + temp_distance_vertex[1]*temp_distance_vertex[1] + temp_distance_vertex[2]*temp_distance_vertex[2]);
				temp_phi0 = phi[IX3((int)temp_distance_vertex[0], (int)temp_distance_vertex[1], (int)temp_distance_vertex[2])];

				if(temp_phi < temp_phi0) phi[IX3((int)temp_distance_vertex[0], (int)temp_distance_vertex[1], (int)temp_distance_vertex[2])] = temp_phi;

				//double temp_phi0 = phi[((int)temp_phi[i])+((width+2)*temp_phi[1]+((width+2)*(height+2)*temp_phi[2]))];
			}
		}
	

		fclose(fp);
	}







};








//	두 점 사이의 직선에 해당하는 셀에 값을 채운다
void GRIDY_IMPLICIT_OBJECT::Fill_density_between_two_vertices(double x0, double y0, double z0, double x1, double y1, double z1, double *arr, double mDefValue)
{
	//bool axis_x, axis_y, axis_z;
	//double *value = &arr;
	//array = &arr;
	int cnt_f, cnt_t;
	double delta_x, delta_y, delta_z;


	//axis_x = aixs_y = axis_z = false;

	//if (MAX(x0, x1) > x1) swap(x0, x1);
	//if (MAX(y0, y1) > y1) swap(y0, y1);
	//if (MAX(z0, z1) > z1) swap(z0, z1);

	// 각 축의 길이
	double len_axis_x = abs(x0 - x1);
	double len_axis_y = abs(y0 - y1);
	double len_axis_z = abs(z0 - z1);



	// 어떤 축이 더 긴지 판별
	double max_axis_value = MAX(MAX(len_axis_x, len_axis_y), len_axis_z);
	
	
	delta_x = ((x1 - x0) > 0) ? (x1 - x0) / max_axis_value : 0.0f;
	delta_y = ((y1 - y0) > 0) ? (y1 - y0) / max_axis_value : 0.0f;
	delta_z = ((z1 - z0) > 0) ? (z1 - z0) / max_axis_value : 0.0f;
	
	/*
	delta_x = len_axis_x / max_axis_value;
	delta_y = len_axis_y / max_axis_value;
	delta_z = len_axis_z / max_axis_value;
	*/


	if (len_axis_x >= max_axis_value)
	{
		//axis_x = true;
		cnt_f = x0;		cnt_t = x1;
	}
	else if (len_axis_y >= max_axis_value)
	{
		//axis_y = true;
		cnt_f = y0;		cnt_t = y1;
	}
	else
	{
		//asix_z = true;
		cnt_f = z0;		cnt_t = z1;
	}

	//printf("%d %d\n", cnt_f, cnt_t);
	//printf("%lf %lf %lf : %lf %lf %lf %lf %lf %lf\n", delta_x, delta_y, delta_z, x0, y0, z0, x1, y1, z1);
	int cnt = 0;
	for (int i = 0; i <= (int)abs(cnt_t - cnt_f); i++)
	{
		arr[IX3((int)(x0 + delta_x * (double)i), (int)(y0 + delta_y * (double)i), (int)(z0 + delta_z * (double)i))] = mDefValue;
		//printf("%lf %lf %lf %lf\n", (x0 + delta_x * i), (y0 + delta_y * i), (z0 + delta_z * i), arr[IX3((int)(x0 + delta_x * i), (int)(y0 + delta_y * i), (int)(z0 + delta_z * i))]);
		cnt++; 
	}
	//printf("# %d\n", cnt);
}








class GRIDY_IMPLICIT_OBJECT_SPHERE
{
public:
	double center_x;
	double center_y;
	double center_z;
	double radius;


public:
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

	int Levelset(int x, int y)
	{
		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2));

		if(dist <= radius) return -1;
		else if(abs(dist-radius) <= 1.0f) return 0;
		else return 1;
	}

	int Levelset(int x, int y, int z)
	{
		int r = 1;

		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

		if(dist < radius) r = -1;
		if(abs(dist-radius) <= 1.74f) r = 0;
		//if (abs(dist - radius) <= 0.5f) r = 0;
		
		return r;
	}

	double Get_distance(int x, int y, int z)
	{
		double dist = sqrt(pow((center_x - x), 2) + pow((center_y - y), 2) + pow((center_z - z), 2));

		return abs(dist-radius);		
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
			norm = sqrt(((double)x-center_x)*((double)x-center_x) + ((double)y-center_y)*((double)y-center_y) + ((double)z-center_z)*((double)z-center_z));

			normal[0] = (x-center_x)/norm;
			normal[1] = (y-center_y)/norm;
			normal[2] = (z-center_z)/norm;
		}
		
		/*
		norm = sqrt(
			((double)(x + 0.5) - center_x) * ((double)(x + 0.5) - center_x) + 
			((double)(y + 0.5) - center_y) * ((double)(y + 0.5) - center_y) +
			((double)(z + 0.5) - center_z) * ((double)(z + 0.5) - center_z)
		);

		normal[0] = ((double)(x + 0.5) - center_x) / norm;
		normal[1] = ((double)(y + 0.5) - center_y) / norm;
		normal[2] = ((double)(z + 0.5) - center_z) / norm;
		*/
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

			printf("!!%f %f %f\n", normal[0], normal[1], normal[2]);
		}
	}
};
















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
		if(n < 0 || n >= mNo)
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
		if(n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Levelset).\n";
			exit(0);
		}
		else
		{
			int levelset = 1;

			double dist = sqrt((mCx[n] - x)*(mCx[n] - x) + (mCy[n] - y)*(mCy[n] - y) + (mCz[n] - z)*(mCz[n] - z));

			if(dist < mRadius[n]) levelset = -1;
			else if(abs(dist-mRadius[n]) < 1.74f) levelset = 0;
		
			return levelset;
		}
	}



	double Get_distance(int n, int x, int y, int z)
	{
		if(n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Get_distance).\n";
			exit(0);
		}
		else
		{
			double dist = sqrt(pow((mCx[n] - x), 2) + pow((mCy[n] - y), 2) + pow((mCz[n] - z), 2));

			return abs(dist-mRadius[n]);
		}
	}



	void Get_normal(int n, int x, int y, int z, double *normal)
	{
		if(n < 0 || n >= mNo)
		{
			cout << "\nGRIDY_IMPLICIT_OBJECT_SPHERES ERROR(Get_normal).\n";
			exit(0);
		}
		else
		{
			double dx = (double)x-mCx[n];
			double dy = (double)y-mCy[n];
			double dz = (double)z-mCz[n];
			double norm = sqrt(dx*dx + dy*dy + dz*dz);

			if(norm != 0.0f)
			{
				normal[0] = dx/norm;
				normal[1] = dy/norm;
				normal[2] = dz/norm;
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
		if(n < 0 || n >= mNo)
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
			if(isGrowing[n] == true && isActivate[n] == true)
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







#endif //_GRIDY_IMPLICIT_OBJECT_H_