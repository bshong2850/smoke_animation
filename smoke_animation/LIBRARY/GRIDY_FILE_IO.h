#ifndef _GRIDY_FILE_IO_H_
#define _GRIDY_FILE_IO_H_



#pragma warning(disable:4996)

#include <stdio.h>
#include <direct.h>
#include <math.h>
#include <io.h>
#include <string>
#include <iostream>

#include "GRIDY_COMMON.h"
#include "GRIDY_VECTOR3F.h"
#include "GRIDY_IMPLICIT_OBJECTS.h"






//	#include <Windows.h>
//	#include <LOG.h>
using namespace std;

class GRIDY_FILE_IO
{
protected :


	int mWidth, mHeight, mDepth;
	int mFile_type;

	static int const GRIDY_MAX_FILE_NUMBER_DIGIT = 5;

	GRIDY_COMMON mGRIDY_COMMON;

public:
	FILE *fp;
	void Open_file_stream(string file_path, int file_number, string file_ext, string mode);
	void Open_file_stream(string full_path_file_name, string mode);
	void Close_file_stream(void);




public :
	void Write_data(float *arr, int width, int height, string file_path, int file_number, string file_ext);
	void Write_data(float *arr, int width, int height, int depth, string file_path, int file_number, string file_ext);
	void Write_data(double *arr, int width, int height, string file_path, int file_number, string file_ext);
	void Write_data(double *arr, int width, int height, int depth, string file_path, int file_number, string file_ext);

	void Write_data_binary(float *arr, int size, string file_path, int file_number, string file_ext);	// made for PBD
	void Write_data_binary(float *arr, int width, int height, int depth, string file_path, int file_number, string file_ext);	// 3D only
	void Write_data_binary(double *arr, int width, int height, int depth, string file_path, int file_number, string file_ext);	// 3D only

	void Read_data(float *arr, string file_path, int file_number, string file_ext);
	void Read_data(double *arr, string file_path, int file_number, string file_ext);
	//void Read_data(GRIDY_SCALAR_FIELD &SF, string file_path, int file_number, string file_ext);
	void Read_data_binary(float *arr, string file_path, int file_number, string file_ext);	// 3D only
	void Read_data_binary(double *arr, string file_path, int file_number, string file_ext);	// 3D only
	void Read_data_binary(float *arr, string file_path, int arr_size, int file_number, string file_ext); // made for PBD
	void Read_data_binary(double *arr, string file_path, int arr_size, int file_number, string file_ext);

	void Write_pbd_simulation_info(string file_path, int vertex_num, int face_num, int start_frame, int end_frame);
	void Read_pbd_simulation_info(string file_path, int &psize, int &fsize, int &start_frame, int &end_frame);
	void Write_pbd_couple_simulation_info(string file_path, int vertex_num, int face_num, int width, int height, int depth, vector3f pbd_st, vector3f pbd_ed, vector3f smoke_st, vector3f smoke_ed,int start_frame, int end_frame);
	void Read_pbd_couple_simulation_info(string file_path, int &psize, int &fsize, int &width, int &height, int &depth, vector3f &pbd_st, vector3f &pbd_ed, vector3f &smoke_st, vector3f &smoke_ed,int &start_frame, int &end_frame);
	
	void Write_simulation_info(string file_path, int width, int height, int start_frame, int end_frame);
	void Write_simulation_info(string file_path, int width, int height, int depth, int start_frame, int end_frame, int file_type);

	void Read_simulation_info(string file_path, int& dimension, int& width, int& height, int& start_frame, int& end_frame);
	void Read_simulation_info(string file_path, int& dimension, int& width, int& height, int& depth, int &start_frame, int &end_frame, int &type_type);	// 3D only, type_type : 0-text, 1-binary

	void Write_simulation_scene_file(string file_path);
	
	void Write_sphere_data(int center_x, int center_y, int center_z, float radius, string file_path, int file_number, string file_ext);
	void Read_sphere_data(int& center_x, int& center_y, int& center_z, float& radius, string file_path, int file_number, string file_ext);

	void Write_particle_data(GRIDY_IMPLICIT_OBJECT_SPHERES& particle, string file_path, int file_number, string file_ext);
	void Read_particle_data(GRIDY_IMPLICIT_OBJECT_SPHERES& particle, string file_path, int file_number, string file_ext);

	//	void put_obj_to_grid(float *arr, char* full_file_name, int width, int height, int depth);
	static int const GRIDY_DATA_FILE_TYPE_BINARY = 1;
	static int const GRIDY_DATA_FILE_TYPE_TEXT = 0;
};



void GRIDY_FILE_IO::Open_file_stream(string file_path, int file_number, string file_ext, string mode)
{
	string file_name = "";
	int full_path_len = 0;
	int file_number_digits = 0;
	errno_t err;
	if(_mkdir(file_path.c_str()) == 0)
	{
		if(mode == "w" || mode == "wb")
			printf("Directory %s is created.\n", file_path.c_str());
		else
		{
			string msg = "Directory";
			msg += file_path.c_str();

			mGRIDY_COMMON.Terminate_system(msg);
		}
	}

	if(file_number<10) file_number_digits = 1;
	else if(file_number<100) file_number_digits = 2;
	else if(file_number<1000) file_number_digits = 3;
	else if(file_number<10000) file_number_digits = 4;
	else file_number_digits = 5;
	// TODO : if MAX_FILE_NUMBER_DIGIT is changed, have to be modited.


	for(int i=0; i<GRIDY_MAX_FILE_NUMBER_DIGIT-file_number_digits; i++)
		file_name += '0';

	file_name += to_string(file_number);

	string full_path_file_name(file_path); 
	full_path_file_name += "\\";
	full_path_file_name += file_name;
	full_path_file_name += ".";
	full_path_file_name += file_ext;


	if((err = fopen_s(&fp, full_path_file_name.c_str(), mode.c_str()))!= 0)
	{
		printf("Cannot open %s\n", full_path_file_name.c_str());
		exit(0);
	}
}





void GRIDY_FILE_IO::Open_file_stream(string full_path_file_name, string mode)
{
	errno_t err;

	if((err = fopen_s(&fp, full_path_file_name.c_str(), mode.c_str()))!= 0)
	{
		printf("Cannot open %s\n", full_path_file_name.c_str());
		exit(0);
	}
}



void GRIDY_FILE_IO::Read_data(float *arr, string file_path, int file_number, string file_ext)
{
	int dimension = 0;
	int width, height, depth;
	int mCount = 0;


	Open_file_stream(file_path, file_number, file_ext, "r");

	fscanf_s(fp, "%d", &dimension);
	

	if(dimension == 2)
	{
		fscanf_s(fp, "%d %d", &width, &height);
		
		for(int i=1; i<=width; i++)
			for(int j=1; j<=height; j++)
			{
				fscanf_s(fp, "%f", &arr[(i)+((height+2)*j)]);
				mCount++;
			}
	}
	else
	{
		fscanf_s(fp, "%d %d %d", &width, &height, &depth);

		for(int i=0; i<width+2; i++)
		{
			for(int j=0; j<height+2; j++)
			{
				for(int k=0; k<depth+2; k++)
				{
					fscanf(fp, "%f", &arr[(i)+((width+2)*j)+(width+2)*(height+2)*(k)]);
					mCount++;
				}
			}
		}

	}
	printf("\t%d data read.\n", mCount);
			
	Close_file_stream();	
}




void GRIDY_FILE_IO::Read_data(double *arr, string file_path, int file_number, string file_ext)
{
	int dimension = 0;
	int width, height, depth;
	int mCount = 0;


	Open_file_stream(file_path, file_number, file_ext, "r");

	fscanf_s(fp, "%d", &dimension);


	if (dimension == 2)
	{
		fscanf_s(fp, "%d %d", &width, &height);

		for (int i = 1; i <= width; i++)
			for (int j = 1; j <= height; j++)
			{
				fscanf_s(fp, "%f", &arr[(i)+((height + 2)*j)]);
				mCount++;
			}
	}
	else
	{
		fscanf_s(fp, "%d %d %d", &width, &height, &depth);

		for (int i = 0; i<width + 2; i++)
		{
			for (int j = 0; j<height + 2; j++)
			{
				for (int k = 0; k<depth + 2; k++)
				{
					fscanf(fp, "%f", &arr[(i)+((width + 2)*j) + (width + 2)*(height + 2)*(k)]);
					mCount++;
				}
			}
		}

	}
	printf("\t%d data read.\n", mCount);

	Close_file_stream();
}





void GRIDY_FILE_IO::Read_data_binary(float *arr, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "rb");

	printf("%d %d %d\n", mWidth, mHeight, mDepth);

	fread(arr, sizeof(float)*(mWidth+2)*(mHeight+2)*(mDepth+2), 1, fp);

	Close_file_stream();
}




void GRIDY_FILE_IO::Read_data_binary(double *arr, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "rb");
	fread(arr, sizeof(double)*(mWidth+2)*(mHeight+2)*(mDepth+2), 1, fp);
	Close_file_stream();
}



void GRIDY_FILE_IO::Read_data_binary(float *arr, string file_path, int arr_size, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "rb");

	fread(arr, sizeof(float)*(arr_size), 1, fp);

	Close_file_stream();
}





void GRIDY_FILE_IO::Read_data_binary(double *arr, string file_path, int arr_size, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "rb");

	fread(arr, sizeof(double)*(arr_size), 1, fp);

	Close_file_stream();
}





/*
// NOT USE YET
void GRIDY_FILE_IO::Read_data(GRIDY_SCALAR_FIELD &SF, string file_path, int file_number, string file_ext)
{
	int dimension = 0;
	int width, height, depth;
	int mCount = 0;

	Open_file_stream(file_path, file_number, file_ext, "r");
	//printf("open_file_stream completed.\n");
	fscanf(fp, "%d", &dimension);

	if(dimension == 2)
	{
		fscanf_s(fp, "%d %d", &width, &height);
		
		for(int i=1; i<=width; i++)
			for(int j=1; j<=height; j++)
			{
				fscanf_s(fp, "%f", &SF.value[(i)+((height+2)*j)]);
				mCount++;
			}
	}
	else
	{
		fscanf(fp, "%d %d %d", &width, &height, &depth);
		
		for(int i=0; i<width+2; i++)
		{
			for(int j=0; j<height+2; j++)
			{
				for(int k=0; k<depth+2; k++)
				{
					fscanf(fp, "%f", &SF.value[(i)+((width+2)*j)+(width+2)*(height+2)*(k)]);
					mCount++;
				}
			}
		}
	}
	printf("%d data read\n", mCount);
	Close_file_stream();
}
*/



void GRIDY_FILE_IO::Write_data(float *arr, int width, int height, string file_path, int file_number, string file_ext)
{
	int mCount = 0;

	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d %d %d 0\n", 2, width, height);	// 2 is dimension

	for(int i=1; i<=width; i++)
	{
		for(int j=1; j<=height; j++)
		{
			fprintf(fp, "%f\t", arr[(i)+((height+2)*j)]);
			mCount++;
		}
		
		fprintf(fp, "\n");
	}

	//printf("\t%d data were written.\n", mCount);

	Close_file_stream();

}




void GRIDY_FILE_IO::Write_data(double *arr, int width, int height, string file_path, int file_number, string file_ext)
{
	int mCount = 0;

	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d %d %d 0\n", 2, width, height);	// 2 is dimension

	for (int i = 1; i <= width; i++)
	{
		for (int j = 1; j <= height; j++)
		{
			fprintf(fp, "%lf\t", arr[(i)+((height + 2)*j)]);
			mCount++;
		}

		fprintf(fp, "\n");
	}

	//printf("\t%d data were written.\n", mCount);

	Close_file_stream();

}




void GRIDY_FILE_IO::Write_data(float *arr, int width, int height, int depth, string file_path, int file_number, string file_ext)
{
	int mCount = 0;

	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d %d %d %d\n", 3, width, height, depth);	// 3 is dimension

	for(int i=0; i<width+2; i++)
	{
		for(int j=0; j<height+2; j++)
		{
			for(int k=0; k<depth+2; k++)
			{
				fprintf(fp, "%f\t", arr[(i)+((width+2)*j+((width+2)*(height+2)*k))]);

				mCount++;
			}
			fprintf(fp, "\n");
		}
	}

	Close_file_stream();
}




void GRIDY_FILE_IO::Write_data(double *arr, int width, int height, int depth, string file_path, int file_number, string file_ext)
{
	int mCount = 0;

	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d %d %d %d\n", 3, width, height, depth);	// 3 is dimension

	for (int i = 0; i<width + 2; i++)
	{
		for (int j = 0; j<height + 2; j++)
		{
			for (int k = 0; k<depth + 2; k++)
			{
				fprintf(fp, "%f\t", arr[(i)+((width + 2)*j + ((width + 2)*(height + 2)*k))]);

				mCount++;
			}
			fprintf(fp, "\n");
		}
	}

	Close_file_stream();
}



void GRIDY_FILE_IO::Write_data_binary(float *arr, int width, int height, int depth, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "wb");

	fwrite(arr, sizeof(float)*(width+2)*(height+2)*(depth+2), 1, fp);

	Close_file_stream();
}

void GRIDY_FILE_IO::Write_data_binary(double *arr, int width, int height, int depth, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "wb");

	//fwrite(arr, sizeof(double)*(width + 2)*(height + 2)*(depth + 2), 1, fp);
	fwrite(arr, sizeof(double)*(width+2)*(height+2)*(depth+2), 1, fp);

	Close_file_stream();
}

void GRIDY_FILE_IO::Write_data_binary(float *arr, int size, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "wb");

	fwrite(arr, sizeof(float)*size, 1, fp);

	Close_file_stream();
}


void GRIDY_FILE_IO::Write_sphere_data(int center_x, int center_y, int center_z, float radius, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d %d %d %f", center_x, center_y, center_z, radius);

	Close_file_stream();
}



void GRIDY_FILE_IO::Read_sphere_data(int& center_x, int& center_y, int& center_z, float& radius, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "r");

	fscanf(fp, "%d %d %d %f", &center_x, &center_y, &center_z, &radius);

	Close_file_stream();
}




void GRIDY_FILE_IO::Write_particle_data(GRIDY_IMPLICIT_OBJECT_SPHERES& particle, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "w");

	fprintf(fp, "%d\n", particle.mNo);

	for (int l = 0; l<particle.mNo; l++)
	{
		if (particle.isActivate[l] == true)
			fprintf(fp, "%f %f %f %f\n", particle.mCx[l], particle.mCy[l], particle.mCz[l], particle.mRadius[l]);
	}

	Close_file_stream();
}





void GRIDY_FILE_IO::Read_particle_data(GRIDY_IMPLICIT_OBJECT_SPHERES& particle, string file_path, int file_number, string file_ext)
{
	Open_file_stream(file_path, file_number, file_ext, "r");
	
	fscanf(fp, "%d\n", &particle.mNo);
	particle.Init(particle.mNo);
	
	for (int l = 0; l<particle.mNo; l++)
	{
//		fscanf(fp, "%lf %lf %lf %lf\n", &particle.mCx[l], &particle.mCy[l], &particle.mCz[l], &particle.mRadius[l]);
//		double uu,vv,ww;
//		fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &particle.mCx[l], &particle.mCy[l], &particle.mCz[l], &uu, &vv, &ww,&particle.mNx[l], &particle.mNy[l], &particle.mNz[l], &particle.mAngle[l],&particle.mRadius[l], &particle.mTemperature[l]);
		
		//	original
		fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &particle.mCx[l], &particle.mCy[l], &particle.mCz[l], &particle.mNx[l], &particle.mNy[l], &particle.mNz[l], &particle.mAngle[l], &particle.mRadius[l], &particle.mTemperature[l]);
		
//		fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &particle.mCx[l], &particle.mCy[l], &particle.mCz[l], &particle.mNx[l], &particle.mNy[l], &particle.mNz[l]);
//		particle.mRadius[l] = 0.5;
//		particle.mAngle[l] = 0;
//		particle.mTemperature[l] = 3000;
//		printf("\t%lf %lf %lf %lf\n", particle.mCx[l], particle.mCy[l], particle.mCz[l], particle.mRadius[l]);
	}

	Close_file_stream();
}




void GRIDY_FILE_IO::Write_simulation_info(string file_path, int width, int height, int start_frame, int end_frame)
{
	Open_file_stream(file_path, 0, "GRIDY_INIT", "w");

	fprintf(fp, "%d %d %d 0\n", 2, width, height);	// 2 is dimension
	fprintf(fp, "%d %d", start_frame, end_frame);	

	Close_file_stream();
}
void GRIDY_FILE_IO::Write_pbd_simulation_info(string file_path, int vertex_num, int face_num, int start_frame, int end_frame)
{
	Open_file_stream(file_path, 0, "GRIDY_PBD_INIT", "w");

	fprintf(fp, "%d %d\n", vertex_num, face_num);	// 2 is dimension
	fprintf(fp, "%d %d", start_frame, end_frame);	

	Close_file_stream();
}




void GRIDY_FILE_IO::Write_simulation_info(string file_path, int width, int height, int depth, int start_frame, int end_frame, int file_type)
{
	Open_file_stream(file_path, 0, "GRIDY_INIT", "w");

	fprintf(fp, "%d %d %d %d\n", 3, width, height, depth);	// 3 is dimension
	fprintf(fp, "%d %d\n", start_frame, end_frame);	
	fprintf(fp, "%d", file_type);

	Close_file_stream();
}


void GRIDY_FILE_IO::Read_simulation_info(string file_path, int &dimension, int &width, int &height, int &start_frame, int &end_frame)
{
	printf("Reading INIT file from %s\n", file_path.c_str());


	Open_file_stream(file_path, 0, "GRIDY_INIT", "r");

	
	fscanf_s(fp, "%d %d %d\n", &dimension, &width, &height);	
	fscanf_s(fp, "%d %d", &start_frame, &end_frame);	

	Close_file_stream();
}
void GRIDY_FILE_IO::Read_pbd_simulation_info(string file_path, int &psize, int &fsize, int &start_frame, int &end_frame)
{
	printf("Reading INIT file from %s\n", file_path.c_str());


	Open_file_stream(file_path, 0, "GRIDY_PBD_INIT", "r");

	
	fscanf_s(fp, "%d %d\n", &psize, &fsize);
	fscanf_s(fp, "%d %d", &start_frame, &end_frame);	

	Close_file_stream();
}



void GRIDY_FILE_IO::Read_simulation_info(string file_path, int &dimension, int &width, int &height, int &depth, int &start_frame, int &end_frame, int &file_type)
{
	printf("Reading INIT file from %s\n", file_path.c_str());

	Open_file_stream(file_path, 0, "GRIDY_INIT", "r");
	
	fscanf_s(fp, "%d %d %d %d\n", &dimension, &width, &height, &depth);	
	fscanf_s(fp, "%d %d", &start_frame, &end_frame);	
	fscanf_s(fp, "%d", &file_type);	

	mWidth = width;
	mHeight = height;
	mDepth = depth;
	mFile_type = file_type;

	Close_file_stream();
}


void GRIDY_FILE_IO::Write_pbd_couple_simulation_info(string file_path, int vertex_num, int face_num, int width, int height, int depth, vector3f pbd_st, vector3f pbd_ed, vector3f smoke_st, vector3f smoke_ed,int start_frame, int end_frame)
{
	
	Open_file_stream(file_path, 0, "GRIDY_PBD_COUPLE_INIT", "w");
	
	fprintf(fp,"%d %d\n",vertex_num, face_num);
	fprintf(fp,"%d %d %d\n",width,height,depth);

	fprintf(fp,"%.9lf %.9lf %.9lf  %.9lf %.9lf %.9lf\n",pbd_st.x,pbd_st.y,pbd_st.z,pbd_ed.x,pbd_ed.y,pbd_ed.z);
	fprintf(fp,"%.9lf %.9lf %.9lf  %.9lf %.9lf %.9lf\n",smoke_st.x,smoke_st.y,smoke_st.z,smoke_ed.x,smoke_ed.y,smoke_ed.z);

	fprintf(fp,"%d %d\n",start_frame, end_frame);
	Close_file_stream();

}
void GRIDY_FILE_IO::Read_pbd_couple_simulation_info(string file_path, int &psize, int &fsize, int &width, int &height, int &depth, vector3f &pbd_st, vector3f &pbd_ed, vector3f &smoke_st, vector3f &smoke_ed,int &start_frame, int &end_frame)
{
	printf("Reading INIT file from %s\n", file_path.c_str());
	Open_file_stream(file_path, 0, "GRIDY_PBD_COUPLE_INIT", "r");

	fscanf_s(fp, "%d %d",&psize, &fsize);
	fscanf_s(fp, "%d %d %d",&width, &height, &depth);
	fscanf_s(fp, "%lf %lf %lf  %lf %lf %lf\n",&pbd_st.x,&pbd_st.y,&pbd_st.z,&pbd_ed.x,&pbd_ed.y,&pbd_ed.z);
	fscanf_s(fp, "%lf %lf %lf  %lf %lf %lf\n",&smoke_st.x,&smoke_st.y,&smoke_st.z,&smoke_ed.x,&smoke_ed.y,&smoke_ed.z);

	fscanf_s(fp, "%d %d",&start_frame, &end_frame);
	Close_file_stream();
}




void GRIDY_FILE_IO::Close_file_stream(void)
{
	if(fp) fclose(fp);
}





#endif