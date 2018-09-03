#ifndef _GRIDY_SIGNED_DISTANCE_FIELD_H_
#define _GRIDY_SIGNED_DISTANCE_FIELD_H_


#include "GRIDY_SCALAR_FIELD.h"
#include "GRIDY_ADT.h"
#include "GRIDY_COMMON.h"


using namespace std;

class GRIDY_SIGNED_DISTANCE_FIELD : GRIDY_SCALAR_FIELD
{
public:
	GRIDY_SIGNED_DISTANCE_FIELD(void) {};
	~GRIDY_SIGNED_DISTANCE_FIELD(void) 
	{
		Terminate();
	};

	//void Make_SDF_from_obj(GRIDY_SCALAR_FIELD SDF, string full_file_name, int width, int height, int depth);
	void Make_SDF_from_obj(string full_file_name, int width, int height, int depth, double resize_ratio = 1.0f, int margin0 = 0, int margin1 = 0, int margin2 = 0, int option = 0); // option은 임시
	void Input_SDF(string full_file_name, int width, int height, int depth);
	void Input_SDF(string full_file_name, int width, int height, int depth, double ratio);

	double *phi;
	double **normal;
	short *sign;

	static int const GRIDY_UNKNOWN_GRID = 999999;
	static short const GRIDY_SDF_INTERFACE = 0;
	static short const GRIDY_SDF_INNER = -1;
	static short const GRIDY_SDF_OUTER = 1;

	inline double Get_grid_cell_size(void) { return mGrid_cell_size; }

private:
	GRIDY_COMMON mUTIL;

	int mWidth, mHeight, mDepth;
	string mOBJ_path_and_name;
	int mSize;
	double mGrid_cell_size;
	double mResize_ratio;	// whether resize needed. ratio = 1/original size
	int mMargin_on_grid[3];

	void Init(int width, int height, int depth);
	void Terminate(void);
	void Clear_data(double default_value = 0.0f);

	void Find_bounding_box(double* bounding_box_AA, double* bounding_box_BB);
	void Find_nearist_point_from_obj(void);
	void Find_interface_and_allocate_sign(void);
	void Loop_over_the_unknown_grid_points(void);
	void Set_interface_btw_edges(double *edgeAA, double *edgeBB, double *edgeCC, double *delta_grid, double *bounding_box_AA, double *margin);
	void Set_interface_btw_points(double *edgeAA, double *edgeBB, double *delta_grid, double *bounding_box_AA, double *margin);

	void Make_unsigned_distance_field(double *phi, double h);
	void Inner_flood_fill(int x, int y, int z);
	//void Find_nearist_point_from_obj(double *phi, double **normal);
	//void Find_interface_and_allocate_sign(short *sign, double *phi);
	//void Loop_over_the_unknown_grid_points(short *sign, double *phi);

	inline double MAX3(double a, double b, double c)
	{
		double r = a;

		if(r < b) r = b;
		if(r < c) r = c;

		return r;
	}

	inline double MIN3(double a, double b, double c)
	{
		double r = a;

		if(r > b) r = b;
		if(r > c) r = c;

		return r;
	}

	
};




///////사용//////////////
void GRIDY_SIGNED_DISTANCE_FIELD::Init(int width, int height, int depth)
{
	//cout << "Init 시작" << endl;
	mSize = (width+2)*(height+2)*(depth+2);
	mResize_ratio = 1.0f;
	mWidth = width;
	mHeight = height;
	mDepth = depth;

	mMargin_on_grid[0] = mMargin_on_grid[1] = mMargin_on_grid[2] = 0;

	//	allocate array
	phi = new double[mSize];
	sign = new short[mSize];
	normal = new double*[mSize];

	for(int i=0; i<mSize; i++)
		normal[i] = new double[3];


	Clear_data();
	//cout << "Init 종료" << endl;
}




void GRIDY_SIGNED_DISTANCE_FIELD::Terminate(void)
{
	if(phi) delete[] phi;
	if(sign) delete[] sign;
	if(normal) {
		for(int i=0; i<mSize; i++)
			delete[] normal[i];

		delete[] normal;
	}
}




///////사용//////////////
void GRIDY_SIGNED_DISTANCE_FIELD::Clear_data(double default_value)
{
	#pragma omp parallel for
	for(int i=0; i<mSize; i++)
	{
		phi[i] = (double)GRIDY_UNKNOWN_GRID;
		normal[i][0] = normal[i][1] = normal[i][2] = 0.0f;
		sign[i] = GRIDY_SDF_OUTER;
	}
}




void GRIDY_SIGNED_DISTANCE_FIELD::Input_SDF(string full_file_name, int width, int height, int depth)
{
	FILE *fp;												//	file pointer
	errno_t err;



	Init(width, height, depth);
	mOBJ_path_and_name = full_file_name;



	//	open file stream
	if((err = fopen_s(&fp, mOBJ_path_and_name.c_str(), "r"))!= 0)
	{
		printf("Cannot open %s\n", mOBJ_path_and_name.c_str());
		exit(0);
	}


	//	input resolution of the sdf file
	int sdf_width, sdf_height, sdf_depth;
	fscanf(fp, "%d %d %d", &sdf_width, &sdf_height, &sdf_depth);
		
	if((sdf_width > width) || (sdf_height > height) || (sdf_depth > depth))
	{
		printf("%d %d %d > %d %d %d\n", sdf_width, sdf_height, sdf_depth, width, height, depth);
		printf("SDF file resolution error.\n");
		exit(0);
	}

	
	
	//	input bounding position and grid cell size
	double sdf_bounding_x, sdf_bounding_y, sdf_bounding_z;
//	double sdf_grid_cell_size;
	
	fscanf(fp, "%lf %lf %lf", &sdf_bounding_x, &sdf_bounding_y, &sdf_bounding_z);
	fscanf(fp, "%lf", &mGrid_cell_size);
	


	for(int k=0; k<sdf_depth; k++)
	{
		for(int j=0; j<sdf_height; j++)
		{
			for(int i=0; i<sdf_width; i++)
			{
				fscanf(fp, "%lf", &phi[IX3(i,j,k)]);
				//printf("%lf", phi[IX3(i, j, k)]);
			}
		}
	}

	
	
	fclose(fp);
}





void GRIDY_SIGNED_DISTANCE_FIELD::Input_SDF(string full_file_name, int width, int height, int depth, double ratio)
{
	FILE *fp;												//	file pointer
	errno_t err;



	Init(width, height, depth);
	mOBJ_path_and_name = full_file_name;



	//	open file stream
	if((err = fopen_s(&fp, mOBJ_path_and_name.c_str(), "r"))!= 0)
	{
		printf("Cannot open %s\n", mOBJ_path_and_name.c_str());
		exit(0);
	}


	//	input resolution of the sdf file
	int sdf_width, sdf_height, sdf_depth;
	fscanf(fp, "%d %d %d", &sdf_width, &sdf_height, &sdf_depth);
		
	if((sdf_width*ratio > width) || (sdf_height*ratio > height) || (sdf_depth*ratio > depth))
	{
		printf("%d %d %d > %d %d %d\n", sdf_width, sdf_height, sdf_depth, width, height, depth);
		printf("SDF file resolution error.\n");
		exit(0);
	}

	
	
	//	input bounding position and grid cell size
	double sdf_bounding_x, sdf_bounding_y, sdf_bounding_z;

	
	fscanf(fp, "%f %f %f", &sdf_bounding_x, &sdf_bounding_y, &sdf_bounding_z);
	fscanf(fp, "%f", &mGrid_cell_size);
	
	printf("%d %d %d\n", sdf_width, sdf_height, sdf_depth);

	for(int k=0; k<sdf_depth; k++)
	{
		for(int j=0; j<sdf_height; j++)
		{
			for(int i=0; i<sdf_width; i++)
			{
				int ti = (int)(i*ratio);
				int tj = (int)(j*ratio);
				int tk = (int)(k*ratio);

				double tphi = 0.0f;
				fscanf(fp, "%f", &tphi);
				if(tphi != 9.6f) printf("%f ", tphi);
				//tphi *= ratio;

				//if(phi[IX3(ti, tj, tk)] > tphi) phi[IX3(ti, tj, tk)] = tphi;
				phi[IX3(ti, tj, tk)] = tphi;

				//if(tphi < 0.96f)	printf("%d %d %d %f\n", ti, tj, tk, tphi);
			}
		}
	}

	
	
	fclose(fp);
}


///////사용//////////////
void GRIDY_SIGNED_DISTANCE_FIELD::Make_SDF_from_obj(string full_file_name, int width, int height, int depth, double resize_ratio, int margin0, int margin1, int margin2, int option)
{
	if(resize_ratio > 1.0f) mUTIL.Terminate_system("Resize ratio is must be larger than 1.0f.");


	Init(width, height, depth);

	mOBJ_path_and_name = full_file_name;
	mResize_ratio = resize_ratio;
	mMargin_on_grid[0] = margin0;	mMargin_on_grid[1] = margin1;	mMargin_on_grid[2] = margin2;	


	
	Find_nearist_point_from_obj();
	if(option == 1) Inner_flood_fill(width/2, height/3, depth/2);
	if(option == 2) Make_unsigned_distance_field(phi, 1.0f);	// doesn't work well

	//Find_interface_and_allocate_sign();
	//Find_nearist_point_from_obj(phi, normal);
	////Set_interface_along_edge(sign, phi);
	//Find_interface_and_allocate_sign(sign, phi);
	//Loop_over_the_unknown_grid_points();
}





void GRIDY_SIGNED_DISTANCE_FIELD::Find_nearist_point_from_obj(void)
{
	char id[256];											//	id (v, e, f ...)
	char buf_to_drop_string[256];							//	a buffer for dropping useless strings
	double vertex_position[3];								//	position of vertex where read from obj file
	double face_position[3];								//	position of face where read from obj file
	double vertex_normal[3];
	double bounding_box_AA[3], bounding_box_BB[3];			//	boundary of obj file; AA is the smallest values, BB is the largest values
	double edge_length[3];									//	the length of each edges (x, y, z axis)
	double delta_grid[3];									//	the size of a cell
	double margin[3];

	int no_faces = 0;
	int no_vertexs = 0;

	double **array_vertex;
	int **array_face;

	FILE *fp;												//	file pointer
	errno_t err;

	double const LARGEST_VALUE = 9999999.0f, SMALLEST_VALUE = -9999999.0f;



	//	init arr mAA, mBB
	for(int i=0; i<3; i++)
	{
		bounding_box_AA[i] = LARGEST_VALUE;
		bounding_box_BB[i] = SMALLEST_VALUE;
	}

	//	open file stream
	if((err = fopen_s(&fp, mOBJ_path_and_name.c_str(), "r"))!= 0)
	{
		printf("Cannot open %s\n", mOBJ_path_and_name.c_str());
		exit(0);
	}
		
	//	set bounding box
	while(!feof(fp))
	{
		//fscanf(fp, "%s %f %f %f", id, &vertex_position[0], &vertex_position[1], &vertex_position[2]);
		fscanf(fp, "%s", id);

		if(id[0] == '#')	// # have to drop next strings
		{
			fgets(buf_to_drop_string, 255, fp);
		}
		else if(id[0] == 'v' && id[1] != 'n')	// v is vertex, and vn is vertex normal
		{
			no_vertexs++;

			fscanf(fp, "%lf %lf %lf", &vertex_position[0], &vertex_position[1], &vertex_position[2]);

			//cout << vertex_position[0] << vertex_position[1] << vertex_position[2] << endl;

			for(int i=0; i<3; i++)
			{
				if(bounding_box_AA[i] > vertex_position[i]) bounding_box_AA[i] = vertex_position[i];		// make AA get the smallest value
				if(bounding_box_BB[i] < vertex_position[i]) bounding_box_BB[i] = vertex_position[i];	// make BB get the largest value

				//cout << "bounding_box_AA : " << bounding_box_AA[i] << endl;
				//cout << "bounding_box_BB : " << bounding_box_BB[i] << endl;
			}
		}
		else if(id[0] == 'f') no_faces++;
		else
		{
			fgets(buf_to_drop_string, 255, fp);
		}
	}
	//cout << "no_faces : " << no_faces << endl << "no_vertexs : " << no_vertexs << endl;


	//	----------------------------
	//	reload the obj file
	//	----------------------------
	rewind(fp);



	//	allocate arrays for save vertex and face info
	array_vertex = new double*[no_vertexs];
	array_face = new int*[no_faces];

	for(int i=0; i<no_vertexs; i++)
		array_vertex[i] = new double[3];

	for(int i=0; i<no_faces; i++)
		array_face[i] = new int[3];



	//	get the length of each edge
	for(int i=0; i<3; i++)
	{
		edge_length[i] = margin[i] = (double)abs(bounding_box_BB[i]-bounding_box_AA[i]);

		//cout << edge_length[i] << endl;

		edge_length[i] /= mResize_ratio;
		margin[i] -= edge_length[i];
		margin[i] /= 2.0f;
	}
	
	delta_grid[0] = (double)edge_length[0]/mWidth;
	delta_grid[1] = (double)edge_length[1]/mHeight;
	delta_grid[2] = (double)edge_length[2]/mDepth;


	/*cout << "edge_length[0] : " << edge_length[0] << endl;
	cout << "edge_length[1] : " << edge_length[1] << endl;
	cout << "edge_length[2] : " << edge_length[2] << endl;
	cout << "delta_grid[0] : " << delta_grid[0] << endl;
	cout << "delta_grid[1] : " << delta_grid[1] << endl;
	cout << "delta_grid[2] : " << delta_grid[2] << endl;*/

	//	set delta for each edges to keep ratio of obj file
	double max_delta_grid = MAX3(delta_grid[0], delta_grid[1], delta_grid[2]);
	//delta_grid[0] = delta_grid[1] = delta_grid[2] = max_delta_grid;


	int cnt_vertex_no = 0, cnt_face_no = 0;

	//	put DefValue on grid
	while(!feof(fp))
	{
		//cout << "cnt_vertex_no : " << cnt_vertex_no << endl; // test
		char temp_input[3][100];
		fscanf(fp, "%s %s %s %s", id, temp_input[0], temp_input[1], temp_input[2]);

		if(id[0] == '#')	// have to drop next strings
		{
			fgets(buf_to_drop_string, 255, fp);
		}
		else if(id[0] == 'v')
		{
			//fscanf(fp, "%f %f %f", &vertex_position[0], &vertex_position[1], &vertex_position[2]);
			for (int i = 0; i < 3; i++)
			{
				//cout << temp_input[i] << endl;
				vertex_position[i] = (double)atof(temp_input[i]);
			}

			if (id[1] == 'n') //	vertex normal
			{
				for (int i = 0; i<3; i++)
					vertex_normal[i] = vertex_position[i];
			}
			else //	vertex
			{
				//	save vertex position to vertex_array
				for(int i=0; i<3; i++)
					array_vertex[cnt_vertex_no][i] = vertex_position[i];
				
				cnt_vertex_no++;


				//	get vertex_position_on_grid and phi value
				double vertex_position_on_grid[3], distance_vertex_btw_grid[3], temp_phi, temp_phi0;

				for(int i=0; i<3; i++)
				{
					vertex_position_on_grid[i] = abs((vertex_position[i] - bounding_box_AA[i] - margin[i]) / delta_grid[i]);
					distance_vertex_btw_grid[i] = (vertex_position_on_grid[i] - (int)vertex_position_on_grid[i]);

					vertex_position_on_grid[i] += mMargin_on_grid[i];
					//cout << "vertex_position_on_grid[" << i << "] : " << vertex_position_on_grid[i] << endl;
				}

				temp_phi0 = abs(phi[IX3((int)vertex_position_on_grid[0], (int)vertex_position_on_grid[1], (int)vertex_position_on_grid[2])]);
				temp_phi = sqrt(distance_vertex_btw_grid[0]*distance_vertex_btw_grid[0] + distance_vertex_btw_grid[1]*distance_vertex_btw_grid[1] + distance_vertex_btw_grid[2]*distance_vertex_btw_grid[2]);
			
				if(temp_phi < temp_phi0) 
				{
					phi[IX3((int)vertex_position_on_grid[0], (int)vertex_position_on_grid[1], (int)vertex_position_on_grid[2])] = temp_phi;
					sign[IX3((int)vertex_position_on_grid[0], (int)vertex_position_on_grid[1], (int)vertex_position_on_grid[2])] = GRIDY_SDF_INTERFACE;

					double norm_vertex_normal = sqrt(vertex_normal[0]*vertex_normal[0] + vertex_normal[1]*vertex_normal[1] + vertex_normal[2]*vertex_normal[2]);

					for(int i=0; i<3; i++)
						normal[IX3((int)vertex_position_on_grid[0], (int)vertex_position_on_grid[1], (int)vertex_position_on_grid[2])][i] = vertex_normal[i]/norm_vertex_normal;
				}
			}
		} 
		////////각 vertex에 대해서 grid로 나눈 곳에 sign을 주었음 inter값으로
		
		// face
		else if(id[0] == 'f')
		{
			string str_vertex_no[3];
			
			for(int i=0; i<3; i++)
			{
				str_vertex_no[i] = temp_input[i];

				for (int i = 0; i<3; i++)
					face_position[i] = (double)atof(temp_input[i]);
					array_face[cnt_face_no][i] = face_position[i];

				//scanf("%d", &array_face[cnt_face_no][i]);

				/*
				// 만약 중간에 //가 있는 형태라면
				while((cut_at = str_vertex_no[i].find_first_of("//")) != str_vertex_no[i].npos)
				{
					if(cut_at > 0)
					{
						str_vertex_no[i] = str_vertex_no[i].substr(0, cut_at);
						array_face[cnt_face_no][i] = atoi(str_vertex_no[i].c_str());

						break;
					}
				}
				*/

				/*
				// 주석
				if(str_vertex_no[i].find("//") == string::npos) 
				{
					array_face[cnt_face_no][i] = atof(temp_input[i]);
				}
				*/
			}
			cnt_face_no++;
		}
		else fgets(buf_to_drop_string, 255, fp);
	}

	fclose(fp);

	//cout << "test " << endl;
	

	//	find interface along edges
	double pointA[3], pointB[3], pointC[3];
	

	//printf("여기까지는 ok\n");
	for(int i=0; i<no_faces; i++)
	{
		for(int j=0; j<3; j++)
		{
			pointA[j] = array_vertex[array_face[i][0]-1][j];
			pointB[j] = array_vertex[array_face[i][1]-1][j];
			pointC[j] = array_vertex[array_face[i][2]-1][j];
			/*cout << "pointA[" << j << "] : " << pointA[j] << endl;
			cout << "pointB[" << j << "] : " << pointB[j] << endl;
			cout << "pointC[" << j << "] : " << pointC[j] << endl;*/
		}
		//cout << "한면 봤다" << endl;
		/*
		Set_interface_btw_points(pointA, pointB, delta_grid, bounding_box_AA, margin);
		Set_interface_btw_points(pointB, pointC, delta_grid, bounding_box_AA, margin);
		Set_interface_btw_points(pointA, pointC, delta_grid, bounding_box_AA, margin);
		*/
		/*cout << "delta_grid[0] : " << delta_grid[0] << endl;
		cout << "delta_grid[1] : " << delta_grid[1] << endl;
		cout << "delta_grid[2] : " << delta_grid[2] << endl;
		cout << "bounding_box_AA[0] : " << bounding_box_AA[0] << endl;
		cout << "bounding_box_AA[1] : " << bounding_box_AA[1] << endl;
		cout << "bounding_box_AA[2] : " << bounding_box_AA[2] << endl;*/
		Set_interface_btw_edges(pointA, pointB, pointC, delta_grid, bounding_box_AA, margin);

			/*
			for(int k=0; k<3; k++)
			{
				// get edge info
				if(j == 2)
				{
					edgeAA[k] = array_vertex[array_face[i][j]-1][k];		
					edgeBB[k] = array_vertex[array_face[i][0]-1][k];		
				}
				else
				{
					edgeAA[k] = array_vertex[array_face[i][j]-1][k];	
					edgeBB[k] = array_vertex[array_face[i][j+1]-1][k];
				}
			}
			*/
		

		
			
			/*
			//double edge_distance = (int)(sqrt((edgeBB_on_grid[0]-edgeAA_on_grid[0])*(edgeBB_on_grid[0]-edgeAA_on_grid[0]) + (edgeBB_on_grid[1]-edgeAA_on_grid[1])*(edgeBB_on_grid[1]-edgeAA_on_grid[1]) + (edgeBB_on_grid[2]-edgeAA_on_grid[2])*(edgeBB_on_grid[2]-edgeAA_on_grid[2])))+1;
			double edge_distance = (sqrt((edgeBB[0]-edgeAA[0])*(edgeBB[0]-edgeAA[0]) + (edgeBB[1]-edgeAA[1])*(edgeBB[1]-edgeAA[1]) + (edgeBB[2]-edgeAA[2])*(edgeBB[2]-edgeAA[2])));
			//edge_distance = 1.0f;


			double step_size = edge_distance / min_delta_grid;
			//printf("edge_distance : %f, min_delta_grid : %f, step size : %f\n", edge_distance, min_delta_grid, step_size);

			for(int k=0; k<3; k++)
				edge_step_size[k] = (edgeBB[k] - edgeAA[k])/step_size;
				//edge_step_size[k] = edgeBB[k] - edgeAA[k]/edge_distance;
				//edge_step_size[k] = abs(edgeBB[k] - edgeAA[k])/step_size;



			double current_position[3];
			int current_position_on_grid[3];
		
			current_position[0] = edgeAA[0];	current_position[1] = edgeAA[1];	current_position[2] = edgeAA[2];

			printf("%f %f / %f %f / %f %f\n", abs(current_position[0] - edgeBB[0]), abs(edge_step_size[0]), abs(current_position[1] - edgeBB[1]), abs(edge_step_size[1]), abs(current_position[2] - edgeBB[2]), abs(edge_step_size[2]));

			//while(current_position[0] <= edgeBB[0] && current_position[1] <= edgeBB[1] && current_position[2] <= edgeBB[2])
			while(abs(current_position[0] - edgeBB[0]) >= abs(edge_step_size[0]) && abs(current_position[1] - edgeBB[1]) >= abs(edge_step_size[1]) && abs(current_position[2] - edgeBB[2]) >= abs(edge_step_size[2]))
			{
				for(int k=0; k<3; k++)
					current_position_on_grid[k] = abs((current_position[k] - bounding_box_AA[k] - margin[k])/delta_grid[k]+1);

				sign[IX3((int)current_position_on_grid[0], (int)current_position_on_grid[1], (int)current_position_on_grid[2])] = GRIDY_SDF_INTERFACE;
				
				printf("\t%d : %d %d %d, %f %f %f, %d\n", i, current_position_on_grid[0], current_position_on_grid[1], current_position_on_grid[2], edge_step_size[0], edge_step_size[1], edge_step_size[2], sign[IX3((int)current_position_on_grid[0], (int)current_position_on_grid[1], (int)current_position_on_grid[2])]);

				for(int k=0; k<3; k++)
					current_position[k] += edge_step_size[k];

				if(edge_step_size[0] == 0 && edge_step_size[1] == 0 && edge_step_size[2] == 0) break;
			}
			*/
		
	}
	
}



void GRIDY_SIGNED_DISTANCE_FIELD::Set_interface_btw_edges(double *edgeAA, double *edgeBB, double *edgeCC, double *delta_grid, double *bounding_box_AA, double *margin)
{
	double edge_step_size[3];
	double min_delta_grid = MIN3(delta_grid[0], delta_grid[1], delta_grid[2]);

	double edge_distance = (sqrt((edgeBB[0]-edgeAA[0])*(edgeBB[0]-edgeAA[0]) + (edgeBB[1]-edgeAA[1])*(edgeBB[1]-edgeAA[1]) + (edgeBB[2]-edgeAA[2])*(edgeBB[2]-edgeAA[2])));
	//edge_distance = 1.0f;


	double step_size = edge_distance / min_delta_grid * 5.0f;
	//printf("edge_distance : %f, min_delta_grid : %f, step size : %f\n", edge_distance, min_delta_grid, step_size);

	for(int k=0; k<3; k++)
		edge_step_size[k] = (edgeBB[k] - edgeAA[k])/step_size;
		//edge_step_size[k] = edgeBB[k] - edgeAA[k]/edge_distance;
		//edge_step_size[k] = abs(edgeBB[k] - edgeAA[k])/step_size;



	double current_position[3];
	//int current_position_on_grid[3];
		
	current_position[0] = edgeAA[0];	current_position[1] = edgeAA[1];	current_position[2] = edgeAA[2];

	//printf("%f %f / %f %f / %f %f\n", abs(current_position[0] - edgeBB[0]), abs(edge_step_size[0]), abs(current_position[1] - edgeBB[1]), abs(edge_step_size[1]), abs(current_position[2] - edgeBB[2]), abs(edge_step_size[2]));

	//while(current_position[0] <= edgeBB[0] && current_position[1] <= edgeBB[1] && current_position[2] <= edgeBB[2])
	while(abs(current_position[0] - edgeBB[0]) >= abs(edge_step_size[0]) && abs(current_position[1] - edgeBB[1]) >= abs(edge_step_size[1]) && abs(current_position[2] - edgeBB[2]) >= abs(edge_step_size[2]))
	{
		Set_interface_btw_points(current_position, edgeCC, delta_grid, bounding_box_AA, margin);

		for(int k=0; k<3; k++)
			current_position[k] += edge_step_size[k];

		if(edge_step_size[0] == 0 && edge_step_size[1] == 0 && edge_step_size[2] == 0) break;
	}
}




void GRIDY_SIGNED_DISTANCE_FIELD::Set_interface_btw_points(double *edgeAA, double *edgeBB, double *delta_grid, double *bounding_box_AA, double *margin)
{
	double edge_step_size[3];
	double min_delta_grid = MIN3(delta_grid[0], delta_grid[1], delta_grid[2]);

	double edge_distance = (sqrt((edgeBB[0]-edgeAA[0])*(edgeBB[0]-edgeAA[0]) + (edgeBB[1]-edgeAA[1])*(edgeBB[1]-edgeAA[1]) + (edgeBB[2]-edgeAA[2])*(edgeBB[2]-edgeAA[2])));
	//edge_distance = 1.0f;


	double step_size = edge_distance / min_delta_grid * 5.0f;
	//printf("edge_distance : %f, min_delta_grid : %f, step size : %f\n", edge_distance, min_delta_grid, step_size);

	for(int k=0; k<3; k++)
		edge_step_size[k] = (edgeBB[k] - edgeAA[k])/step_size;
		//edge_step_size[k] = edgeBB[k] - edgeAA[k]/edge_distance;
		//edge_step_size[k] = abs(edgeBB[k] - edgeAA[k])/step_size;



	double current_position[3];
	int current_position_on_grid[3], distance_vertex_btw_grid[3];
		
	current_position[0] = edgeAA[0];	current_position[1] = edgeAA[1];	current_position[2] = edgeAA[2];

	//printf("%f %f / %f %f / %f %f\n", abs(current_position[0] - edgeBB[0]), abs(edge_step_size[0]), abs(current_position[1] - edgeBB[1]), abs(edge_step_size[1]), abs(current_position[2] - edgeBB[2]), abs(edge_step_size[2]));

	while(abs(current_position[0] - edgeBB[0]) >= abs(edge_step_size[0]) && abs(current_position[1] - edgeBB[1]) >= abs(edge_step_size[1]) && abs(current_position[2] - edgeBB[2]) >= abs(edge_step_size[2]))
	{
		for(int k=0; k<3; k++)
		{
			current_position_on_grid[k] = (int)abs((current_position[k] - bounding_box_AA[k] - margin[k]) / delta_grid[k] + 1);
			distance_vertex_btw_grid[k] = (int)(current_position[k] - (int)current_position_on_grid[k]);
		}

		sign[IX3((int)current_position_on_grid[0], (int)current_position_on_grid[1], (int)current_position_on_grid[2])] = GRIDY_SDF_INTERFACE;
		
		/*
		double temp_phi0 = abs(phi[IX3(current_position_on_grid[0], current_position_on_grid[1], current_position_on_grid[2])]);
		double temp_phi = sqrt(distance_vertex_btw_grid[0]*distance_vertex_btw_grid[0] + distance_vertex_btw_grid[1]*distance_vertex_btw_grid[1] + distance_vertex_btw_grid[2]*distance_vertex_btw_grid[2]);

		if(temp_phi < temp_phi0) 
		{
		//	phi[IX3(current_position_on_grid[0], current_position_on_grid[1], current_position_on_grid[2])] = temp_phi;

			// normal 근사화 모듈이 필요함
			
			//double norm_vertex_normal = sqrt(vertex_normal[0]*vertex_normal[0] + vertex_normal[1]*vertex_normal[1] + vertex_normal[2]*vertex_normal[2]);

			//for(int i=0; i<3; i++)
			//	normal[IX3((int)vertex_position_on_grid[0], (int)vertex_position_on_grid[1], (int)vertex_position_on_grid[2])][i] = vertex_normal[i]/norm_vertex_normal;
			
		}
		*/

		for(int k=0; k<3; k++)
			current_position[k] += edge_step_size[k];

		if(edge_step_size[0] == 0 && edge_step_size[1] == 0 && edge_step_size[2] == 0) break;
	}

}





//void GRIDY_SIGNED_DISTANCE_FIELD::Find_interface_and_allocate_sign(short *sign, double *phi)
void GRIDY_SIGNED_DISTANCE_FIELD::Find_interface_and_allocate_sign(void)
{
	short current_sign;

	for(int i=0; i<=mWidth+1; i++)
		for(int j=0; j<=mHeight+1; j++)
		{
			current_sign = GRIDY_SDF_OUTER;

			for(int k=0; k<=mDepth+1; k++)
			{
				if(sign[IX3(i,j,k)] == GRIDY_SDF_INTERFACE) 
				{
					if(current_sign == GRIDY_SDF_OUTER) current_sign = GRIDY_SDF_INNER;
					else current_sign = GRIDY_SDF_OUTER;
				}
				else sign[IX3(i,j,k)] = current_sign;
			}
		}
}
	




//void GRIDY_SIGNED_DISTANCE_FIELD::Loop_over_the_unknown_grid_points(short *sign, double *phi)
void GRIDY_SIGNED_DISTANCE_FIELD::Loop_over_the_unknown_grid_points(void)
{



}




void GRIDY_SIGNED_DISTANCE_FIELD::Make_unsigned_distance_field(double *phi, double h)
{
	bool *flag = (bool*)malloc(sizeof(bool)*(mWidth+2)*(mHeight+2)*(mDepth+2));		
	double smallist_phi = (double)GRIDY_UNKNOWN_GRID;

	GRIDY_LIST_STACK3<int> stack;
	int mPos[3];						// temp arr for stack

	//cout << "Make unsigned distance field : " << (mWidth+2)*(mHeight+2)*(mDepth+2) << endl;

	int cnt_not_unknown = 0;

	// init arr and find smallest phi value
	for(int i=0; i<=mWidth+1; i++)
		for(int j=0; j<=mHeight+1; j++)
			for(int k=0; k<=mDepth+1; k++)
			{
				if(phi[IX3(i,j,k)] < 1.7f) flag[IX3(i,j,k)] = true;
				else flag[IX3(i,j,k)] = false;

				if(phi[IX3(i,j,k)] < smallist_phi) 
				{
					smallist_phi = phi[IX3(i,j,k)];
					mPos[0] = i;	mPos[1] = j;	mPos[2] = k;
				}
				cnt_not_unknown ++;

				/*
				//if(arr[IX3(i,j,k)] == 999999) flag[IX3(i,j,k)] = false;
				if(phi[IX3(i,j,k)] >= (double)GRIDY_UNKNOWN_GRID-1) flag[IX3(i,j,k)] = false;
				else
				{
					flag[IX3(i,j,k)] = true;

					if(phi[IX3(i,j,k)] < smallist_phi) 
					{
						smallist_phi = phi[IX3(i,j,k)];
						mPos[0] = i;	mPos[1] = j;	mPos[2] = k;
					}
					cnt_not_unknown ++;
				}
				*/
			}

	
	//printf("%d %d %d %f\n", mPos[0], mPos[1], mPos[2], phi[IX3(mPos[0], mPos[1], mPos[2])]);
	flag[IX3(mPos[0], mPos[1], mPos[2])] = true;
	stack.push(mPos[0]-1, mPos[1], mPos[2]);
	stack.push(mPos[0]+1, mPos[1], mPos[2]);
	stack.push(mPos[0], mPos[1]-1, mPos[2]);
	stack.push(mPos[0], mPos[1]+1, mPos[2]);
	stack.push(mPos[0], mPos[1], mPos[2]-1);
	stack.push(mPos[0], mPos[1], mPos[2]+1);
		
	double mDist = (double)GRIDY_UNKNOWN_GRID;

	while(stack.getSize() > 0)
	{
		mPos[2] = stack.pop();	mPos[1] = stack.pop();	mPos[0] = stack.pop();	
		mDist = phi[IX3(mPos[0], mPos[1], mPos[2])];
		double mDist_min = 999999;
		//printf("POP : %d %d %d(%f), %d left.\n", mPos[0], mPos[1], mPos[2], mDist, stack.getSize()/3);

		// check neighbor
		if(mPos[0] >= 1)
		{
			if(flag[IX3(mPos[0]-1, mPos[1], mPos[2])] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0]-1, mPos[1], mPos[2])]+h, mDist_min);
			}
			else stack.push(mPos[0]-1, mPos[1], mPos[2]);
		}
		if(mPos[0] <= mWidth)
		{
			if(flag[IX3(mPos[0]+1, mPos[1], mPos[2])] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0]+1, mPos[1], mPos[2])]+h, mDist_min);
			}
			else stack.push(mPos[0]+1, mPos[1], mPos[2]);
		}
		if(mPos[1] >= 1)
		{
			if(flag[IX3(mPos[0], mPos[1]-1, mPos[2])] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0], mPos[1]-1, mPos[2])]+h, mDist_min);
			}
			else stack.push(mPos[0], mPos[1]-1, mPos[2]);
		}
		if(mPos[1] <= mHeight)
		{
			if(flag[IX3(mPos[0], mPos[1]+1, mPos[2])] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0], mPos[1]+1, mPos[2])]+h, mDist_min);
			}
			else stack.push(mPos[0], mPos[1]+1, mPos[2]);
		}			
		if(mPos[2] >= 1)
		{
			if(flag[IX3(mPos[0], mPos[1], mPos[2]-1)] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0], mPos[1], mPos[2]-1)]+h, mDist_min);
			}
			else stack.push(mPos[0], mPos[1], mPos[2]-1);
		}
		if(mPos[2] <= mDepth)
		{
			if(flag[IX3(mPos[0], mPos[1], mPos[2]+1)] == true)
			{
				mDist_min = MIN(phi[IX3(mPos[0], mPos[1], mPos[2]+1)]+h, mDist_min);
			}
			else stack.push(mPos[0], mPos[1], mPos[2]+1);
		}		
		

		phi[IX3(mPos[0], mPos[1], mPos[2])] = mDist_min;
		flag[IX3(mPos[0], mPos[1], mPos[2])] = true;
	}

	
	free(flag);
}






void GRIDY_SIGNED_DISTANCE_FIELD::Inner_flood_fill(int x, int y, int z)
{

	GRIDY_LIST_STACK3<int> stack;
	int mPos[3];						// temp arr for stack
	mPos[0] = x;	mPos[1] = y;	mPos[2] = z;

	//std::cout << "Inner_flood_fill " << (mWidth+2)*(mHeight+2)*(mDepth+2) << endl;

	sign[IX3(mPos[0], mPos[1], mPos[2])] = GRIDY_SDF_INNER;
	stack.push(mPos[0]-1, mPos[1], mPos[2]);
	stack.push(mPos[0]+1, mPos[1], mPos[2]);
	stack.push(mPos[0], mPos[1]-1, mPos[2]);
	stack.push(mPos[0], mPos[1]+1, mPos[2]);
	stack.push(mPos[0], mPos[1], mPos[2]-1);
	stack.push(mPos[0], mPos[1], mPos[2]+1);


	
	while(stack.getSize() > 0) 
	{
		mPos[2] = stack.pop();	mPos[1] = stack.pop();	mPos[0] = stack.pop();
		//cout << "mPos[2] :" << mPos[2] << endl;
		//cout << "mPos[1] :" << mPos[1] << endl;
		//cout << "mPos[0] :" << mPos[0] << endl << endl;
		//printf("POP : %d %d %d, %d left.\n", mPos[0], mPos[1], mPos[2], stack.getSize()/3);

		// check neighbor
		if(mPos[0] >= 1)
		{
			if(sign[IX3(mPos[0]-1, mPos[1], mPos[2])] == GRIDY_SDF_OUTER) stack.push(mPos[0]-1, mPos[1], mPos[2]);
		}
		if(mPos[0] <= mWidth)
		{
			if(sign[IX3(mPos[0]+1, mPos[1], mPos[2])] == GRIDY_SDF_OUTER) stack.push(mPos[0]+1, mPos[1], mPos[2]);
		}			
		if(mPos[1] >= 1)
		{
			if(sign[IX3(mPos[0], mPos[1]-1, mPos[2])] == GRIDY_SDF_OUTER) stack.push(mPos[0], mPos[1]-1, mPos[2]);
		}
		if(mPos[1] <= mHeight)
		{
			if(sign[IX3(mPos[0], mPos[1]+1, mPos[2])] == GRIDY_SDF_OUTER) stack.push(mPos[0], mPos[1]+1, mPos[2]);
		}			
		if(mPos[2] >= 1)
		{
			if(sign[IX3(mPos[0], mPos[1], mPos[2]-1)] == GRIDY_SDF_OUTER) stack.push(mPos[0], mPos[1], mPos[2]-1);
		}
		if(mPos[2] <= mDepth)
		{
			if(sign[IX3(mPos[0], mPos[1], mPos[2]+1)] == GRIDY_SDF_OUTER) stack.push(mPos[0], mPos[1], mPos[2]+1);
		}		

		sign[IX3(mPos[0], mPos[1], mPos[2])] = GRIDY_SDF_INNER;
	}
	
}





#endif //_GRIDY_SIGNED_DISTANCE_FIELD