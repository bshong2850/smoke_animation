#include "../LIBRARY/GRIDY_SMOKE_3D.h"
#include <fstream>

class GRIDY_SMOKE_ANIMATED_MESH : public GRIDY_SMOKE_3D {
public:


protected:
	GRIDY_VECTOR_FIELD potential;
	GRIDY_SCALAR_FIELD sdf;

	float grid_cell_size = 0.2;

	int x_pos_1 = 0;
	int x_pos_2 = 0;
	int x_pos_3 = 0;

	string save_path;
	float *dens;

	void make_folder(string path, string ray_path)
	{
		char *tmp_1;
		char *tmp_2;

		tmp_1 = new char[path.size() + 1];
		tmp_2 = new char[ray_path.size() + 1];
		copy(path.begin(), path.end(), tmp_1);
		copy(ray_path.begin(), ray_path.end(), tmp_2);
		_mkdir(tmp_1);
		_mkdir(tmp_2);
	}

	void dens_init()
	{
		int dens_size = (mWidth + 2) * (mHeight + 2) * (mDepth + 2);
		dens = (float *)malloc(dens_size * sizeof(float));
		for (int i = 0; i<dens_size; i++)
		{
			dens[i] = 0.0f;
		}
	}

	void dens_free()
	{
		free(dens);
	}

	void output_data(int frame, string ray_path)
	{
		int size = mWidth * mHeight * mDepth;

		char tmp[1024];
		FILE *fp;

		save_path = ray_path;

		sprintf(tmp, "%s/%d.d", save_path.c_str(), frame);
		fp = fopen(tmp, "wb");
		fwrite(dens, sizeof(float) * size, 1, fp);
		fclose(fp);
	}



	void Compute_potential_force(GRIDY_VECTOR_FIELD &potential, GRIDY_SCALAR_FIELD &sdf, float grid_cell_size, float threshold)
	{
		mLOG.In("Compute potential force");

		float h = 1.0f / grid_cell_size;

#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++) {
					potential.u[IX3(i, j, k)] = -0.5f * h *(sdf.value[IX3(i + 1, j, k)] - sdf.value[IX3(i - 1, j, k)]);
					potential.v[IX3(i, j, k)] = -0.5f * h *(sdf.value[IX3(i, j + 1, k)] - sdf.value[IX3(i, j - 1, k)]);
					potential.w[IX3(i, j, k)] = -0.5f * h *(sdf.value[IX3(i, j, k + 1)] - sdf.value[IX3(i, j, k - 1)]);


					//	error correction 
					//	너무 튀는 값은 0으로 조정
					//if (abs(potential.u[IX3(i, j, k)]) > threshold || abs(potential.v[IX3(i, j, k)]) > threshold || abs(potential.w[IX3(i, j, k)]) > threshold) {
					//potential.u[IX3(i, j, k)] = potential.v[IX3(i, j, k)] = potential.w[IX3(i, j, k)] = 0.0f;
					//}


					//	normalize
					float norm = sqrt(potential.u[IX3(i, j, k)] * potential.u[IX3(i, j, k)] + potential.v[IX3(i, j, k)] * potential.v[IX3(i, j, k)] + potential.w[IX3(i, j, k)] * potential.w[IX3(i, j, k)]);

					if (norm > 0.0f) {
						potential.u[IX3(i, j, k)] = potential.u[IX3(i, j, k)] / norm;
						potential.v[IX3(i, j, k)] = potential.v[IX3(i, j, k)] / norm;
						potential.w[IX3(i, j, k)] = potential.w[IX3(i, j, k)] / norm;
					}

					/*potential.u[IX3(i, j, k)] = 1;
					potential.v[IX3(i, j, k)] = 0;
					potential.w[IX3(i, j, k)] = 0;*/

					//printf("%f %f %f\n", potential.u[IX3(i,j,k)], potential.v[IX3(i,j,k)], potential.w[IX3(i,j,k)] );
				}

		mLOG.Out();
	}

	void Add_potential_force(GRIDY_VECTOR_FIELD &vel, GRIDY_VECTOR_FIELD &potential, GRIDY_SCALAR_FIELD& sdf, GRIDY_SCALAR_FIELD& dens, float grid_cell_size, float coef_potential, float D_SP, float dt)
	{
		mLOG.In("Add potential force");


		//	add forces
#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++) {
					//if(abs(sdf.phi[IX3(i,j,k)]) < 2.0f * sdf.Get_grid_cell_size())

					if (sdf.value[IX3(i, j, k)] < 0.0f) {
						//float coef_dens = abs(dens.value[IX3(i, j, k)]);
						float coef_dens = 1.7f;
						vel.u[IX3(i, j, k)] -= coef_potential * coef_dens * potential.u[IX3(i, j, k)] * dt;
						vel.v[IX3(i, j, k)] -= coef_potential * coef_dens * potential.v[IX3(i, j, k)] * dt;
						vel.w[IX3(i, j, k)] -= coef_potential * coef_dens * potential.w[IX3(i, j, k)] * dt;
						/*
						vel.u[IX3(i-1,j,k)] -= coef_potential * potential.u[IX3(i,j,k)] * dt;
						vel.u[IX3(i+1,j,k)] -= coef_potential * potential.u[IX3(i,j,k)] * dt;
						vel.v[IX3(i,j-1,k)] -= coef_potential * potential.v[IX3(i,j,k)] * dt;
						vel.v[IX3(i,j+1,k)] -= coef_potential * potential.v[IX3(i,j,k)] * dt;
						vel.w[IX3(i,j,k-1)] -= coef_potential * potential.w[IX3(i,j,k)] * dt;
						vel.w[IX3(i,j,k+1)] -= coef_potential * potential.w[IX3(i,j,k)] * dt;
						*/
						//cout << coef_potential * potential.u[IX3(i,j,k)] * dt << endl;
					}
					else {
						//if(abs(sdf.phi[IX3(i,j,k)]) <= abs(sdf.Get_grid_cell_size() * mD_SP))
						//if (abs(sdf.value[IX3(i, j, k)]) <= abs(grid_cell_size * D_SP))
						{
							float coef_dens = abs(dens.value[IX3(i, j, k)]);
							//coef_dens = 1.0f;
							// 여기에 density를 넣어보면 어떨까?
							float coef_surface = 15.0f * coef_potential;

							vel.u[IX3(i, j, k)] += coef_surface * coef_dens * potential.u[IX3(i, j, k)] * dt;
							vel.v[IX3(i, j, k)] += coef_surface * coef_dens * potential.v[IX3(i, j, k)] * dt;
							vel.w[IX3(i, j, k)] += coef_surface * coef_dens * potential.w[IX3(i, j, k)] * dt;
							/*
							vel.u[IX3(i-1,j,k)] -= coef_surface * potential.u[IX3(i,j,k)] * dt;
							vel.u[IX3(i+1,j,k)] -= coef_surface * potential.u[IX3(i,j,k)] * dt;
							vel.v[IX3(i,j-1,k)] -= coef_surface * potential.v[IX3(i,j,k)] * dt;
							vel.v[IX3(i,j+1,k)] -= coef_surface * potential.v[IX3(i,j,k)] * dt;
							vel.w[IX3(i,j,k-1)] -= coef_surface * potential.w[IX3(i,j,k)] * dt;
							vel.w[IX3(i,j,k+1)] -= coef_surface * potential.w[IX3(i,j,k)] * dt;
							*/
						}
					}

					/*
					vel.u[IX3(i,j,k)] -= coef_potential * potential.u[IX3(i,j,k)] * dt;
					vel.v[IX3(i,j,k)] -= coef_potential * potential.v[IX3(i,j,k)] * dt;
					vel.w[IX3(i,j,k)] -= coef_potential * potential.w[IX3(i,j,k)] * dt;
					//printf("%f %f %f\n", potential.u[IX3(i,j,k)], potential.v[IX3(i,j,k)], potential.w[IX3(i,j,k)]);
					*/
				}


		mLOG.Out();
	}

	void LoadSDF(GRIDY_SCALAR_FIELD &sdf, std::string filename, int frame, int total_frame)
	{
		mLOG.In("Initializing SDF");

		potential.Init(mWidth, mHeight, mDepth);
		potential.Clear_data();

		sdf.Init(mWidth, mHeight, mDepth);
		sdf.Clear_data(10);

		FILE *fp = fopen(filename.c_str(), "r");

		float dist;

		int NX, NY, NZ;

		fscanf(fp, "%d %d %d ", &NX, &NY, &NZ);

		int XX, YY, ZZ;
		int grid;
		grid = (mWidth - 5 - NX / 2) / total_frame;

		//XX = (mWidth - grid) * (frame - 1) / (total_frame);  // 움직임

		if (frame > 13 && frame < 17)
		{
			XX = x_pos_1 - 8 + (grid + 12) * (frame - 13);    // 2번째
			x_pos_2 = XX;
		}
		else if (frame >= 17 && frame < 35)
		{
			XX = x_pos_2 + (grid - 3) * (frame - 17);    // 2번째
			x_pos_3 = XX;
		}
		else if (frame >= 35)
		{
			XX = x_pos_3 + (grid + 3) * 0.8 * (frame - 35);    // 2번째
		}
		else
		{
			XX = grid * (frame - 1);
			x_pos_1 = XX;
		}
		//XX = grid * (frame - 1);
		YY = mHeight / 2 - NY / 2;
		ZZ = mDepth / 2 - NZ / 2;

		for (int i = 0; i < NX; i++) {
			for (int j = 0; j < NY; j++) {
				for (int k = 0; k < NZ; k++) {
					fscanf(fp, "%e", &dist);

					sdf.value[IX3(XX + i, YY + j, ZZ + k)] = dist;
				}
			}
		}
		fclose(fp);

		Compute_potential_force(potential, sdf, grid_cell_size, 1.0);
		mLOG.Out();
	}

public:
	void obj_to_sdf()
	{
		int x = 78;
		int y = 218;
		int z = 78;

		int start_frame = 1;
		int end_frame = 16;
		int total_frame;
		total_frame = end_frame - start_frame + 1;

		for (int frame = start_frame; frame <= end_frame; frame++)
		{
			int tmp = total_frame;

			int total = (x + 2) * (y + 2) * (z + 2);

			int save_start_frame = frame - start_frame + 1;
			string x_ = to_string(x + 2);
			string y_ = to_string(y + 2);
			string z_ = to_string(z + 2);

			string total_frame_;

			string load_a, load_b, load_full, load_frame;
			string save_a, save_b, save_c, save_folder_a, save_folder, save_full, save_frame;

			load_frame = to_string(frame);
			save_frame = to_string(save_start_frame);
			total_frame_ = to_string(total_frame);

			load_a = "D:\\KU\\IITP\\smoke_animation\\frame_body_edit\\frame_";
			load_b = "_edit.obj";

			save_folder_a = "D:\\KU\\IITP\\smoke_animation\\sdf_frame_body";
			save_a = "\\body_frame_";
			save_b = "_x" + x_ + "_y" + y_ + "_z" + z_;
			save_c = ".sdf";

			save_folder = save_folder_a + total_frame_ + save_b;
			char *folder_dir = new char[save_folder.size() + 1];
			copy(save_folder.begin(), save_folder.end(), folder_dir);
			_mkdir(folder_dir);

			GRIDY_SIGNED_DISTANCE_FIELD test;
			load_full = load_a + load_frame + load_b;

			test.Make_SDF_from_obj(load_full, x, y, z, 1.0f, 0, 0, 0, 1);

			save_full = save_folder + save_a + save_frame + save_b + save_c;

			ofstream writeFile(save_full.data());
			if (writeFile.is_open())
			{
				writeFile << x + 2 << " " << y + 2 << " " << z + 2 << endl;

				for (int i = 0; i < total; i++)
					writeFile << test.sign[i] << endl;
				writeFile.close();
			}
		}
	}


	void Simulate(int width, int height, int depth, std::string path, int start_frame = 1, int end_frame = 200, int end_total_frame = 300)
	{
		obj_to_sdf();

		GRIDY_SMOKE_3D::Init(width, height, depth, path);

		/////////////////////////////////
		dens_init();
		/////////////////////////////////
		mStart_frame = start_frame;
		mEnd_frame = end_frame;
		mTotal_frame = end_total_frame;


		int sdf_frame = 1;

		string prefix_frame = "FRAME #";
		string str_frame, file_frame;

		string size, a, b, full;

		size = "_x80_y220_z80";
		a = "D:\\KU\\IITP\\SDF_file\\body_frame40" + size + "\\body_frame_";
		b = size + ".sdf";

		GRIDY_SIGNED_DISTANCE_FIELD mBunny, mAngel;

		mLOG.In("SIMULATE");

		mLOG.In("Write Info File");
		mFILE_IO.Write_simulation_info(mSimulation_path, mWidth, mHeight, mDepth, mStart_frame, mEnd_frame, mData_file_type); // mEnd_frame 이냐 mTotal_frame 이냐
		mLOG.Out();

		float mSource_x = (mWidth / 2.0f);

		for (int frame = mStart_frame; frame <= mTotal_frame; frame++) {
			str_frame = to_string(frame);
			file_frame = to_string(sdf_frame);
			mLOG.In(prefix_frame + str_frame);
			full = a + file_frame + b;
			cout << "sdf_frame : " << sdf_frame << endl;
			cout << "mEnd_frame : " << mEnd_frame << endl;

			LoadSDF(sdf, full, sdf_frame, mEnd_frame);
#pragma omp parallel for
			for (int i = 1; i <= mWidth; i++)
				for (int j = 1; j <= mHeight; j++)
					for (int k = 1; k <= mDepth; k++) {
						density0.value[IX3(i, j, k)] = (float)0.0f;
					}
			Clear_previous_step_data();
			float xx;
			float yy;
			float zz;

			int vel_point = 0;

			xx = mWidth / 2;
			yy = mHeight / 3;
			zz = mDepth / 2;

			mSource.Init(xx, yy, zz, 3.0);

#pragma omp parallel for
			for (int i = 1; i <= mWidth; i++)
				for (int j = 1; j <= mHeight; j++)
					for (int k = 1; k <= mDepth; k++) {
						if (sdf.value[IX3(i, j, k)] == 0) {
							density0.value[IX3(i, j, k)] = (float)4.0f;
						}
						else if (sdf.value[IX3(i, j, k)] == 1) {
							density0.value[IX3(i, j, k)] = (float)4.0f;
						}
						else{
							density0.value[IX3(i, j, k)] = 0.0f;
						}


					}

#pragma omp parallel for
			for (int i = 1; i <= mWidth; i++)
				for (int j = 1; j <= mHeight; j++)
					for (int k = 1; k <= mDepth; k++)
						if (sdf.value[IX3(i, j, k)] == 0)
						{
							velocity.u[IX3(i, j, k)] = (float)-0.1f;
						}

			//Add_potential_force(velocity0, potential, sdf, density, grid_cell_size, 1.0, 1.0, mDt);

			//Velocity_step();

			Density_step_for_ani(density, density0, velocity, 0.0f, mDt);

			Write_data(frame);

			sdf.~GRIDY_SCALAR_FIELD();

			if (sdf_frame >= mEnd_frame)
			{
				sdf_frame = 1;
			}
			else
			{
				sdf_frame++;
			}
			mLOG.Out();
			potential.~GRIDY_VECTOR_FIELD();
		}

		mLOG.Out();
	}

	void Velocity_step()
	{
		mLOG.In("Velocity step");



		Add_velocity(velocity, velocity0, mDt);		// Add velocity0 to velocity

		//Add_bouyance(velocity, density, 0.04f, mDt);

		SWAP_value(velocity, velocity0);

		//	Diffuse velocity (use only if we need)
		mLOG.In("Diffuse velocity");
		Diffuse(1, velocity.u, velocity0.u, mDt);
		Diffuse(2, velocity.v, velocity0.v, mDt);
		Diffuse(3, velocity.w, velocity0.w, mDt);
		mLOG.Out();



		//Add_normal_force(velocity0, 5.0f);
		//Add_vorticity_confinement(velocity, 1.6f, mDt);

		Project(velocity, velocity0.u, velocity0.v);
		SWAP_value(velocity, velocity0);


		Advect_velocity_for_ani(velocity, velocity0, mDt);
		Project(velocity, velocity0.u, velocity0.v);

		mLOG.Out();
	}
};
