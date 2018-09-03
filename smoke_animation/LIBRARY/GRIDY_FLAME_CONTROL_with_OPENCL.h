#ifndef _FLAME_CONTROL_WITH_OPENCL_H_
#define _FLAME_CONTROL_WITH_OPENCL_H_


#include "../../LIBRARY/GRIDY_FIRE_STANDARD.h"
#include "../../LIBRARY/GRIDY_OPENCL_SOLVER.h"


class GRIDY_FLAME_CONTROL : public GRIDY_FIRE
{
protected:
	int mNo_src_spheres;
	GRIDY_VECTOR_FIELD potential;					// target-shape을 위한 포텐셜 필드
	GRIDY_SCALAR_FIELD soot, soot0;					// 재(soot)을 위한 스칼라 필드

	double mD_TP, mD_SP;								// Turbulence preserving region의 phi 값, Shape preserving region의 phi 값

	// fuel particle 관련 parameter
	int mEFP_seed;									// explosive fuel particle을 얼마 간격으로 seeding 할 것인가?
	double mEFP_growing_radius;						// gr는 timestep마다 성장하는 간격
	double mEFP_max_radius;							// max는 최대 반지름
	double mEFP_init_radius;							// 초기 반지름 크기
	double mEFP_activate_threshold;					// explosive fuel particle이 activate되는 온도 값

	bool mEFP_activate_all;							// 모든 explosive fuel particle을 activate한다

	bool mMakeSoot;									// soot을 만들 것인가
	bool mForComparisionScene;						// fuel particle을 사용하지 않고 기존 기법으로 결과를 만듦

private:
	string mTarget_shape_sdf_file;
	int mTargeting;
	int mFrame;

	GRIDY_OPENCL_SOLVER mCL;

	//void Transfer_buffers_to_GPU(void);
	//void Transfer_buffers_to_CPU(void);


public:
	GRIDY_FLAME_CONTROL(void);
	~GRIDY_FLAME_CONTROL(void);

	void Init(int width, int height, int depth, string path);
	void Setting_for_shape_control(void);

	//void Init_opencl(void);

	void Simulate(int start_frame, int end_frame);
	//void Velocity_step(GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres);
	void Velocity_step(GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, GRIDY_SCALAR_FIELD& lv_spheres);
	void Density_step(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& fuel0, GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double diff, double dt);

	//void Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p);
	void Set_boundary(int b, double *x);

	void Advect_fuel(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt);
	//void Advect_temperature(GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_SCALAR_FIELD& fuel, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt);
	void Advect_temperature(GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& soot, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt);
	void Advect_soot(GRIDY_SCALAR_FIELD& soot, GRIDY_SCALAR_FIELD& soot0, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt);
	void Write_data(int frame);

	void Add_bouyancy(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &temperature, GRIDY_SIGNED_DISTANCE_FIELD& sdf, double coef_bouyance, double dt);
	void Add_source_normal_force(GRIDY_IMPLICIT_OBJECT_SPHERES &spheres, GRIDY_SCALAR_FIELD& lv_spheres, GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &dens, double coef_normal_velocity, double dt);
	void Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, GRIDY_SIGNED_DISTANCE_FIELD &sdf, double coef_vorticity, double dt);
	//void Add_potential_force(GRIDY_VECTOR_FIELD &vel, GRIDY_VECTOR_FIELD &potential, GRIDY_SIGNED_DISTANCE_FIELD& sdf, double coef_potential, double dt);
	void Add_potential_force(GRIDY_VECTOR_FIELD &vel, GRIDY_VECTOR_FIELD &potential, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_SCALAR_FIELD& dens, double coef_potential, double dt);
	void Compute_potential_force(GRIDY_VECTOR_FIELD &potential, GRIDY_SIGNED_DISTANCE_FIELD &sdf, double threshold);
};





GRIDY_FLAME_CONTROL::GRIDY_FLAME_CONTROL(void)
{
	mStart_frame = mEnd_frame = 0;
	mWidth = mHeight = mDepth = 0;

	mSimulation_path = "C:\\";
}




GRIDY_FLAME_CONTROL::~GRIDY_FLAME_CONTROL(void)
{
	// terminate this program
}





void GRIDY_FLAME_CONTROL::Setting_for_shape_control(void)
{
	// fuel control을 위한 parameter setting
	Set_parameters(0.1f, 0.0f, 0.0f);	// GRIDY 기본 parameter set

	//	시뮬레이션 결과 데이터 파일 타입
	mData_file_type = mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY;

	mMakeSoot = false;
	mForComparisionScene = false;
	//mForComparisionScene = true;
	
	//mEFP_activate_all = false;
	mEFP_activate_all = true;

	// Tetra (128) #2
	mEFP_seed = (int)(mWidth/6);
	mEFP_init_radius = 1.0f;
	//mEFP_growing_radius = mWidth/(256.0f);
	mEFP_growing_radius = mWidth / (256.0f * 2.0f);	// 염2
	mEFP_max_radius = mWidth / 32.0f;

	mEFP_activate_threshold = 0.1f;

	/*
	//mTarget_shape_sdf_file = "E:\\GRIDY\\MODEL\\SDF\\cgi (243 120 66) 20b.sdf";
	mEFP_seed = (int)(mWidth / 40);		// cgi 캐잘됨
	mD_SP = (double)mWidth/12.0f;			// 12 cells
	mD_TP = (double)mWidth/50.0f;			// 8 cells 40인가 50인가
	*/

	//mTarget_shape_sdf_file = "E:\\GRIDY\\MODEL\\SDF\\cfire2 (210 209 83) 30b.sdf";
	//mTarget_shape_sdf_file = "E:\\GRIDY\\MODEL\\SDF\\cfire_b (196 187 77) 30b.sdf";
	mTarget_shape_sdf_file = "E:\\GRIDY\\MODEL\\SDF\\taehyeong2 (320 131 53) 20b.sdf";

	//mEFP_seed = (int)(mWidth / 40);			// 염
	//mD_SP = (double)mWidth / 10.0f;			// 12 cells
	//mD_TP = (double)mWidth / 105.0f;			// 8 cells
	
	//mEFP_seed = (int)(mWidth / 35);			// 염2
	//mD_SP = (double)mWidth / 10.0f;			// 12 cells
	//mD_TP = (double)mWidth / 49.0f;			// 8 cells
	//mEFP_seed = (int)(mWidth / 1);		

	mEFP_seed = (int)(mWidth / 80);			// 김태형2
	mD_SP = (double)mWidth / 10.0f;			// 12 cells
	mD_TP = (double)mWidth / 320.0f;			// 8 cells
	
}






void GRIDY_FLAME_CONTROL::Simulate(int start_frame, int end_frame)
{
	mLOG.Set_file_out(mSimulation_path, "GRIDY_FLAME_CONTROL", 4, 2);

	mLOG.In("Init Simulation");

	mStart_frame = start_frame;
	mEnd_frame = end_frame;

	string prefix_frame = "FRAME #";
	char str_frame[10] = {};


	print_platforms_devices();
	Init_opencl();

	string mTarget_shape_file;



	//GRIDY_SIGNED_DISTANCE_FIELD mTarget_shape;
	GRIDY_SIGNED_DISTANCE_FIELD mTarget_shape;	// SDF of target shape
	GRIDY_IMPLICIT_OBJECT_SPHERES mFuel_particle((int)(mWidth / mEFP_seed * mHeight / mEFP_seed * mDepth / mEFP_seed), (double)mEFP_max_radius);
	//GRIDY_IMPLICIT_OBJECT_SPHERES mFuel_particle;
	GRIDY_SCALAR_FIELD mFuel_particle_index(mWidth, mHeight, mDepth);

	mLOG.Out();


	
	// 귀는 64, 42, 24에서 잘 나왔음...

	mLOG.In("SIMULATE");


	mLOG.In("Write Info File");
	mFILE_IO.Write_simulation_info(mSimulation_path, mWidth, mHeight, mDepth, mStart_frame, mEnd_frame, mData_file_type);
	mLOG.Out();

	mLOG.In("Input Object from SDF file");
	mTarget_shape.Input_SDF(mTarget_shape_sdf_file, mWidth, mHeight, mDepth);
	mLOG.Out();

	mLOG.In("Compute potential force from SDF");
	Compute_potential_force(potential, mTarget_shape, 1.0f);
	mLOG.Out();

	
	mLOG.In("Seed Fuel particles");
	int cnt_fp = 0;	// cnt fuel particles

	// 여기는 병렬 처리하면 안됨. 카운터로 배열 위치를 지정하기 때문.
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				/*
				// 주의 주의 주의!!
				if(k > mDepth * 3.0f / 4.0f) mEFP_seed = (int)(mWidth / 20);
				else mEFP_seed = (int)(mWidth / 6);
				// 주의 주의 주의!!
				*/
				if ((i % mEFP_seed == 0) && (j % mEFP_seed == 0) && (k % mEFP_seed == 0) && (mTarget_shape.phi[IX3(i, j, k)] < -1.0f * (mTarget_shape.Get_grid_cell_size() * (mD_TP + mEFP_max_radius / 5.0f))))
				{
					double temp_rnd_x = (double)((rand() % (int)(mEFP_seed * 2 / 6.0f)) - mEFP_seed / 6.0f);
					double temp_rnd_y = (double)((rand() % (int)(mEFP_seed * 2 / 6.0f)) - mEFP_seed / 6.0f);
					double temp_rnd_z = (double)((rand() % (int)(mEFP_seed * 2 / 6.0f)) - mEFP_seed / 6.0f);

					int temp_x = (int)(i + temp_rnd_x);
					int temp_y = (int)(j + temp_rnd_y);
					int temp_z = (int)(k + temp_rnd_z);

					if (mTarget_shape.phi[IX3(temp_x, temp_y, temp_z)] < -1.0f * (mTarget_shape.Get_grid_cell_size() * (mD_TP + mEFP_max_radius / 5.0f)))
						//if(mTarget_shape.phi[IX3(temp_x, temp_y, temp_z)] < -1.0f * (mTarget_shape.Get_grid_cell_size() * (mD_TP + mFP_max_radius * 2.0f) ))
					{
						mFuel_particle.Set(cnt_fp, (double)(i + temp_rnd_x), (double)(j + temp_rnd_y), (double)(k + temp_rnd_z), mEFP_init_radius);
					}
					else
					{
						mFuel_particle.Set(cnt_fp, (double)i, (double)j, (double)k, mEFP_init_radius);
					}

					printf("cnt_fuel_particles #%d : %f %f %f\n", cnt_fp, mFuel_particle.mCx[cnt_fp], mFuel_particle.mCy[cnt_fp], mFuel_particle.mCz[cnt_fp]);

					if (mEFP_activate_all == true)
					{
						mFuel_particle.isActivate[cnt_fp] = true;
						mFuel_particle.isGrowing[cnt_fp] = true;
					}
					else
					{
						mFuel_particle.isActivate[cnt_fp] = false;
						mFuel_particle.isGrowing[cnt_fp] = false;
					}
					cnt_fp++;
				}
			}

	mFuel_particle.mNo = cnt_fp;

	cout << "seeded particles # : " << cnt_fp << endl;
	//cout << (mWidth / mEFP_seed) * (mHeight / mEFP_seed) * (mDepth / mEFP_seed) << endl;

	mLOG.Out();
	


	if (mEFP_activate_all == false)
	{
		mLOG.In("Select initial-seeded fuel particle");


		// center (for tetra)
		double initX = mWidth / 2.0f;
		double initY = mHeight / 2.0f;
		double initZ = mDepth / 2.0f;



		int min_idx = -1;
		double min_value = 99999999.0f;

		for (int l = 0; l<mFuel_particle.mNo; l++)
		{
			double t_value;
			t_value = sqrt((mFuel_particle.mCx[l] - initX) * (mFuel_particle.mCx[l] - initX) + (mFuel_particle.mCy[l] - initY) * (mFuel_particle.mCy[l] - initY) + (mFuel_particle.mCz[l] - initZ) * (mFuel_particle.mCz[l] - initZ));

			//cout << t_value << endl;
			if (t_value < min_value)
			{
				min_value = t_value;
				min_idx = l;
			}
		}

		mFuel_particle.isActivate[min_idx] = true;
		mFuel_particle.isGrowing[min_idx] = true;

		cout << "initial-seed particle : " << min_idx << endl;

		mLOG.Out();
	}


	
	//for (int i = 0; i < mFuel_particle.mNo; i++)
	//		mFuel_particle.isActivate[i] = false;
	
	/*
	float mSphere_x = 210.0f;
	float mSphere_y = mHeight;

	float mSphere_x0 = mSphere_x;
	float mSphere_y0 = mSphere_y;

	double mSphere_vx = 0.0f;
	double mSphere_vy = 0.0f;

	double mSphere_ax = 0.0f;
	double mSphere_ay = -9.8f * 0.2f;

	double mSphere_stiff = 0.7f;

	mSrc_sphere.Init(mWidth / 2.0f, mHeight / 4.0f, mDepth / 2.0f, 5.0f);
	*/

	for (int frame = mStart_frame; frame <= mEnd_frame; frame++)
	{
		_itoa_s(frame, str_frame, 10);
		mLOG.In(prefix_frame + str_frame);

		mFrame = frame;

		//Clear_previous_step_data();


		//if(frame >= 150)	mTargeting = 1;
		mTargeting = 1;
		//mTargeting = 0;
		//if(frame > 40)	mTargeting = 2;

		//mSrc_sphere.Init(mWidth/2.0f, mHeight/2.0f, mDepth/2.0f, mWidth/10.0f);

		/*
		if(mTargeting > 0) mFP_max_radius = mWidth/10.0f;
		else mFP_max_radius = mWidth/15.0f;
		*/

		mFuel_particle.mMax_radius = mEFP_max_radius;

		mFuel_particle_index.Clear_data(-1.0f);

		/*
		// turn off activate
		if (frame > 150)
		{
			for (int i = 0; i < mFuel_particle.mNo; i++)
				mFuel_particle.isActivate[i] = false;
		}
		*/

		/*
		// CGI 한 글자씩
		for (int i = 0; i < mFuel_particle.mNo; i++)
		{
			mFuel_particle.isActivate[i] = false;

			if (frame < 30)
			{
				if (mFuel_particle.mCx[i] < 100) mFuel_particle.isActivate[i] = true;
			}
			else if (frame < 60)
			{
				if (mFuel_particle.mCx[i] < 200) mFuel_particle.isActivate[i] = true;
			}
			else if (frame < 90)
			{
				mFuel_particle.isActivate[i] = true;
			}
			else if (frame > 170 && frame <= 200)
			{
				if (mFuel_particle.mCx[i] < 100) mFuel_particle.isActivate[i] = false;
				else mFuel_particle.isActivate[i] = true;
			}
			else if (frame > 200 && frame <= 230)
			{
				if (mFuel_particle.mCx[i] < 200) mFuel_particle.isActivate[i] = false;
				else mFuel_particle.isActivate[i] = true;
			}
			else if (frame > 230)
			{
				mFuel_particle.isActivate[i] = false;
			}
			else 
				mFuel_particle.isActivate[i] = true;
		}
		*/

		// dissipate
		/*
		if (frame > 100)
		{
			for (int i = 0; i < mFuel_particle.mNo; i++)
			{
				//if (mFuel_particle.mCy[i] < (frame + 10) && mFuel_particle.mCy[i] > frame) mFuel_particle.isActivate[i] = true;
				//else mFuel_particle.isActivate[i] = false;

				if (mFuel_particle.mCy[i] < (frame - 100))// && mFuel_particle.mCy[i] > (frame/2 - 45)) 
				{
					mFuel_particle.isActivate[i] = false;
					mFuel_particle.mRadius[i] = 1.0f;
				}
				//else mFuel_particle.isActivate[i] = true;
			}			
		}
		*/



		/*
		// Moving sphere
		mSphere_x0 = mSphere_x;
		mSphere_y0 = mSphere_y;

		mSphere_vx += mSphere_ax * mDt;
		mSphere_vy += mSphere_ay * mDt;

		mSphere_x += mSphere_vx * mDt;
		mSphere_y += mSphere_vy * mDt;


		// boundary check
		if (mSphere_y < mSrc_sphere.radius)
		{
		mSphere_vy *= -1.0f * mSphere_stiff;
		mSphere_vx *= mSphere_stiff;
		}
		if ( mSphere_x > (mWidth - mSrc_sphere.radius) || mSphere_x < (mSrc_sphere.radius) )
		{
		mSphere_vx *= -1.0f * mSphere_stiff;
		mSphere_vy *= mSphere_stiff;
		}

		mSphere_x = mSphere_x0 + mSphere_vx * mDt;
		mSphere_y = mSphere_y0 + mSphere_vy * mDt;

		mSrc_sphere.Init(mSphere_x, mSphere_y, mDepth / 2.0f, mSrc_sphere.radius);	// update sphere position
		*/





		mLOG.In("Calculate Fuel particle area");
		#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++)
				{
					for (int l = 0; l<mFuel_particle.mNo; l++)
					{
						if (mFuel_particle.Levelset(l, i, j, k) <= 0)
						{
							if (mFuel_particle.isActivate[l] == true)
							{
								int current_cx = (int)mFuel_particle.mCx[l];
								int current_cy = (int)mFuel_particle.mCy[l];
								int current_cz = (int)mFuel_particle.mCz[l];

								if (mTarget_shape.phi[IX3(i, j, k)] < -1.0f * abs(mTarget_shape.Get_grid_cell_size() * mD_TP)) // fuel propagation region
								{
									mFuel_particle_index.value[IX3(i, j, k)] = (double)l;
									break;
									//l = mSrc_spheres.mNo;	// for exit loop
								}
								else
								{
									//mLv_spheres.value[IX3(i,j,k)] = l;
									mFuel_particle.isGrowing[l] = false;
								}
							}
							else if (mEFP_activate_all == false && mTargeting > 0)
							{
								if (temperature0.value[IX3(i, j, k)] >= mEFP_activate_threshold)
								{
									mFuel_particle.isActivate[l] = true;
									mFuel_particle.isGrowing[l] = true;
									//cout << "ACTIVATE! : " << density0.value[IX3(i,j,k)] << endl;
									printf("ACTIVATE : #%d %f(%d, %d, %d)\n", l, density0.value[IX3(i, j, k)], i, j, k, mFuel_particle);
								}
							}
						}
					}
				}

		mLOG.Out();



		Clear_previous_step_data();


		/*
		#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++)
				{
					if (mSrc_sphere.Levelset(i, j, k) <= 0) density0.value[IX3(i, j, k)] += 0.5f;
				}
		*/

		double reduce_coef = 1.0f;
		//if (mFrame > 50) reduce_coef = (double)((80 - mFrame) / (80 - 50));
		//if (mFrame > 50) reduce_coef = 0.0f;
		//if (reduce_coef < 0.0f) reduce_coef = 0.0f;
		//reduce_coef = powf(reduce_coef, 3);


		#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++)
				{
					if (mForComparisionScene == true)
					{
						if (mTarget_shape.phi[IX3(i, j, k)] < 0.0f)				// Check Rapid consumption area
						{
							density0.value[IX3(i, j, k)] = 1.0f;
						}
					}
					else
					{
						if (mFuel_particle_index.value[IX3(i, j, k)] >= 0.0f)
						{
							//cout << "index : " << mFuel_particle_index.value[IX3(i,j,k)] << endl;

							//density0.value[IX3(i,j,k)] += mWidth/240.0f;
							
							//density0.value[IX3(i, j, k)] += reduce_coef * (mWidth / 120.0f);	// for bunny
							density0.value[IX3(i, j, k)] += (mWidth / 120.0f);	// for bunny
							
							
							//density0.value[IX3(i,j,k)] = 1.0f;
							if (density0.value[IX3(i, j, k)] > 1.0f) density0.value[IX3(i, j, k)] = 1.0f;
						}
					}
				}

		Velocity_step(mTarget_shape, mFuel_particle, mFuel_particle_index);

		Density_step(density, density0, temperature, temperature0, velocity, mTarget_shape, mFuel_particle, mDiffuse, mDt);

		Write_data(frame);


		mLOG.In("Fuel particle growing");
		printf("Growing : ");
		for (int l = 0; l<mFuel_particle.mNo; l++)
		{

			//mSrc_spheres.Growing_radius(l, mSrc_spheres.mMax_radius / 5.0f);
			if (mFuel_particle.isGrowing[l] == true)
			{
				mFuel_particle.Growing_radius(l, mEFP_growing_radius);
				printf("%d ", l);
			}
		}
		printf("\n");
		mLOG.Out();

		mLOG.Out();
	}

	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Init(int width, int height, int depth, string path = "C:\\")
{
	mLOG.In("INITIALIZE SIMULATION");

	mWidth = width;
	mHeight = height;
	mDepth = depth;

	mSimulation_path = path;
	mTargeting = 0;

	density.Init(mWidth, mHeight, mDepth);
	density0.Init(mWidth, mHeight, mDepth);

	temperature.Init(mWidth, mHeight, mDepth);
	temperature0.Init(mWidth, mHeight, mDepth);

	velocity.Init(mWidth, mHeight, mDepth);
	velocity0.Init(mWidth, mHeight, mDepth);

	potential.Init(mWidth, mHeight, mDepth);

	//soot.Init(mWidth, mHeight, mDepth);
	//soot0.Init(mWidth, mHeight, mDepth);

	mLOG.In("Setting for shape control");
	Setting_for_shape_control();
	mLOG.Out();

	mLOG.Out();
}







void GRIDY_FLAME_CONTROL::Velocity_step(GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, GRIDY_SCALAR_FIELD& lv_spheres)
{
	mLOG.In("Velocity step w/ Fuel-control");




	// 0.1 0.02 0.05
	// 20b에서는 0.2, 0.05, 0.02
	Add_velocity(velocity, velocity0, mDt);										// Add velocity0 to velocity
	////Add_source_normal_force(src_spheres, velocity, 3.0f, mDt);				// originally 5.0f, 128 bunny 10.0f;
	////Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 5.00f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f; 

	/*
	// CGI
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.1f * (64.0f / (double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f; 
	if (mTargeting > 0)
	{
		Add_potential_force(velocity, potential, sdf, density, 0.03f*(64.0f / (double)mWidth), mDt); // 0.08?
	}
	Add_bouyancy(velocity, temperature, sdf, 0.01f*(64.0f / (double)mWidth), mDt);	// 옆구리 뾰죽이 0.3? 0.15는 안정적
	*/

	/*
	// Tetra (128) #2
	//Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.2f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f;
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.15f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f;
	if(mTargeting > 0)
	{
	//Add_potential_force(velocity, potential, sdf, density, 0.13f*(64.0f/(double)mWidth), mDt);
	Add_potential_force(velocity, potential, sdf, density, 0.10f*(64.0f/(double)mWidth), mDt);
	}
	////Add_potential_force(velocity, potential, sdf, 0.5f, mDt);				// tetra(60)
	Add_bouyancy(velocity, temperature, sdf, 0.01f*(64.0f/(double)mWidth), mDt);	// 옆구리 뾰죽이 0.3? 0.15는 안정적
	*/

	/*
	// Teddy Bear (128)
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.4f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f;
	if(mTargeting > 0)
	{
	Add_potential_force(velocity, potential, sdf, density, 0.1f*(64.0f/(double)mWidth), mDt); // 0.08?
	}
	Add_bouyancy(velocity, temperature, sdf, 0.02f*(64.0f/(double)mWidth), mDt);	// 옆구리 뾰죽이 0.3? 0.15는 안정적
	*/

	/*
	// bunny (128)
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.5f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f;
	if(mTargeting > 0)
	{
	Add_potential_force(velocity, potential, sdf, density, 0.08f*(64.0f/(double)mWidth), mDt); // 0.08?
	}
	Add_bouyancy(velocity, temperature, sdf, 0.02f*(64.0f/(double)mWidth), mDt);	// 옆구리 뾰죽이 0.3? 0.15는 안정적
	*/

	/*
	// cgi 243
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.4f * (64.0f / (double)mWidth), mDt);
	if (mTargeting > 0)
	{
		Add_potential_force(velocity, potential, sdf, density, 0.03f*(64.0f / (double)mWidth), mDt);
	}
	Add_bouyancy(velocity, temperature, sdf, 0.03f*(64.0f / (double)mWidth), mDt);
	*/


	
	// 염
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.7f * (64.0f / (double)mWidth), mDt);
	if (mTargeting > 0)
	{
		Add_potential_force(velocity, potential, sdf, density, 0.03f*(64.0f / (double)mWidth), mDt);
	}
	Add_bouyancy(velocity, temperature, sdf, 0.03f*(64.0f / (double)mWidth), mDt);
	
	/*
	// dragon
	Add_source_normal_force(src_spheres, lv_spheres, velocity, density, 0.5f * (64.0f/(double)mWidth), mDt);		// originally 5.0f, 128 bunny 10.0f;
	if(mTargeting > 0)
	{
	Add_potential_force(velocity, potential, sdf, density, 0.1f*(64.0f/(double)mWidth), mDt); // 0.08?
	}
	Add_bouyancy(velocity, temperature, sdf, 0.03f*(64.0f/(double)mWidth), mDt);	// 옆구리 뾰죽이 0.3? 0.15는 안정적
	*/
	/*
	//	비교 장면을 위한 구문
	mLOG.In("Add potential force for comparision");

	double coef_potential = 0.075f;

	#pragma omp parallel for
	for(int i=1; i<=mWidth; i++)
	for(int j=1; j<=mHeight; j++)
	for(int k=1; k<=mDepth; k++)
	{
	//if(abs(sdf.phi[IX3(i,j,k)]) < 2.0f * sdf.Get_grid_cell_size())
	if(sdf.phi[IX3(i,j,k)] <= 0.0f)
	{
	velocity.u[IX3(i-1,j,k)] -= coef_potential * potential.u[IX3(i,j,k)] * mDt;
	velocity.u[IX3(i+1,j,k)] -= coef_potential * potential.u[IX3(i,j,k)] * mDt;
	velocity.v[IX3(i,j-1,k)] -= coef_potential * potential.v[IX3(i,j,k)] * mDt;
	velocity.v[IX3(i,j+1,k)] -= coef_potential * potential.v[IX3(i,j,k)] * mDt;
	velocity.w[IX3(i,j,k-1)] -= coef_potential * potential.w[IX3(i,j,k)] * mDt;
	velocity.w[IX3(i,j,k+1)] -= coef_potential * potential.w[IX3(i,j,k)] * mDt;

	//cout << coef_potential * potential.u[IX3(i,j,k)] * dt << endl;
	}
	}
	mLOG.Out();

	//Add_bouyancy(velocity, temperature, sdf, (0.15f-coef_potential)*(64.0f/(double)mWidth), mDt);	// ForComparision small flame
	Add_bouyancy(velocity, temperature, sdf, (0.30f-coef_potential)*(64.0f/(double)mWidth), mDt);	// ForComparision small flame
	*/


	// 왼쪽 위로 올라가는 현상?


	// 3.0f, 0.5f, 0.3f
	// 속도는 1.5 0.25 0.15가 가장 나은 듯
	// 귀는 1.5 0.25 0.15에서 나왔음
	// 보장된 버젼은 2.0, 0.25, 0.15



	/*
	SWAP_value(velocity, velocity0);
	//	Diffuse velocity (use only if we need)
	mLOG.In("Diffuse velocity");
	Diffuse(1, velocity.u, velocity0.u, mDt);
	Diffuse(2, velocity.v, velocity0.v, mDt);
	Diffuse(3, velocity.w, velocity0.w, mDt);
	mLOG.Out();
	SWAP_value(velocity, velocity0);
	*/



	////Add_vorticity_confinement(velocity, 5.0f, mDt); 
	if (mTargeting > 0) Add_vorticity_confinement(velocity, sdf, 20.0f*(64.0f / (double)mWidth), mDt);	// 3이 좋았는데 5도 버티는 듯
	else Add_vorticity_confinement(velocity, sdf, 20.0f*(64.0f / (double)mWidth), mDt);


	Project_using_GPU(velocity, velocity0.u, velocity0.v);
	SWAP_value(velocity, velocity0);

	Advect_velocity(velocity, velocity0, mDt);
	Project_using_GPU(velocity, velocity0.u, velocity0.v);

	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Density_step(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& fuel0, GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double diff, double dt)
{
	mLOG.In("Density step w/ Fuel-control");

	Add_source(fuel, fuel0, dt);
	SWAP_value(fuel0, fuel);
	//Diffuse(0, fuel.value, fuel0.value, dt);
	//SWAP_value(fuel0, fuel);
	Advect_fuel(fuel, fuel0, vel, sdf, src_spheres, dt);



	Add_source(temperature, temperature0, dt);
	SWAP_value(temperature0, temperature);
	//Diffuse(0, temperature.value, temperature0.value, dt);
	//SWAP_value(temperature0, temperature);
	Advect_temperature(temperature, temperature0, fuel, soot, vel, sdf, src_spheres, dt);



	if (mMakeSoot == true)
	{
		Add_source(soot, soot0, dt);
		SWAP_value(soot0, soot);
		//Diffuse(0, soot.value, soot0.value, dt);
		//SWAP_value(soot0, soot);
		Advect_fuel(soot, soot0, vel, sdf, src_spheres, dt);
	}

	mLOG.Out();
}





//	Add bouyance must be execute before projection
void GRIDY_FLAME_CONTROL::Add_bouyancy(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &temperature, GRIDY_SIGNED_DISTANCE_FIELD &sdf, double coef_bouyance, double dt)
{
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				// is inside?
				double coef_sign = 1.0f;

				double min_threshold = 0.1f;
				double temp_dens = temperature.value[IX3(i, j, k)];

				if (temp_dens < min_threshold && temp_dens > 0.0f) temp_dens = min_threshold;

				if (mTargeting > 0)
				{
					if (abs(sdf.phi[IX3(i, j, k)]) <= sdf.Get_grid_cell_size() * 2.0f)	// 바운더리 근처라면!
					{
						vel.v[IX3(i, j - 1, k)] += 0.5f * coef_sign * coef_bouyance * temp_dens * mDt;
						vel.v[IX3(i, j + 1, k)] += 0.5f * coef_sign * coef_bouyance * temp_dens * mDt;
					}
					else if (sdf.phi[IX3(i, j, k)] > 0) // && (abs(sdf.phi[IX3(i,j,k)]) > sdf.Get_grid_cell_size() * mD_RCool))
					{
						vel.v[IX3(i, j - 1, k)] += 0.5f * coef_sign * coef_bouyance * temp_dens * mDt;
						vel.v[IX3(i, j + 1, k)] += 0.5f * coef_sign * coef_bouyance * temp_dens * mDt;
					}
					else
					{
						vel.v[IX3(i, j - 1, k)] += 0.05f * coef_sign * coef_bouyance * temp_dens * mDt;
						vel.v[IX3(i, j + 1, k)] += 0.05f * coef_sign * coef_bouyance * temp_dens * mDt;
					}
				}
				else
				{
					vel.v[IX3(i, j - 1, k)] += 0.5f * coef_sign * 1.0f * coef_bouyance * temperature.value[IX3(i, j, k)] * mDt;
					vel.v[IX3(i, j + 1, k)] += 0.5f * coef_sign * 1.0f * coef_bouyance * temperature.value[IX3(i, j, k)] * mDt;
				}
				/*
				if(mSDF.sign[IX3(i,j,k)] == mSDF.GRIDY_SDF_INNER && mTargeting == 1)
				{
				coef_sign = 0.1f;
				}
				else if(mTarget_shape.sign[IX3(i,j,k)] == mSDF.GRIDY_SDF_INNER && mTargeting == 2) coef_sign = 0.1f;
				*/
				/*
				else if(mSDF.sign[IX3(i,j,k)] == mSDF.GRIDY_SDF_OUTER && density0.value[IX3(i,j,k)] > 0.5f && mTargeting == 1)
				{
				double circle_normal[3];
				double coef_normal_velocity = 10.0f;


				mSource.Get_normal(circle_normal, i, j, k);		// TO DO : 수정 요망

				velocity0.u[IX3(i-1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * mDt;
				velocity0.v[IX3(i,j-1,k)] += 2.0f * coef_normal_velocity * circle_normal[1]/2.0f * mDt;
				velocity0.w[IX3(i,j,k-1)] += coef_normal_velocity * circle_normal[2]/2.0f * mDt;

				velocity0.u[IX3(i+1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * mDt;
				velocity0.v[IX3(i,j+1,k)] += 2.0f * coef_normal_velocity * circle_normal[1]/2.0f * mDt;
				velocity0.w[IX3(i,j,k+1)] += coef_normal_velocity * circle_normal[2]/2.0f * mDt;

				}
				*/

				/*
				double coef = 0.9f;
				if(mSDF.sign[IX3(i,j,k)] == mSDF.GRIDY_SDF_INNER) coef = 1.1f;

				velocity0.u[IX3(i,j,k)] = coef * velocity.u[IX3(i,j,k)];
				velocity0.v[IX3(i,j,k)] = coef * velocity.v[IX3(i,j,k)];
				velocity0.w[IX3(i,j,k)] = coef * velocity.w[IX3(i,j,k)];
				*/
			}

}







void GRIDY_FLAME_CONTROL::Advect_fuel(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt)
{
	mLOG.In("Advect fuel");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				double x = i - dt0*vel.u[IX3(i, j, k)];
				double y = j - dt0*vel.v[IX3(i, j, k)];
				double z = k - dt0*vel.w[IX3(i, j, k)];

				if (x < 0.5)	x = 0.5f;	if (x > mWidth + 0.5) x = mWidth + 0.5f;
				if (y < 0.5)	y = 0.5f;	if (y > mHeight + 0.5) y = mHeight + 0.5f;
				if (z < 0.5)	z = 0.5f;	if (z > mDepth + 0.5) z = mDepth + 0.5f;

				int i0 = (int)x;	int i1 = i0 + 1;
				int j0 = (int)y;	int j1 = j0 + 1;
				int k0 = (int)z;	int k1 = k0 + 1;

				double s1 = x - i0;	double s0 = 1 - s1;
				double t1 = y - j0;	double t0 = 1 - t1;
				double r1 = z - k0;	double r0 = 1 - r1;

				double coef = 1.0f;
				//double coef = 0.5f;
				//double coef = (double)(2/mWidth);

				/*
				if(mTargeting > 0)
				{
				//coef = 0.999f;

				if(sdf.phi[IX3(i,j,k)] < 0) coef = 0.99f;
				else if(abs(sdf.phi[IX3(i,j,k)]) < (sdf.Get_grid_cell_size() * 2.0f)) coef = 0.99f;
				//else coef = 0.5f;
				}
				//else if(mSrc_sphere.Levelset(i,j,k) > 0) coef = 0.8f;
				*/


				/*
				//	진도보고용
				if(mCircle.Levelset(i,j,k) > 0 && mTargeting == 1) coef = 0.7f;
				else if(mSrc_sphere.Levelset(i,j,k) > 0) coef = 0.95f;
				*/
				/*
				if(mSDF.sign[IX3(i,j,k)] != mSDF.GRIDY_SDF_OUTER && mTargeting == 1) coef = 0.85;
				if(mTarget_shape.sign[IX3(i,j,k)] != mSDF.GRIDY_SDF_OUTER && mTargeting == 2) coef = 0.85;
				*/
				dens.value[IX3(i, j, k)] = r0*(s0*(t0*dens0.value[IX3(i0, j0, k0)] + t1*dens0.value[IX3(i0, j1, k0)]) + s1*(t0*dens0.value[IX3(i1, j0, k0)] + t1*dens0.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*dens0.value[IX3(i0, j0, k1)] + t1*dens0.value[IX3(i0, j1, k1)]) + s1*(t0*dens0.value[IX3(i1, j0, k1)] + t1*dens0.value[IX3(i1, j1, k1)]));
				//dens.value[IX3(i,j,k)] *= coef;
			}

	Set_boundary(0, dens.value);

	/*
	Set_boundary(1, dens.value);
	Set_boundary(2, dens.value);
	Set_boundary(3, dens.value);
	*/

	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Set_boundary(int b, double *x)
{
	bool is_top_opened = true;
	bool is_side_opened = true;


#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
	{
		for (int j = 1; j <= mDepth; j++)
		{
			x[IX3(i, 0, j)] = (b == 2 ? -x[IX3(i, 1, j)] : x[IX3(i, 1, j)]);

			if (is_top_opened) x[IX3(i, mHeight + 1, j)] = x[IX3(i, mHeight, j)] = 0.0f;
			else x[IX3(i, mHeight + 1, j)] = (b == 2 ? -x[IX3(i, mHeight, j)] : x[IX3(i, mHeight, j)]);
		}
	}

#pragma omp parallel for
	for (int i = 1; i <= mHeight; i++)
	{
		for (int j = 1; j <= mDepth; j++)
		{
			if (is_side_opened)
			{
				x[IX3(0, i, j)] = x[IX3(1, i, j)] = x[IX3(mWidth + 1, i, j)] = x[IX3(mWidth, i, j)] = 0.0f;
			}
			else
			{
				x[IX3(0, i, j)] = (b == 1 ? -x[IX3(1, i, j)] : x[IX3(1, i, j)]);
				x[IX3(mWidth + 1, i, j)] = (b == 1 ? -x[IX3(mWidth, i, j)] : x[IX3(mWidth, i, j)]);
			}
		}
	}

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
	{
		for (int j = 1; j <= mHeight; j++)
		{
			if (is_side_opened)
			{
				x[IX3(i, j, 0)] = x[IX3(i, j, 1)] = x[IX3(i, j, mDepth + 1)] = x[IX3(i, j, mDepth)] = 0.0f;
			}
			else
			{
				x[IX3(i, j, 0)] = (b == 3 ? -x[IX3(i, j, 1)] : x[IX3(i, j, 1)]);
				x[IX3(i, j, mDepth + 1)] = (b == 3 ? -x[IX3(i, j, mDepth)] : x[IX3(i, j, mDepth)]);
			}
		}
	}

	if (is_side_opened)
	{
		x[IX3(0, 0, 0)] = x[IX3(mWidth + 1, 0, 0)] = x[IX3(0, 0, mDepth + 1)] = x[IX3(mWidth + 1, 0, mDepth + 1)] = 0.0f;
	}
	else
	{
		x[IX3(0, 0, 0)] = (x[IX3(1, 0, 0)] + x[IX3(0, 1, 0)] + x[IX3(0, 0, 1)]) / 3.0f;
		x[IX3(mWidth + 1, 0, 0)] = (x[IX3(mWidth, 0, 0)] + x[IX3(mWidth + 1, 1, 0)] + x[IX3(mWidth + 1, 0, 1)]) / 3.0f;
		x[IX3(0, 0, mDepth + 1)] = (x[IX3(1, 0, mDepth + 1)] + x[IX3(0, 1, mDepth + 1)] + x[IX3(0, 0, mDepth)]) / 3.0f;
		x[IX3(mWidth + 1, 0, mDepth + 1)] = (x[IX3(mWidth, 0, mDepth + 1)] + x[IX3(mWidth + 1, 1, mDepth + 1)] + x[IX3(mWidth + 1, 0, mDepth)]) / 3.0f;
	}


	if (is_top_opened)
	{
		x[IX3(0, mHeight + 1, 0)] = x[IX3(0, mHeight + 1, mDepth + 1)] = x[IX3(mWidth + 1, mHeight + 1, 0)] = x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = 0.0f;
	}
	else
	{
		x[IX3(0, mHeight + 1, 0)] = (x[IX3(1, mHeight + 1, 0)] + x[IX3(0, mHeight, 0)] + x[IX3(0, mHeight + 1, 1)]) / 3.0f;
		x[IX3(0, mHeight + 1, mDepth + 1)] = (x[IX3(1, mHeight + 1, mDepth + 1)] + x[IX3(0, mHeight, mDepth + 1)] + x[IX3(0, mHeight + 1, mDepth)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, 0)] = (x[IX3(mWidth, mHeight + 1, 0)] + x[IX3(mWidth + 1, mHeight, 0)] + x[IX3(mWidth + 1, mHeight + 1, 1)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = (x[IX3(mWidth, mHeight + 1, mDepth + 1)] + x[IX3(mWidth + 1, mHeight, mDepth + 1)] + x[IX3(mWidth + 1, mHeight + 1, mDepth)]) / 3.0f;
	}

}




void GRIDY_FLAME_CONTROL::Advect_temperature(GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& soot, GRIDY_VECTOR_FIELD& vel, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_IMPLICIT_OBJECT_SPHERES& src_spheres, double dt)
{
	mLOG.In("Advect temperature");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

	//double C_cooling = 0.985f; //0.975f;	// 귀는 3.0에서 잘 나옴(아닌듯) -일때는 0.01f *일때는 1.0f	/ outer에서 다 쿨링할 때는 0.75였음
	//double C_cooling_outer = 1.0f;	// 테스트
	//double C_Rcool = 0.8f;	// 0.5에서 굿, 그러나 더 살려보려고 실험
	double temperature_max = 1.0f;
	//double C_TPA = 0.1f;	// 귀는 0.1에서 잘 나옴

	double res = pow((64.0f / (double)mWidth), 3);

	//printf("res : %f\n", res);

	//double cooling_coef = 20.0f; // CGI 285 org and cFire
	double cooling_coef = 50.0f; // CGI trade-off
	double reduce_coef = 1.0f;
	//if (mFrame > 50) reduce_coef = 0.0f;
	//if (mFrame > 200) reduce_coef *= (230 - mFrame) / (230 - 200);
	//if (reduce_coef < 0.0f) reduce_coef = 0.0f;

	//double cooling_coef = 10.0f; // tetra
	//double cooling_coef = 30.0f;	// dragon?



	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				double x = i - dt0*vel.u[IX3(i, j, k)];
				double y = j - dt0*vel.v[IX3(i, j, k)];
				double z = k - dt0*vel.w[IX3(i, j, k)];

				if (x < 0.5)	x = 0.5f;	if (x > mWidth + 0.5) x = mWidth + 0.5f;
				if (y < 0.5)	y = 0.5f;	if (y > mHeight + 0.5) y = mHeight + 0.5f;
				if (z < 0.5)	z = 0.5f;	if (z > mDepth + 0.5) z = mDepth + 0.5f;

				int i0 = (int)x;	int i1 = i0 + 1;
				int j0 = (int)y;	int j1 = j0 + 1;
				int k0 = (int)z;	int k1 = k0 + 1;

				double s1 = x - i0;	double s0 = 1 - s1;
				double t1 = y - j0;	double t0 = 1 - t1;
				double r1 = z - k0;	double r0 = 1 - r1;

				double coef = 1.0f;

				temperature.value[IX3(i, j, k)] = r0*(s0*(t0*temperature0.value[IX3(i0, j0, k0)] + t1*temperature0.value[IX3(i0, j1, k0)]) + s1*(t0*temperature0.value[IX3(i1, j0, k0)] + t1*temperature.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*temperature0.value[IX3(i0, j0, k1)] + t1*temperature0.value[IX3(i0, j1, k1)]) + s1*(t0*temperature0.value[IX3(i1, j0, k1)] + t1*temperature0.value[IX3(i1, j1, k1)]));

				if (temperature.value[IX3(i, j, k)] > 1.0f) temperature.value[IX3(i, j, k)] = 1.0f;
				if (fuel.value[IX3(i, j, k)] > 1.0f) fuel.value[IX3(i, j, k)] = 1.0f;

				if (mTargeting > 0)
				{
					if (sdf.phi[IX3(i, j, k)] < 0)
					{
						//if(abs(sdf.phi[IX3(i,j,k)]) <= sdf.Get_grid_cell_size() * mD_SP)
						if (abs(sdf.phi[IX3(i, j, k)]) <= sdf.Get_grid_cell_size() * mD_TP)
						{
							// 가까울수록 값이 강하도록, 멀수록 값이 약하도록, 최소 값은 0, 최대 값은 1
							double dist_RCom = ((abs(sdf.Get_grid_cell_size() * mD_TP) - abs(sdf.phi[IX3(i, j, k)])) / abs(sdf.Get_grid_cell_size() * mD_TP));


							double coef_RCom = (2.0f * dist_RCom * res);	// bunny
							//double coef_RCom = (0.5f * dist_RCom * res);	// 뭐지
							//if(coef_RCom > 1.0f) cout << "coef_RCom overflow" << endl;

							temperature.value[IX3(i, j, k)] += reduce_coef * (fuel.value[IX3(i, j, k)] * coef_RCom);
							//temperature.value[IX3(i, j, k)] += (fuel.value[IX3(i, j, k)] * coef_RCom);
							fuel.value[IX3(i, j, k)] *= (1.0f - coef_RCom);

							//temperature.value[IX3(i, j, k)] = 0.8f * temperature.value[IX3(i, j, k)] + pow(((temperature.value[IX3(i, j, k)] * 0.5f) / temperature_max), 4.0f);	// TETRA 128 #2
							if (temperature.value[IX3(i, j, k)] > 1.0f) temperature.value[IX3(i, j, k)] = 1.0f;
							temperature.value[IX3(i, j, k)] = pow(((1 / pow(temperature.value[IX3(i, j, k)], 3)) + (3 * cooling_coef * dt)), -1.0f / 3.0f);




							/*
							//dragon
							//temperature.value[IX3(i,j,k)] = 0.95f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// small flame
							temperature.value[IX3(i,j,k)] = 1.5f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// BUNNY
							//temperature.value[IX3(i,j,k)] = 0.8f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// TETRA 128 #2
							*/
							

							// others
							//temperature.value[IX3(i,j,k)] = 0.95f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// small flame
							//temperature.value[IX3(i, j, k)] = 1.0f * temperature.value[IX3(i, j, k)] + pow(((temperature.value[IX3(i, j, k)] * 0.5f) / temperature_max), 4.0f);	// BUNNY
							////temperature.value[IX3(i,j,k)] = 0.95f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// TETRA 128 #2 no.
							//temperature.value[IX3(i,j,k)] = 0.8f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// TETRA 128 #2
							//temperature.value[IX3(i,j,k)] = 0.95f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4.0f);	// teddy bear




							//temperature.value[IX3(i,j,k)] = 0.5f; // check defined region

							// 아래는 안씁니다
							//double dist_RCom = 1.0f;
							/*
							//temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.5f * res);
							//fuel.value[IX3(i,j,k)] *= 1.0f - (0.5f * res);


							//temperature.value[IX3(i,j,k)] -= 20.0f * pow(((temperature.value[IX3(i,j,k)])/ temperature_max), 4);	// 이건 이제 안씁니다

							//if(temperature.value[IX3(i,j,k)] < 0.0f) temperature.value[IX3(i,j,k)] = 0.0f;

							// 멀수록 값이 크도록, 가까울수록 값이 작도록, 최소 값은 0, 최대 값은 1
							//double dist_RCom = (abs(sdf.phi[IX3(i,j,k)])-abs(sdf.Get_grid_cell_size() * mD_TP)) / abs(sdf.Get_grid_cell_size() * mD_TP);

							//temperature.value[IX3(i,j,k)] -= 50.0f * (dist_RCom*temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)])/ temperature_max), 4));

							//temperature.value[IX3(i,j,k)] = (dist_RCom + 0.5f) * temperature.value[IX3(i,j,k)];
							//temperature.value[IX3(i,j,k)] = dist_RCom * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)] * 0.5f)/ temperature_max), 4);	// small flame
							// 0.9f
							//if(temperature.value[IX3(i,j,k)] < 0) temperature.value[IX3(i,j,k)] = 0.0f;
							*/
						}
						else // 여기가 core 부분
						{
							//temperature.value[IX3(i,j,k)] = 0.0f; // check defined region
							/*
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.5f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.5f * res);
							*/

							/*
							// cgi
							temperature.value[IX3(i, j, k)] += fuel.value[IX3(i, j, k)] * (0.2f * res);
							fuel.value[IX3(i, j, k)] *= 1.0f - (0.2f * res);
							*/

							/*
							//	cgi 287 되게 잘 됨
							double temp_coef = 25.0f;
							temp_coef *= reduce_coef;
							temperature.value[IX3(i, j, k)] += fuel.value[IX3(i, j, k)] * (temp_coef * res);
							fuel.value[IX3(i, j, k)] *= (1.0f - (temp_coef * res));
							*/
							
							//	염
							//double temp_coef = 30.0f;
							//	염2
							double temp_coef = 5.0f;
							temp_coef *= reduce_coef;

							temperature.value[IX3(i,j,k)] += 10.0f * fuel.value[IX3(i,j,k)] * (temp_coef * res);
							fuel.value[IX3(i, j, k)] *= (1.0f - (temp_coef * res));
							// 여기까지
							
							
							temperature.value[IX3(i, j, k)] = 1.0f / pow(((1 / pow(temperature.value[IX3(i, j, k)], 3)) + (3 * cooling_coef * dt)), 1.0f / 3.0f);

							//temperature.value[IX3(i,j,k)] *= 0.95f;
						}

						//temperature.value[IX3(i,j,k)] = 1.0f;
					}
					else // 쿨링을 합니다
					{
						//temperature.value[IX3(i,j,k)] = 0.0f;	//	싹 다 지웁니다
						//temperature.value[IX3(i,j,k)] = 1.0f/pow(((1/pow(temperature.value[IX3(i,j,k)], 3)) + (3 * c * dt)), 1.0f/3.0f);


						if (abs(sdf.phi[IX3(i, j, k)]) <= sdf.Get_grid_cell_size() * mD_SP)
						{
							//if(i < mWidth/2.0f) temperature.value[IX3(i,j,k)] = 1.0f; // check define region


							double dist_TPA = pow(((abs(sdf.Get_grid_cell_size() * mD_SP) - abs(sdf.phi[IX3(i, j, k)])) / abs(sdf.Get_grid_cell_size() * mD_SP)), 3.0f);
							if (dist_TPA > 1.0f) printf("dist_TPA :  %f\n", dist_TPA);

							
							// 여기가 trade-off control
							if (temperature.value[IX3(i, j, k)] > 1.0f) temperature.value[IX3(i, j, k)] = 1.0f;
							temperature.value[IX3(i, j, k)] = (0.95f * temperature.value[IX3(i, j, k)]) + (0.05f * temperature.value[IX3(i, j, k)] * dist_TPA); // MAX
							//temperature.value[IX3(i, j, k)] = (0.8f * temperature.value[IX3(i, j, k)]) + (0.2f * temperature.value[IX3(i, j, k)] * dist_TPA); //org
							//temperature.value[IX3(i, j, k)] = (0.65f * temperature.value[IX3(i, j, k)]) + (0.35f * temperature.value[IX3(i, j, k)] * dist_TPA); //middle?
							//temperature.value[IX3(i, j, k)] = (0.3f * temperature.value[IX3(i, j, k)]) + (0.7f * temperature.value[IX3(i, j, k)] * dist_TPA); // MIN
							//temperature.value[IX3(i, j, k)] = (0.95f * temperature.value[IX3(i, j, k)]) + (0.05f * temperature.value[IX3(i, j, k)] * dist_TPA); // MAX
							
							
							/*
							// cgi 283
							if (temperature.value[IX3(i, j, k)] > 1.0f) temperature.value[IX3(i, j, k)] = 1.0f;
							temperature.value[IX3(i, j, k)] = (0.8f * temperature.value[IX3(i, j, k)]) + (0.2f * temperature.value[IX3(i, j, k)] * dist_TPA);
							*/

							/*
							//	CGI
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.1f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.1f * res);

							if(temperature.value[IX3(i,j,k)] > 1.0f) temperature.value[IX3(i,j,k)] = 1.0f;
							temperature.value[IX3(i,j,k)] = (0.2f * temperature.value[IX3(i,j,k)]) + (0.8f * temperature.value[IX3(i,j,k)] * dist_TPA);
							*/

							// Tetra 128 #2
							//temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.1f * res);
							//fuel.value[IX3(i,j,k)] *= 1.0f - (0.1f * res);


							
							/*
							// Teddy bear 128
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.1f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.1f * res);

							if(temperature.value[IX3(i,j,k)] > 1.0f) temperature.value[IX3(i,j,k)] = 1.0f;
							temperature.value[IX3(i,j,k)] = (0.45f * temperature.value[IX3(i,j,k)]) + (0.55f * temperature.value[IX3(i,j,k)] * dist_TPA);
							*/

							/*
							// bunny 128 (정확히 이건 아님, 튜닝 필요)
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.1f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.1f * res);

							if(temperature.value[IX3(i,j,k)] > 1.0f) temperature.value[IX3(i,j,k)] = 1.0f;
							temperature.value[IX3(i,j,k)] = (0.1f * temperature.value[IX3(i,j,k)]) + (0.9f * temperature.value[IX3(i,j,k)] * dist_TPA);
							*/

							/*
							// dragon 128
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.1f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.1f * res);

							if(temperature.value[IX3(i,j,k)] > 1.0f) temperature.value[IX3(i,j,k)] = 1.0f;
							temperature.value[IX3(i,j,k)] = (0.5f * temperature.value[IX3(i,j,k)]) + (0.5f * temperature.value[IX3(i,j,k)] * dist_TPA);
							*/
							//temperature.value[IX3(i,j,k)] = (0.2f * temperature.value[IX3(i,j,k)]) + (0.8f * temperature.value[IX3(i,j,k)] * dist_TPA);
							//temperature.value[IX3(i,j,k)] = (1.0f * temperature.value[IX3(i,j,k)] * dist_TPA);

							//temperature.value[IX3(i,j,k)] = (0.8f * temperature.value[IX3(i,j,k)]) + (temperature.value[IX3(i,j,k)] * 0.2f * dist_TPA);
							//temperature.value[IX3(i,j,k)] *= 0.9f;
							//temperature.value[IX3(i,j,k)] = 0.9f * ((0.5f * temperature.value[IX3(i,j,k)]) + (temperature.value[IX3(i,j,k)] * 0.5f * dist_TPA));
							//temperature.value[IX3(i,j,k)] = 0.5f * (temperature.value[IX3(i,j,k)] * dist_TPA);
							//

							// 쿨링 해야 함
							temperature.value[IX3(i, j, k)] = 1.0f / pow(((1 / pow(temperature.value[IX3(i, j, k)], 3)) + (3 * cooling_coef * dt)), 1.0f / 3.0f);

							//if(i < mWidth/2.0f) temperature.value[IX3(i,j,k)] = 0.0f; // check define region

							//cooling_coef *= 500.0f;

							//temperature.value[IX3(i,j,k)] = 0.0f;	//	싹 다 지웁니다
						}
						else
						{
							temperature.value[IX3(i, j, k)] = 1.0f / pow(((1 / pow(temperature.value[IX3(i, j, k)], 3)) + (3 * cooling_coef * dt)), 1.0f / 3.0f);

							// bunny only
							//temperature.value[IX3(i,j,k)] *= 0.5f;

							//temperature.value[IX3(i,j,k)] = 0.0f;	//	싹 다 지웁니다
							//double c = 500.0f;


							/*
							if(fuel.value[IX3(i,j,k)] > 0.4f)
							{
							temperature.value[IX3(i,j,k)] += fuel.value[IX3(i,j,k)] * (0.2f * res);
							fuel.value[IX3(i,j,k)] *= 1.0f - (0.2f * res);
							}
							*/
						}

					}


					if (temperature.value[IX3(i, j, k)] <= 0.0000001f) temperature.value[IX3(i, j, k)] = 0.0f;
				}
				else // targeting이 없이는 무의미한 부분
				{
					if (sdf.phi[IX3(i, j, k)] < 0)
					{
						temperature.value[IX3(i, j, k)] += fuel.value[IX3(i, j, k)] * (0.5f * res);
						fuel.value[IX3(i, j, k)] *= (1.0f - (0.5f * res));
					}
					//else temperature.value[IX3(i,j,k)] *= 0.75f;
					else
					{
						temperature.value[IX3(i, j, k)] += fuel.value[IX3(i, j, k)] * (0.5f * res);
						fuel.value[IX3(i, j, k)] *= (1.0f - (0.5f * res));

						// cooling
						//temperature.value[IX3(i,j,k)] = 0.5f * temperature.value[IX3(i,j,k)] + pow(((temperature.value[IX3(i,j,k)]*0.5f)/ temperature_max), 4);
					}
					//if(temperature.value[IX3(i,j,k)] > 1.0f) temperature.value[IX3(i,j,k)] = 1.0f;
					//if(temperature.value[IX3(i,j,k)] < 0.0f) temperature.value[IX3(i,j,k)] = 0.0f;
				}


				if (temperature.value[IX3(i, j, k)] < 0.000001f) temperature.value[IX3(i, j, k)] = 0.0f;


				if (mMakeSoot == true)
					if (temperature.value[IX3(i, j, k)] < 1.0f) soot.value[IX3(i, j, k)] += 0.1f * temperature.value[IX3(i, j, k)];



			}

	Set_boundary(0, temperature.value);

	/*
	Set_boundary(1, dens.value);
	Set_boundary(2, dens.value);
	Set_boundary(3, dens.value);
	*/

	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Add_source_normal_force(GRIDY_IMPLICIT_OBJECT_SPHERES &spheres, GRIDY_SCALAR_FIELD &lv_spheres, GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &dens, double coef_normal_velocity, double dt)
{
	mLOG.In("Add source normal force");


	//#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				if (lv_spheres.value[IX3(i, j, k)] > 0)
				{
					double normal_sum[3], circle_normal[3];
					int sum_cnt = 0;

					for (int l = 0; l<3; l++)
						normal_sum[l] = 0.0f;

					for (int l = 0; l<spheres.mNo; l++)
					{
						if (spheres.Levelset(l, i, j, k) <= 0 && spheres.isActivate[l] == true)
						{
							//normal_sum[0] = normal_sum[2] = 0.1f;
							//normal_sum[1] = 0.1f;

							spheres.Get_normal(l, i, j, k, circle_normal);	// sphere의 중점에서 normalize된 방향 벡터를 얻어와서 적용함

							normal_sum[0] += circle_normal[0];
							normal_sum[1] += circle_normal[1];
							normal_sum[2] += circle_normal[2];

							sum_cnt++;

						}
						//else if(spheres.Levelset(l,i,j,k) == 0) cout << "untriggered sphere detected." << endl;
					}
					if (sum_cnt > 0)
					{
						//sum_cnt = 1;
						//if(sum_cnt > 1) printf("%d : %f %f %f \n", sum_cnt, normal_sum[0], normal_sum[1], normal_sum[2]);
						//printf(" %f %f %f ", vel.u[IX3(i,j,k)], vel.v[IX3(i,j,k)], vel.w[IX3(i,j,k)]);

						/*
						vel.u[IX3(i,j,k)] += coef_normal_velocity * normal_sum[0] * dt;
						vel.v[IX3(i,j,k)] += coef_normal_velocity * normal_sum[1] * dt;
						vel.w[IX3(i,j,k)] += coef_normal_velocity * normal_sum[2] * dt;
						*/
						//vel.v[IX3(i,j,k)] += 0f * dens.value[IX3(i,j,k)] * dt;
						//printf("(%f) ", coef_normal_velocity * normal_sum[0] * dt);
						double coef_dens = abs(dens.value[IX3(i, j, k)]);
						coef_dens = 1.0f;

						vel.u[IX3(i, j, k)] += coef_normal_velocity * coef_dens * normal_sum[0] / (double)sum_cnt * dt;
						vel.v[IX3(i, j, k)] += coef_normal_velocity * coef_dens * normal_sum[1] / (double)sum_cnt * dt;
						vel.w[IX3(i, j, k)] += coef_normal_velocity * coef_dens * normal_sum[2] / (double)sum_cnt * dt;

						/*
						vel.u[IX3(i-1,j,k)] += coef_normal_velocity * normal_sum[0] / 2.0f * dt;
						vel.v[IX3(i,j-1,k)] += coef_normal_velocity * normal_sum[1] / 2.0f * dt;
						vel.w[IX3(i,j,k-1)] += coef_normal_velocity * normal_sum[2] / 2.0f * dt;

						vel.u[IX3(i+1,j,k)] += coef_normal_velocity * normal_sum[0] / 2.0f * dt;
						vel.v[IX3(i,j+1,k)] += coef_normal_velocity * normal_sum[1] / 2.0f * dt;
						vel.w[IX3(i,j,k+1)] += coef_normal_velocity * normal_sum[2] / 2.0f * dt;
						*/

						//printf("-> %f %f %f \n", vel.u[IX3(i,j,k)], vel.v[IX3(i,j,k)], vel.w[IX3(i,j,k)]);

						/*
						// Add normal velocity from circle
						vel0.u[IX3(i-1,j,k)] += coef_normal_velocity * normal_sum[0]/2.0f/(double)sum_cnt * dt;
						vel0.v[IX3(i,j-1,k)] += coef_normal_velocity * normal_sum[1]/2.0f/(double)sum_cnt * dt;
						vel0.w[IX3(i,j,k-1)] += coef_normal_velocity * normal_sum[2]/2.0f/(double)sum_cnt * dt;

						vel0.u[IX3(i+1,j,k)] += coef_normal_velocity * normal_sum[0]/2.0f/(double)sum_cnt * dt;
						vel0.v[IX3(i,j+1,k)] += coef_normal_velocity * normal_sum[1]/2.0f/(double)sum_cnt * dt;
						vel0.w[IX3(i,j,k+1)] += coef_normal_velocity * normal_sum[2]/2.0f/(double)sum_cnt * dt;
						*/
					}
				}

				/*
				if(lv_spheres.value[IX3(i,j,k)] > 0)
				{
				int l = (int)lv_spheres.value[IX3(i,j,k)];

				if(spheres.Levelset(l,i,j,k) == 0 && spheres.mTrigger[l] > 0)
				{
				double circle_normal[3];

				spheres.Get_normal(l, i, j, k, circle_normal);

				// Add normal velocity from circle
				vel0.u[IX3(i-1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * dt;
				vel0.v[IX3(i,j-1,k)] += coef_normal_velocity * circle_normal[1]/2.0f * dt;
				vel0.w[IX3(i,j,k-1)] += coef_normal_velocity * circle_normal[2]/2.0f * dt;

				vel0.u[IX3(i+1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * dt;
				vel0.v[IX3(i,j+1,k)] += coef_normal_velocity * circle_normal[1]/2.0f * dt;
				vel0.w[IX3(i,j,k+1)] += coef_normal_velocity * circle_normal[2]/2.0f * dt;
				}
				}
				*/

			}



	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Compute_potential_force(GRIDY_VECTOR_FIELD &potential, GRIDY_SIGNED_DISTANCE_FIELD &sdf, double threshold)
{
	mLOG.In("Compute potential force");

	double h = 1.0f / sdf.Get_grid_cell_size();

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				potential.u[IX3(i, j, k)] = -0.5f * h *(sdf.phi[IX3(i + 1, j, k)] - sdf.phi[IX3(i - 1, j, k)]);
				potential.v[IX3(i, j, k)] = -0.5f * h *(sdf.phi[IX3(i, j + 1, k)] - sdf.phi[IX3(i, j - 1, k)]);
				potential.w[IX3(i, j, k)] = -0.5f * h *(sdf.phi[IX3(i, j, k + 1)] - sdf.phi[IX3(i, j, k - 1)]);


				//	error correction 
				//	너무 튀는 값은 0으로 조정
				if (abs(potential.u[IX3(i, j, k)]) > threshold || abs(potential.v[IX3(i, j, k)]) > threshold || abs(potential.w[IX3(i, j, k)]) > threshold)
				{
					potential.u[IX3(i, j, k)] = potential.v[IX3(i, j, k)] = potential.w[IX3(i, j, k)] = 0.0f;
				}


				//	normalize
				double norm = sqrt(potential.u[IX3(i, j, k)] * potential.u[IX3(i, j, k)] + potential.v[IX3(i, j, k)] * potential.v[IX3(i, j, k)] + potential.w[IX3(i, j, k)] * potential.w[IX3(i, j, k)]);

				if (norm > 0.0f)
				{
					potential.u[IX3(i, j, k)] = potential.u[IX3(i, j, k)] / norm;
					potential.v[IX3(i, j, k)] = potential.v[IX3(i, j, k)] / norm;
					potential.w[IX3(i, j, k)] = potential.w[IX3(i, j, k)] / norm;
				}

				//printf("%f %f %f\n", potential.u[IX3(i,j,k)], potential.v[IX3(i,j,k)], potential.w[IX3(i,j,k)] );
			}

	mLOG.Out();
}





void GRIDY_FLAME_CONTROL::Add_potential_force(GRIDY_VECTOR_FIELD &vel, GRIDY_VECTOR_FIELD &potential, GRIDY_SIGNED_DISTANCE_FIELD& sdf, GRIDY_SCALAR_FIELD& dens, double coef_potential, double dt)
{
	mLOG.In("Add potential force");


	//	add forces
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				//if(abs(sdf.phi[IX3(i,j,k)]) < 2.0f * sdf.Get_grid_cell_size())

				if (sdf.phi[IX3(i, j, k)] < 0.0f)
				{
					double coef_dens = abs(dens.value[IX3(i, j, k)]);
					coef_dens = 1.0f;
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
				else
				{
					//if(abs(sdf.phi[IX3(i,j,k)]) <= abs(sdf.Get_grid_cell_size() * mD_SP))
					if (abs(sdf.phi[IX3(i, j, k)]) <= abs(sdf.Get_grid_cell_size() * mD_SP))
					{
						double coef_dens = abs(dens.value[IX3(i, j, k)]);
						coef_dens = 1.0f;
						// 여기에 density를 넣어보면 어떨까?
						double coef_surface = 2.0f * coef_potential;

						vel.u[IX3(i, j, k)] -= coef_surface * coef_dens * potential.u[IX3(i, j, k)] * dt;
						vel.v[IX3(i, j, k)] -= coef_surface * coef_dens * potential.v[IX3(i, j, k)] * dt;
						vel.w[IX3(i, j, k)] -= coef_surface * coef_dens * potential.w[IX3(i, j, k)] * dt;
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





void GRIDY_FLAME_CONTROL::Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, GRIDY_SIGNED_DISTANCE_FIELD &sdf, double coef_vorticity, double dt)
{
	mLOG.In("Add vorticity confinement");

	GRIDY_VECTOR_FIELD vortex(vel.Get_width_size(), vel.Get_height_size(), vel.Get_depth_size());

	//	compute vortex
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				double DwDy = (vel.w[IX3(i, j + 1, k)] - vel.w[IX3(i, j - 1, k)]) * 0.5f;	// central difference, with cell length of 1
				double DvDz = (vel.v[IX3(i, j, k + 1)] - vel.v[IX3(i, j, k - 1)]) * 0.5f;
				double DuDz = (vel.u[IX3(i, j, k + 1)] - vel.u[IX3(i, j, k - 1)]) * 0.5f;
				double DwDx = (vel.w[IX3(i + 1, j, k)] - vel.w[IX3(i - 1, j, k)]) * 0.5f;
				double DvDx = (vel.v[IX3(i + 1, j, k)] - vel.v[IX3(i - 1, j, k)]) * 0.5f;
				double DuDy = (vel.u[IX3(i, j + 1, k)] - vel.u[IX3(i, j - 1, k)]) * 0.5f;

				vortex.u[IX3(i, j, k)] = DwDy - DvDz;
				vortex.v[IX3(i, j, k)] = DuDz - DwDx;
				vortex.w[IX3(i, j, k)] = DvDx - DuDy;
			}



	//	compute vorticity confinement
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				//if(sdf.phi[IX3(i,j,k)] > 0) // || abs(sdf.phi[IX3(i,j,k)]) <= mD_SP)
				{
					double DoDx = (vortex.Get_norm2(i + 1, j, k) - vortex.Get_norm2(i - 1, j, k)) * 0.5f;
					double DoDy = (vortex.Get_norm2(i, j + 1, k) - vortex.Get_norm2(i, j - 1, k)) * 0.5f;
					double DoDz = (vortex.Get_norm2(i, j, k + 1) - vortex.Get_norm2(i, j, k - 1)) * 0.5f;


					double len = sqrt(DoDx*DoDx + DoDy*DoDy + DoDz*DoDz);

					if (len != 0)
					{
						DoDx /= len;
						DoDy /= len;
						DoDz /= len;
					}

					vel.u[IX3(i, j, k)] += (DoDy*vortex.w[IX3(i, j, k)] - vortex.v[IX3(i, j, k)] * DoDz) * coef_vorticity * dt;
					vel.v[IX3(i, j, k)] += (DoDz*vortex.u[IX3(i, j, k)] - vortex.w[IX3(i, j, k)] * DoDx) * coef_vorticity * dt;
					vel.w[IX3(i, j, k)] += (DoDx*vortex.v[IX3(i, j, k)] - vortex.u[IX3(i, j, k)] * DoDy) * coef_vorticity * dt;
				}
			}



	mLOG.Out();
}





/*
void GRIDY_FLAME_CONTROL::Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p)
{
mLOG.In("Project");

double h = (double)(1.0f / mWidth);
if((double)(1.0f/mHeight) > h) h = (double)(1.0f / mHeight);
if((double)(1.0f/mDepth) > h) h = (double)(1.0f / mDepth);


//	compute divergence
#pragma omp parallel for
for(int i=1; i<=mWidth; i++)
for(int j=1; j<=mHeight; j++)
for(int k=1; k<=mDepth; k++)
{
div[IX3(i,j,k)] = -0.5f * h *(vel.u[IX3(i+1,j,k)] - vel.u[IX3(i-1,j,k)] + vel.v[IX3(i,j+1,k)] - vel.v[IX3(i,j-1,k)] + vel.w[IX3(i,j,k+1)] - vel.w[IX3(i,j,k-1)] );
p[IX3(i,j,k)] = 0;
}

Set_boundary(0, div);
Set_boundary(0, p);


Gause_seidel_solver(0, p, div, 1, 6);


//	subtract gradient field
#pragma omp parallel for
for(int i=1; i<=mWidth; i++)
for(int j=1; j<=mHeight; j++)
for(int k=1; k<=mDepth; k++)
{
vel.u[IX3(i,j,k)] -= 0.5f *(p[IX3(i+1, j, k)] - p[IX3(i-1, j, k)])/h;
vel.v[IX3(i,j,k)] -= 0.5f *(p[IX3(i, j+1, k)] - p[IX3(i, j-1, k)])/h;
vel.w[IX3(i,j,k)] -= 0.5f *(p[IX3(i, j, k+1)] - p[IX3(i, j, k-1)])/h;
}

Set_boundary(1, vel.u);
Set_boundary(2, vel.v);
Set_boundary(3, vel.w);

mLOG.Out();

}





void GRIDY_FLAME_CONTROL::Set_boundary(int b, double *x)
{
bool is_top_opened = true;



#pragma omp parallel for
for(int i=1; i<=mWidth; i++)
{
for(int j=1; j<=mDepth; j++)
{

//x[IX3(i,0,j)] = (b==2 ? -x[IX3(i,1,j)] : x[IX3(i,1,j)]);

//if(is_top_opened) x[IX3(i,mHeight+1,j)] = x[IX3(i,mHeight,j)]= 0;
//else x[IX3(i,mHeight+1,j)] = (b==2 ? -x[IX3(i,mHeight,j)] : x[IX3(i,mHeight,j)]);

x[IX3(i,0,j)] = x[IX3(i,mHeight+1,j)] = 0;
}
}

#pragma omp parallel for
for(int i=1; i<=mHeight; i++)
{
for(int j=1; j<=mDepth; j++)
{

//x[IX3(0,i,j)] = (b==1 ? -x[IX3(1,i,j)] : x[IX3(1,i,j)]);
//x[IX3(mWidth+1,i,j)] = (b==1 ? -x[IX3(mWidth,i,j)] : x[IX3(mWidth,i,j)]);

x[IX3(0,i,j)] = x[IX3(mWidth+1,i,j)] = 0;
}
}

#pragma omp parallel for
for(int i=1; i<=mWidth; i++)
{
for(int j=1; j<=mHeight; j++)
{

//x[IX3(i,j,0)] = (b==3 ? -x[IX3(i,j,1)] : x[IX3(i,j,1)]);
//x[IX3(i,j,mDepth+1)] = (b==3 ? -x[IX3(i,j,mDepth)] : x[IX3(i,j,mDepth)]);

x[IX3(i,j,0)] = x[IX3(i,j,mDepth+1)] = 0;
}
}

x[IX3(0,mHeight+1,0)] = x[IX3(0,mHeight+1,mDepth+1)] = x[IX3(mWidth+1,mHeight+1,0)] = x[IX3(mWidth+1,mHeight+1,mDepth+1)] = 0;

x[IX3(0,0,0)] = x[IX3(mWidth+1,0,0)] = x[IX3(0,0,mDepth+1)] = x[IX3(mWidth+1,0,mDepth+1)] = 0;
}
*/





void GRIDY_FLAME_CONTROL::Write_data(int frame)
{
	mLOG.In("Write simulated data");


	mLOG.In("Writing Density");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(density.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	else mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	mLOG.Out();

	mLOG.In("Writing Temperature");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(temperature.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "temperature");
	else mFILE_IO.Write_data(temperature.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "temperature");
	mLOG.Out();

	if (mMakeSoot == true)
	{
		mLOG.In("Writing Soot");
		if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(soot.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "soot");
		else mFILE_IO.Write_data(soot.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "soot");
		mLOG.Out();
	}

	mLOG.In("Writing Velocity");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY)
	{
		mFILE_IO.Write_data_binary(velocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data_binary(velocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data_binary(velocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	}
	else
	{
		mFILE_IO.Write_data(velocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data(velocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data(velocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	}
	mLOG.Out();

	/*
	mLOG.In("Writing Sphere Position Data");
	mFILE_IO.Write_sphere_data((int)mCircle.center_x, (int)mCircle.center_y, (int)mCircle.center_z, mCircle.radius, mSimulation_path, frame, "sphere");
	mLOG.Out();
	*/

	mLOG.Out();
}










#endif




