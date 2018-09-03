//	
//	GRIDY_SIMPLIFIED_FIRE.h
//
//	This file is a part of the GRIDY library.
//	Designed by TaeHyeong Kim (usemagic@gmail.com)
//	and Contributed by Eunki Hong and Nuri Park.
//





#ifndef GRIDY_SIMPLIFIED_FIRE_H
#define GRIDY_SIMPLIFIED_FIRE_H





#include "GRIDY_SCALAR_FIELD.h"
#include "GRIDY_VECTOR_FIELD.h"
#include "GRIDY_FILE_IO.h"
#include "GRIDY_LOG.h"
#include "GRIDY_SIGNED_DISTANCE_FIELD.h"
#include "GRIDY_IMPLICIT_OBJECTS.h"
#include "GRIDY_PARSER.h"




class GRIDY_FIRE 
{
protected:
	int mStart_frame;
	int mEnd_frame;

	double mDt, mDiffuse, mViscosity;

	int mWidth, mHeight, mDepth;

	string mSimulation_path;
	int mData_file_type;

	double mCombustion_ratio;
	double mCoef_normal_force, mCoef_buoyancy, mCoef_vorticity, mCoef_cooling;

	double mMaxMax_radius = 5.0f;





protected:
	GRIDY_LOG mLOG;
	GRIDY_FILE_IO mFILE_IO;
	GRIDY_COMMON mUTIL;
	GRIDY_PARSER mPARSER;

	/*GRIDY_IMPLICIT_OBJECT_SPHERES mSrc_sphere;*/
	GRIDY_IMPLICIT_OBJECT_SPHERES mSrc_sphere;

	GRIDY_SCALAR_FIELD mDensity, mDensity0;
	GRIDY_SCALAR_FIELD mTemperature, mTemperature0;

	GRIDY_VECTOR_FIELD mVelocity, mVelocity0;

	//GRIDY_VECTOR_FIELD potential;
	//GRIDY_SCALAR_FIELD phi;





public:
	GRIDY_FIRE(void);
	~GRIDY_FIRE(void);

	void Set_parameters(double dt, double diffuse, double viscosity, double combustion_ratio = 0.1f, double normal_force = 0.5f, double buoyancy = 0.1f, double vorticity = 1.0f, double cooling = 0.2f);

	void Init(int width, int height, int depth, string path);
	void Init(int equally_size, string path);
	void Init_by_file(string config_file);

	void Simulate(int start_frame, int end_frame);



protected:
	void Clear_previous_step_data(void);

	void Add_source(GRIDY_SCALAR_FIELD& target, GRIDY_SCALAR_FIELD& source, double dt);
	void Add_velocity(GRIDY_VECTOR_FIELD& target, GRIDY_VECTOR_FIELD& source, double dt);

	void Advect_velocity(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt);
	void Advect_density(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt);
	void Advect_temperature(GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_SCALAR_FIELD& fuel, GRIDY_VECTOR_FIELD& vel, double dt);

	void Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p);
	void Project_using_GPU(GRIDY_VECTOR_FIELD& vel, double *div, double *p);

	void Diffuse(int b, double *x, double *x0, double dt);
	void Gause_seidel_solver(int b, double *x, double *x0, double a, double c);

	void Set_boundary(int b, double *x);

	void Velocity_step(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, GRIDY_SCALAR_FIELD& tem, double dt);
	void Density_step(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& fuel0, GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_VECTOR_FIELD& vel, double diff, double dt);
	void Combustion(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& temperature, double ratio, double dt);
	void Cooling(GRIDY_SCALAR_FIELD& temperature, double coef_cooling, double dt);

	void Write_data(int frame);

	void Add_normal_force(GRIDY_VECTOR_FIELD &vel, double coef_normal_velocity, double dt);
	void Add_buoyancy(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &tem, double coef_buoyance, double dt);

	void Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, double epsilon, double dt);





	inline void SWAP_value(GRIDY_SCALAR_FIELD& target, GRIDY_SCALAR_FIELD& source)
	{
		double *temp = target.value;
		target.value = source.value;
		source.value = temp;
	};

	inline void SWAP_value(GRIDY_VECTOR_FIELD& target, GRIDY_VECTOR_FIELD& source)
	{
		double *temp;
		temp = target.u;	target.u = source.u;	source.u = temp;
		temp = target.v;	target.v = source.v;	source.v = temp;
		temp = target.w;	target.w = source.w;	source.w = temp;
	};
};







GRIDY_FIRE::GRIDY_FIRE(void)
{
	//	set parameters
	Set_parameters(0.1f, 0.0f, 0.0f, 0.3f, 8.0f, 0.815f, 0.0f, 12.0f);

	mStart_frame = mEnd_frame = 0;
	mWidth = mHeight = mDepth = 0;

	mSimulation_path = "C:\\";

	mData_file_type = mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY;
	//mData_file_type = mFILE_IO.GRIDY_DATA_FILE_TYPE_TEXT;
}




GRIDY_FIRE::~GRIDY_FIRE(void)
{
	// terminate this program.
}





void GRIDY_FIRE::Set_parameters(double dt, double diffuse, double viscosity, double combustion_ratio, double normal_force, double buoyancy, double vorticity, double cooling)
{
	mDt = dt;
	mDiffuse = diffuse;
	mViscosity = viscosity;

	mCombustion_ratio = combustion_ratio;
	mCoef_normal_force = normal_force;
	mCoef_buoyancy = buoyancy;
	mCoef_vorticity = vorticity;
	mCoef_cooling = cooling;
}





void GRIDY_FIRE::Init(int width, int height, int depth, string path = "C:\\")
{
	mLOG.Set_file_out(path, "GRIDY_FIRE_STANDARD", 4, 2);

	mLOG.In("INITIALIZE SIMULATION");

	mWidth = width;
	mHeight = height;
	mDepth = depth;

	mSimulation_path = path;

	mDensity.Init(mWidth, mHeight, mDepth);
	mDensity0.Init(mWidth, mHeight, mDepth);
	mTemperature.Init(mWidth, mHeight, mDepth);
	mTemperature0.Init(mWidth, mHeight, mDepth);

	mVelocity.Init(mWidth, mHeight, mDepth);
	mVelocity0.Init(mWidth, mHeight, mDepth);

	mLOG.Out();
}





void GRIDY_FIRE::Init(int equally_size, string path = "C:\\")
{
	Init(equally_size, equally_size, equally_size, path);
}




void GRIDY_FIRE::Init_by_file(string config_file)
{
	mPARSER.parse_setup(config_file);

	mPARSER.get("width", &mWidth);
	mPARSER.get("height", &mHeight);
	mPARSER.get("depth", &mDepth);

	mPARSER.get("start_frame", &mStart_frame);
	mPARSER.get("end_frame", &mEnd_frame);

	mPARSER.get("path", &mSimulation_path);

	mPARSER.get("combustion_ratio", &mCombustion_ratio);
	mPARSER.get("coef_normal_force", &mCoef_normal_force);
	mPARSER.get("coef_buoyancy", &mCoef_buoyancy);
	mPARSER.get("coef_vorticity", &mCoef_vorticity);
	mPARSER.get("coef_cooling", &mCoef_cooling);
	mPARSER.get("dt", &mDt);

	mDensity.Init(mWidth, mHeight, mDepth);
	mDensity0.Init(mWidth, mHeight, mDepth);

	mVelocity.Init(mWidth, mHeight, mDepth);
	mVelocity0.Init(mWidth, mHeight, mDepth);

	mTemperature.Init(mWidth, mHeight, mDepth);
	mTemperature0.Init(mWidth, mHeight, mDepth);

	string log_path;
	mPARSER.get("log_path", &log_path);

	mLOG.Set_file_out(log_path, "GRIDY_SIMPLIFIED_FIRE", 4, 2);

	mPARSER.show_all_items();
}





void GRIDY_FIRE::Clear_previous_step_data(void)
{
	mLOG.In("Clear data");

	mDensity0.Clear_data();
	mTemperature0.Clear_data();
	mVelocity0.Clear_data();

	mLOG.Out();
}





void GRIDY_FIRE::Add_source(GRIDY_SCALAR_FIELD& target, GRIDY_SCALAR_FIELD& source, double dt)
{
	mLOG.In("Add source");

	if (target.Get_size() != source.Get_size())
	{
		cout << "\n\nERROR : Add_source\n\n";

	}

	#pragma omp parallel for
	for (int i = 0; i<target.Get_size(); i++)
	{
		target.value[i] += source.value[i];
	}

	mLOG.Out();
}





void GRIDY_FIRE::Add_velocity(GRIDY_VECTOR_FIELD& target, GRIDY_VECTOR_FIELD& source, double dt)
{
	mLOG.In("Add velocity");

	if (target.Get_size() != source.Get_size())
	{
		string err_msg = "ERROR : Add_source";
		mUTIL.Terminate_system(err_msg);
	}

	#pragma omp parallel for
	for (int i = 0; i<target.Get_size(); i++)
	{
		target.u[i] += dt * source.u[i];
		target.v[i] += dt * source.v[i];
		target.w[i] += dt * source.w[i];
	}

	mLOG.Out();
}





void GRIDY_FIRE::Advect_velocity(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt)
{
	mLOG.In("Advect velocity");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				double x = i - dt0*vel0.u[IX3(i, j, k)];
				double y = j - dt0*vel0.v[IX3(i, j, k)];
				double z = k - dt0*vel0.w[IX3(i, j, k)];

				if (x < 0.5)	x = 0.5f; if (x > mWidth + 0.5) x = mWidth + 0.5f;
				if (y < 0.5)	y = 0.5f; if (y > mHeight + 0.5) y = mHeight + 0.5f;
				if (z < 0.5)	z = 0.5f; if (z > mDepth + 0.5) z = mDepth + 0.5f;

				int i0 = (int)x;	int i1 = i0 + 1;
				int j0 = (int)y;	int j1 = j0 + 1;
				int k0 = (int)z;	int k1 = k0 + 1;

				double s1 = x - i0;	double s0 = 1 - s1;
				double t1 = y - j0;	double t0 = 1 - t1;
				double r1 = z - k0;	double r0 = 1 - r1;

				vel.u[IX3(i, j, k)] = r0*(s0*(t0*vel0.u[IX3(i0, j0, k0)] + t1*vel0.u[IX3(i0, j1, k0)]) + s1*(t0*vel0.u[IX3(i1, j0, k0)] + t1*vel0.u[IX3(i1, j1, k0)])) + r1*(s0*(t0*vel0.u[IX3(i0, j0, k1)] + t1*vel0.u[IX3(i0, j1, k1)]) + s1*(t0*vel0.u[IX3(i1, j0, k1)] + t1*vel0.u[IX3(i1, j1, k1)]));
				vel.v[IX3(i, j, k)] = r0*(s0*(t0*vel0.v[IX3(i0, j0, k0)] + t1*vel0.v[IX3(i0, j1, k0)]) + s1*(t0*vel0.v[IX3(i1, j0, k0)] + t1*vel0.v[IX3(i1, j1, k0)])) + r1*(s0*(t0*vel0.v[IX3(i0, j0, k1)] + t1*vel0.v[IX3(i0, j1, k1)]) + s1*(t0*vel0.v[IX3(i1, j0, k1)] + t1*vel0.v[IX3(i1, j1, k1)]));
				vel.w[IX3(i, j, k)] = r0*(s0*(t0*vel0.w[IX3(i0, j0, k0)] + t1*vel0.w[IX3(i0, j1, k0)]) + s1*(t0*vel0.w[IX3(i1, j0, k0)] + t1*vel0.w[IX3(i1, j1, k0)])) + r1*(s0*(t0*vel0.w[IX3(i0, j0, k1)] + t1*vel0.w[IX3(i0, j1, k1)]) + s1*(t0*vel0.w[IX3(i1, j0, k1)] + t1*vel0.w[IX3(i1, j1, k1)]));
			}


	Set_boundary(1, vel.u);
	Set_boundary(2, vel.v);
	Set_boundary(3, vel.w);

	mLOG.Out();
}





void GRIDY_FIRE::Advect_density(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt)
{
	mLOG.In("Advect density");

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
	
				dens.value[IX3(i, j, k)] = r0*(s0*(t0*dens0.value[IX3(i0, j0, k0)] + t1*dens0.value[IX3(i0, j1, k0)]) + s1*(t0*dens0.value[IX3(i1, j0, k0)] + t1*dens0.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*dens0.value[IX3(i0, j0, k1)] + t1*dens0.value[IX3(i0, j1, k1)]) + s1*(t0*dens0.value[IX3(i1, j0, k1)] + t1*dens0.value[IX3(i1, j1, k1)]));
				dens.value[IX3(i, j, k)] *= coef;
			}

	//Set_boundary(0, dens.value);

	/*
	Set_boundary(1, dens.value);
	Set_boundary(2, dens.value);
	Set_boundary(3, dens.value);
	*/

	mLOG.Out();
}



void GRIDY_FIRE::Gause_seidel_solver(int b, double *x, double *x0, double a, double c)
{
	mLOG.In("Execute Gause-seidel solver");

	for (int l = 0; l<20; l++)
	{
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++)
				{
					x[IX3(i, j, k)] = (x0[IX3(i, j, k)] + a*(x[IX3(i - 1, j, k)] + x[IX3(i + 1, j, k)] + x[IX3(i, j - 1, k)] + x[IX3(i, j + 1, k)] + x[IX3(i, j, k - 1)] + x[IX3(i, j, k + 1)])) / c;
				}

		Set_boundary(b, x);
	}

	mLOG.Out();
}









void GRIDY_FIRE::Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p)
{
	mLOG.In("Project");

	double h = (double)(1.0f / mWidth);
	if ((double)(1.0f / mHeight) > h) h = (double)(1.0f / mHeight);
	if ((double)(1.0f / mDepth) > h) h = (double)(1.0f / mDepth);


	//	compute divergence
	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				div[IX3(i, j, k)] = -0.5f * h *(vel.u[IX3(i + 1, j, k)] - vel.u[IX3(i - 1, j, k)] + vel.v[IX3(i, j + 1, k)] - vel.v[IX3(i, j - 1, k)] + vel.w[IX3(i, j, k + 1)] - vel.w[IX3(i, j, k - 1)]);
				p[IX3(i, j, k)] = 0;
			}

	Set_boundary(0, div);
	Set_boundary(0, p);


	Gause_seidel_solver(0, p, div, 1, 6);


	//	subtract gradient field
	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				vel.u[IX3(i, j, k)] -= 0.5f *(p[IX3(i + 1, j, k)] - p[IX3(i - 1, j, k)]) / h;
				vel.v[IX3(i, j, k)] -= 0.5f *(p[IX3(i, j + 1, k)] - p[IX3(i, j - 1, k)]) / h;
				vel.w[IX3(i, j, k)] -= 0.5f *(p[IX3(i, j, k + 1)] - p[IX3(i, j, k - 1)]) / h;
			}

	Set_boundary(1, vel.u);
	Set_boundary(2, vel.v);
	Set_boundary(3, vel.w);

	mLOG.Out();
}










void GRIDY_FIRE::Velocity_step(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, GRIDY_SCALAR_FIELD& tem, double dt)
{
	mLOG.In("Velocity step");


	Add_velocity(vel, vel0, dt);						// Add velocity0 to velocity
	Add_normal_force(vel, mCoef_normal_force, dt);		
	Add_buoyancy(vel, tem, mCoef_buoyancy, dt);
	Add_vorticity_confinement(vel,mCoef_vorticity, dt);


	/*
	SWAP_value(vel, vel0);
	//	Diffuse velocity (use only if we need)
	mLOG.In("Diffuse velocity");
		Diffuse(1, vel.u, vel0.u, mDt);
		Diffuse(2, vel.v, vel0.v, mDt);
		Diffuse(3, vel.w, vel0.w, mDt);
	mLOG.Out();
	//SWAP_value(vel, vel0);	// if not use diffuse velocity
	*/


	Project(vel, vel0.u, vel0.v);
	SWAP_value(vel, vel0);

	Advect_velocity(vel, vel0, dt);
	Project(vel, vel0.u, vel0.v);

	mLOG.Out();
}











void GRIDY_FIRE::Diffuse(int b, double *x, double *x0, double dt)
{
	double a = mDt * mDiffuse * mWidth * mHeight * mDepth;

	Gause_seidel_solver(b, x, x0, a, 1 + 6 * a);
}




void GRIDY_FIRE::Density_step(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& fuel0, GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_VECTOR_FIELD& vel, double diff, double dt)
{
	mLOG.In("Density step");



	Add_source(fuel, fuel0, dt);
	SWAP_value(fuel0, fuel);
	Diffuse(0, fuel.value, fuel0.value, dt);
	SWAP_value(fuel0, fuel);
	Advect_density(fuel, fuel0, vel, dt);


	
	Add_source(temperature, temperature0, dt);
	SWAP_value(temperature0, temperature);
	Diffuse(0, temperature.value, temperature0.value, dt);
	SWAP_value(temperature0, temperature);
	Advect_temperature(temperature, temperature0, fuel, vel, dt);
	

	mLOG.Out();
}









void GRIDY_FIRE::Set_boundary(int b, double *x)
{
	bool is_top_opened = true;



	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
	{
		for (int j = 1; j <= mDepth; j++)
		{
			x[IX3(i, 0, j)] = (b == 2 ? -x[IX3(i, 1, j)] : x[IX3(i, 1, j)]);

			if (is_top_opened) x[IX3(i, mHeight + 1, j)] = x[IX3(i, mHeight, j)] = 0;
			else x[IX3(i, mHeight + 1, j)] = (b == 2 ? -x[IX3(i, mHeight, j)] : x[IX3(i, mHeight, j)]);
		}
	}

	#pragma omp parallel for
	for (int i = 1; i <= mHeight; i++)
	{
		for (int j = 1; j <= mDepth; j++)
		{
			x[IX3(0, i, j)] = (b == 1 ? -x[IX3(1, i, j)] : x[IX3(1, i, j)]);
			x[IX3(mWidth + 1, i, j)] = (b == 1 ? -x[IX3(mWidth, i, j)] : x[IX3(mWidth, i, j)]);
		}
	}

	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
	{
		for (int j = 1; j <= mHeight; j++)
		{
			x[IX3(i, j, 0)] = (b == 3 ? -x[IX3(i, j, 1)] : x[IX3(i, j, 1)]);
			x[IX3(i, j, mDepth + 1)] = (b == 3 ? -x[IX3(i, j, mDepth)] : x[IX3(i, j, mDepth)]);
		}
	}


	x[IX3(0, 0, 0)] = (x[IX3(1, 0, 0)] + x[IX3(0, 1, 0)] + x[IX3(0, 0, 1)]) / 3.0f;
	x[IX3(mWidth + 1, 0, 0)] = (x[IX3(mWidth, 0, 0)] + x[IX3(mWidth + 1, 1, 0)] + x[IX3(mWidth + 1, 0, 1)]) / 3.0f;
	x[IX3(0, 0, mDepth + 1)] = (x[IX3(1, 0, mDepth + 1)] + x[IX3(0, 1, mDepth + 1)] + x[IX3(0, 0, mDepth)]) / 3.0f;
	x[IX3(mWidth + 1, 0, mDepth + 1)] = (x[IX3(mWidth, 0, mDepth + 1)] + x[IX3(mWidth + 1, 1, mDepth + 1)] + x[IX3(mWidth + 1, 0, mDepth)]) / 3.0f;

	if (is_top_opened)
	{
		x[IX3(0, mHeight + 1, 0)] = x[IX3(0, mHeight + 1, mDepth + 1)] = x[IX3(mWidth + 1, mHeight + 1, 0)] = x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = 0;
	}
	else
	{
		x[IX3(0, mHeight + 1, 0)] = (x[IX3(1, mHeight + 1, 0)] + x[IX3(0, mHeight, 0)] + x[IX3(0, mHeight + 1, 1)]) / 3.0f;
		x[IX3(0, mHeight + 1, mDepth + 1)] = (x[IX3(1, mHeight + 1, mDepth + 1)] + x[IX3(0, mHeight, mDepth + 1)] + x[IX3(0, mHeight + 1, mDepth)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, 0)] = (x[IX3(mWidth, mHeight + 1, 0)] + x[IX3(mWidth + 1, mHeight, 0)] + x[IX3(mWidth + 1, mHeight + 1, 1)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = (x[IX3(mWidth, mHeight + 1, mDepth + 1)] + x[IX3(mWidth + 1, mHeight, mDepth + 1)] + x[IX3(mWidth + 1, mHeight + 1, mDepth)]) / 3.0f;
	}
}









void GRIDY_FIRE::Write_data(int frame)
{
	mLOG.In("Write simulated data");


	mLOG.In("Writing Density");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(mDensity.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	else mFILE_IO.Write_data(mDensity.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	mLOG.Out();

	mLOG.In("Writing Velocity");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY)
	{
		mFILE_IO.Write_data_binary(mVelocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data_binary(mVelocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data_binary(mVelocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	}
	else
	{
		mFILE_IO.Write_data(mVelocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data(mVelocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data(mVelocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	}
	mLOG.Out();

	mLOG.In("Writing Temperature");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(mTemperature.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "temperature");
	else mFILE_IO.Write_data(mDensity.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "temperature");
	mLOG.Out();



	mLOG.Out();
}








void GRIDY_FIRE::Add_normal_force(GRIDY_VECTOR_FIELD &vel, double coef_normal_velocity, double dt)
{
	mLOG.In("Add normal force");

	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				if (mSrc_sphere.Levelset(0, i, j, k) == 0)
				{
					double circle_normal[3];

					for (int l = 0; l<3; l++)
						circle_normal[l] = 0.0f;

					mSrc_sphere.Get_normal(0, i, j, k, circle_normal);
					
					// Add normal velocity from circle
					/*
					vel0.u[IX3(i-1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * dt;
					vel0.v[IX3(i,j-1,k)] += coef_normal_velocity * circle_normal[1]/2.0f * dt;
					vel0.w[IX3(i,j,k-1)] += coef_normal_velocity * circle_normal[2]/2.0f * dt;

					vel0.u[IX3(i+1,j,k)] += coef_normal_velocity * circle_normal[0]/2.0f * dt;
					vel0.v[IX3(i,j+1,k)] += coef_normal_velocity * circle_normal[1]/2.0f * dt;
					vel0.w[IX3(i,j,k+1)] += coef_normal_velocity * circle_normal[2]/2.0f * dt;
					*/

					vel.u[IX3(i, j, k)] += coef_normal_velocity * circle_normal[0] * dt;
					vel.v[IX3(i, j, k)] += coef_normal_velocity * circle_normal[1] * dt;
					vel.w[IX3(i, j, k)] += coef_normal_velocity * circle_normal[2] * dt;
				}
			}

	mLOG.Out();
}











void GRIDY_FIRE::Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, double epsilon, double dt)
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

				vel.u[IX3(i, j, k)] += (DoDy*vortex.w[IX3(i, j, k)] - vortex.v[IX3(i, j, k)] * DoDz) * epsilon * dt;
				vel.v[IX3(i, j, k)] += (DoDz*vortex.u[IX3(i, j, k)] - vortex.w[IX3(i, j, k)] * DoDx) * epsilon * dt;
				vel.w[IX3(i, j, k)] += (DoDx*vortex.v[IX3(i, j, k)] - vortex.u[IX3(i, j, k)] * DoDy) * epsilon * dt;
			}

	mLOG.Out();
}

//	Add bouyance function must be execute before projection
void GRIDY_FIRE::Add_buoyancy(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &tem, double coef_buoyance, double dt)
{
	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				// is in inside?
				double coef_sign = 1.0f;

				//vel.v[IX3(i, j, k)] += 0.5f * coef_sign * coef_buoyance * tem.value[IX3(i, j, k)] * mDt;
				vel.v[IX3(i, j - 1, k)] += 0.5f * coef_sign * coef_buoyance * tem.value[IX3(i, j, k)] * mDt;
				vel.v[IX3(i, j + 1, k)] += 0.5f * coef_sign * coef_buoyance * tem.value[IX3(i, j, k)] * mDt;
			}
}

void GRIDY_FIRE::Advect_temperature(GRIDY_SCALAR_FIELD& temperature, GRIDY_SCALAR_FIELD& temperature0, GRIDY_SCALAR_FIELD& fuel, GRIDY_VECTOR_FIELD& vel, double dt)
{
	mLOG.In("Advect temperature");

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

				temperature.value[IX3(i, j, k)] = r0*(s0*(t0*temperature0.value[IX3(i0, j0, k0)] + t1*temperature0.value[IX3(i0, j1, k0)]) + s1*(t0*temperature0.value[IX3(i1, j0, k0)] + t1*temperature.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*temperature0.value[IX3(i0, j0, k1)] + t1*temperature0.value[IX3(i0, j1, k1)]) + s1*(t0*temperature0.value[IX3(i1, j0, k1)] + t1*temperature0.value[IX3(i1, j1, k1)]));



			}

	Set_boundary(0, temperature.value);


	mLOG.Out();
}

void GRIDY_FIRE::Combustion(GRIDY_SCALAR_FIELD& fuel, GRIDY_SCALAR_FIELD& temperature, double ratio, double dt)
{
	double rest_ratio = (1.0f - ratio);
	double threshold_low_density = 0.01f;


	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				if (fuel.value[IX3(i, j, k)] >= threshold_low_density)	// 실행 속도를 위해 최소한의 값 이상에서만 연소 실행
				{
					temperature.value[IX3(i, j, k)] += fuel.value[IX3(i, j, k)] * ratio;
					fuel.value[IX3(i, j, k)] *= rest_ratio;
				}


				// 만약 최대값을 지정하고 싶다면
				if (temperature.value[IX3(i, j, k)] > 20.0f) temperature.value[IX3(i, j, k)] = 20.0f;
				if (fuel.value[IX3(i, j, k)] > 1.0f) fuel.value[IX3(i, j, k)] = 1.0f;
			}
}

void GRIDY_FIRE::Cooling(GRIDY_SCALAR_FIELD& temperature, double coef_cooling, double dt)
{
	double threshold_low_temperature = 0.0000001f;

	#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				// 쿨링 해야 함
				temperature.value[IX3(i, j, k)] = 1.0f / pow(((1 / pow(temperature.value[IX3(i, j, k)], 3)) + (3 * coef_cooling * dt)), 1.0f / 3.0f);

				// 너무 적은 온도는 날려버리자 (시뮬레이션 시간과 뷰어에서 그릴 때 시간도 줄어듦)
				if (temperature.value[IX3(i, j, k)] <= threshold_low_temperature)
					temperature.value[IX3(i, j, k)] = 0.0f;
			}
}

void GRIDY_FIRE::Simulate(int start_frame = 1, int end_frame = 300)
{
	mStart_frame = start_frame;
	mEnd_frame = end_frame;

	string prefix_frame = "FRAME #";
	char str_frame[10] = {};

	GRIDY_SIGNED_DISTANCE_FIELD mBunny;

	mLOG.In("SIMULATE");


	mLOG.In("Write Info File");
	mFILE_IO.Write_simulation_info(mSimulation_path, mWidth, mHeight, mDepth, mStart_frame, mEnd_frame, mData_file_type);
	mLOG.Out();



	//double mSource_x = (mWidth / 2.0f);
	//double mSource_y = (mHeight / 5.0f);

	mSrc_sphere.Init(1, mMaxMax_radius);
	mSrc_sphere.Set(0, 30.0f, 30.0f, 30.0f, 1.0f);
	
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++)
			{
				mSrc_sphere.isActivate[0] = true;
				mSrc_sphere.isGrowing[0] = true;
			}

	mSrc_sphere.mNo = 1;

	for (int frame = mStart_frame; frame <= mEnd_frame; frame++)
	{
		_itoa_s(frame, str_frame, 10);
		mLOG.In(prefix_frame + str_frame);

		mSrc_sphere.mMax_radius = mMaxMax_radius;

		Clear_previous_step_data();
		

		// 연료 공급
		double max_density = 20.0f;

		#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++)
				{
					// sphere source
					if (mSrc_sphere.Levelset(0, i, j, k) <= 0) mDensity0.value[IX3(i, j, k)] += 1.0f;
					if (mDensity0.value[IX3(i, j, k)] > max_density) mDensity0.value[IX3(i, j, k)] = max_density;
				}

		Velocity_step(mVelocity, mVelocity0, mTemperature, mDt);

		Density_step(mDensity, mDensity0, mTemperature, mTemperature0, mVelocity, mDiffuse, mDt);
		
		Combustion(mDensity, mTemperature, mCombustion_ratio, mDt);
		Cooling(mTemperature, mCoef_cooling, mDt);



		Write_data(frame);


		if (mSrc_sphere.isGrowing[0]==true)
			mSrc_sphere.Growing_radius(0, /*mEFP_growing_radius*/ 0.1f);


		mLOG.Out();
	}

	mLOG.Out();
}










#endif