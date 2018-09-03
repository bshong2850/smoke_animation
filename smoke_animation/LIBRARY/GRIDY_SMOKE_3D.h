/******************************************************************************
*	GRIDY_SMOKE_3D.h
*
*														Designed by Taehyeong
*
******************************************************************************/

#ifndef _GRIDY_SMOKE_3D_H_
#define _GRIDY_SMOKE_3D_H_





#include "GRIDY_SCALAR_FIELD.h"
#include "GRIDY_VECTOR_FIELD.h"
#include "GRIDY_FILE_IO.h"
#include "GRIDY_LOG.h"
#include "GRIDY_SIGNED_DISTANCE_FIELD.h"
//#include "GRIDY_IMPLICIT_OBJECT.h"




class GRIDY_SMOKE_3D : public GRIDY_SCALAR_FIELD
{
protected:
	int mStart_frame;
	int mEnd_frame;
	int mTotal_frame;

	double mDt, mDiffuse, mViscosity;

	int mWidth, mHeight, mDepth;
	string mSimulation_path;
	int mData_file_type;

	int mTargeting;

	double mCircle_x_dir, mCircle_y_dir;

	double mMaxMax_radius = 5.0f;



	GRIDY_LOG mLOG;
	GRIDY_FILE_IO mFILE_IO;
	GRIDY_COMMON mCOMMON;

	GRIDY_SIGNED_DISTANCE_FIELD mSDF, mBunny;

	GRIDY_SCALAR_FIELD density, density0;
	GRIDY_SCALAR_FIELD temperature, temperature0;

	GRIDY_VECTOR_FIELD velocity, velocity0;

	//GRIDY_SCALAR_FIELD phi;

	GRIDY_IMPLICIT_OBJECT_SPHERE mCircle, mSource;



public:
	GRIDY_SMOKE_3D(void);
	~GRIDY_SMOKE_3D(void);

	void Set_parameters(double dt, double diffuse, double viscosity, double force, double source);

	void Init(int width, int height, int depth, string path);
	void Init(int equally_size, string path);
	void Init(int width, int height, int depth, string path, string kinect_file);

	void Clear_previous_step_data(void);
	void Add_source(GRIDY_SCALAR_FIELD& target, GRIDY_SCALAR_FIELD& source, double dt);
	void Add_velocity(GRIDY_VECTOR_FIELD& target, GRIDY_VECTOR_FIELD& source, double dt);

	void Advect_velocity(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt);
	void Advect_velocity_for_ani(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt);
	void Advect_density(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt);
	void Advect_density_for_ani(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt);

	void Advect_velocity_with_obstacle(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, GRIDY_SIGNED_DISTANCE_FIELD& sdf, double dt);

	void Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p);

	void Diffuse(int b, double *x, double *x0, double dt);

	void Gause_seidel_solver(int b, double *x, double *x0, double a, double c);



	void Set_boundary(int b, double *x);

	void Velocity_step(void);
	void Density_step(GRIDY_SCALAR_FIELD& x, GRIDY_SCALAR_FIELD& x0, GRIDY_VECTOR_FIELD& vel, double diff, double dt);
	void Write_data(int frame);

	void Set_boundary_by_obstacle(double *flag, double *x, int direction);

	void Add_normal_force(GRIDY_VECTOR_FIELD &vel, double coef_normal_velocity, double dt);
	void Add_bouyance(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &dens, double coef_bouyance, double dt);
	void Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, double epsilon, double dt);


	void Diffuse_for_ani(int b, double *x, double *x0, double dt);
	void Density_step_for_ani(GRIDY_SCALAR_FIELD& x, GRIDY_SCALAR_FIELD& x0, GRIDY_VECTOR_FIELD& vel, double diff, double dt);


	void Simulate(int start_frame, int end_frame);


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







GRIDY_SMOKE_3D::GRIDY_SMOKE_3D(void)
{
	//	set parameters
	Set_parameters(0.1f, 0.0000005f, 0.0f, 5.0f, 100.0f);

	mStart_frame = mEnd_frame = 0;
	mWidth = mHeight = mDepth = 0;

	mSimulation_path = "C:\\";

	mData_file_type = mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY;
	//mData_file_type = mFILE_IO.GRIDY_DATA_FILE_TYPE_TEXT;
}




GRIDY_SMOKE_3D::~GRIDY_SMOKE_3D(void)
{
	//	terminate the program.
}





void GRIDY_SMOKE_3D::Set_parameters(double dt, double diffuse, double viscosity, double force, double source)
{
	mDt = dt;
	mDiffuse = diffuse;
	mViscosity = viscosity;
	//mForce = force;			// not use now
	//mSource = source;			// not use now
}





void GRIDY_SMOKE_3D::Init(int width, int height, int depth, string path = "C:\\")
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

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Init(int width, int height, int depth, string path, string kinect_file)
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

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Init(int equally_size, string path = "C:\\")
{
	Init(equally_size, equally_size, equally_size, path);
}





void GRIDY_SMOKE_3D::Clear_previous_step_data(void)
{
	mLOG.In("Clear data");

	density0.Clear_data();
	temperature0.Clear_data();
	velocity0.Clear_data();

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Add_source(GRIDY_SCALAR_FIELD& target, GRIDY_SCALAR_FIELD& source, double dt)
{
	mLOG.In("Add source");

	if (target.Get_size() != source.Get_size()) {
		cout << "\n\nERROR : Add_source\n\n";
	}

#pragma omp parallel for
	for (int i = 0; i<target.Get_size(); i++) {
		target.value[i] += dt * source.value[i];
	}

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Add_velocity(GRIDY_VECTOR_FIELD& target, GRIDY_VECTOR_FIELD& source, double dt)
{
	mLOG.In("Add velocity");

	if (target.Get_size() != source.Get_size()) {
		string err_msg = "ERROR : Add_source";
		mUTIL.Terminate_system(err_msg);
	}

#pragma omp parallel for
	for (int i = 0; i<target.Get_size(); i++) {
		target.u[i] += dt * source.u[i];
		target.v[i] += dt * source.v[i];
		target.w[i] += dt * source.w[i];
	}

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Advect_velocity(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt)
{
	mLOG.In("Advect velocity");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
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


void GRIDY_SMOKE_3D::Advect_velocity_for_ani(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, double dt)
{
	mLOG.In("Advect velocity");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
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



void GRIDY_SMOKE_3D::Advect_density_for_ani(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt)
{
	mLOG.In("Advect density for ani");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
				double x = i - dt0*vel.u[IX3(i, j, k)];
				double y = j - dt0*vel.v[IX3(i, j, k)];
				double z = k - dt0*vel.w[IX3(i, j, k)];

				if (x < 0.5)	x = 0.5f; if (x > mWidth + 0.5) x = mWidth + 0.5f;
				if (y < 0.5)	y = 0.5f; if (y > mHeight + 0.5) y = mHeight + 0.5f;
				if (z < 0.5)	z = 0.5f; if (z > mDepth + 0.5) z = mDepth + 0.5f;

				int i0 = (int)x;	int i1 = i0 + 1;
				int j0 = (int)y;	int j1 = j0 + 1;
				int k0 = (int)z;	int k1 = k0 + 1;

				double s1 = x - i0;	double s0 = 1 - s1;
				double t1 = y - j0;	double t0 = 1 - t1;
				double r1 = z - k0;	double r0 = 1 - r1;

				dens.value[IX3(i, j, k)] = r0*(s0*(t0*dens0.value[IX3(i0, j0, k0)] + t1*dens0.value[IX3(i0, j1, k0)]) + s1*(t0*dens0.value[IX3(i1, j0, k0)] + t1*dens0.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*dens0.value[IX3(i0, j0, k1)] + t1*dens0.value[IX3(i0, j1, k1)]) + s1*(t0*dens0.value[IX3(i1, j0, k1)] + t1*dens0.value[IX3(i1, j1, k1)]));
			}

	Set_boundary(0, dens.value);

	mLOG.Out();
}




void GRIDY_SMOKE_3D::Advect_density(GRIDY_SCALAR_FIELD& dens, GRIDY_SCALAR_FIELD& dens0, GRIDY_VECTOR_FIELD& vel, double dt)
{
	mLOG.In("Advect density");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
				double x = i - dt0*vel.u[IX3(i, j, k)];
				double y = j - dt0*vel.v[IX3(i, j, k)];
				double z = k - dt0*vel.w[IX3(i, j, k)];

				if (x < 0.5)	x = 0.5f; if (x > mWidth + 0.5) x = mWidth + 0.5f;
				if (y < 0.5)	y = 0.5f; if (y > mHeight + 0.5) y = mHeight + 0.5f;
				if (z < 0.5)	z = 0.5f; if (z > mDepth + 0.5) z = mDepth + 0.5f;

				int i0 = (int)x;	int i1 = i0 + 1;
				int j0 = (int)y;	int j1 = j0 + 1;
				int k0 = (int)z;	int k1 = k0 + 1;

				double s1 = x - i0;	double s0 = 1 - s1;
				double t1 = y - j0;	double t0 = 1 - t1;
				double r1 = z - k0;	double r0 = 1 - r1;

				dens.value[IX3(i, j, k)] = r0*(s0*(t0*dens0.value[IX3(i0, j0, k0)] + t1*dens0.value[IX3(i0, j1, k0)]) + s1*(t0*dens0.value[IX3(i1, j0, k0)] + t1*dens0.value[IX3(i1, j1, k0)])) + r1*(s0*(t0*dens0.value[IX3(i0, j0, k1)] + t1*dens0.value[IX3(i0, j1, k1)]) + s1*(t0*dens0.value[IX3(i1, j0, k1)] + t1*dens0.value[IX3(i1, j1, k1)]));
			}

	Set_boundary(0, dens.value);

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Gause_seidel_solver(int b, double *x, double *x0, double a, double c)
{
	mLOG.In("Execute Gause-seidel solver");

	for (int l = 0; l<20; l++) {
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++) {
					x[IX3(i, j, k)] = (x0[IX3(i, j, k)] + a*(x[IX3(i - 1, j, k)] + x[IX3(i + 1, j, k)] + x[IX3(i, j - 1, k)] + x[IX3(i, j + 1, k)] + x[IX3(i, j, k - 1)] + x[IX3(i, j, k + 1)])) / c;
				}

		Set_boundary(b, x);
	}

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Project(GRIDY_VECTOR_FIELD& vel, double *div, double *p)
{
	mLOG.In("Project");

	double h = (double)(1.0f / mWidth);
	if ((double)(1.0f / mHeight) > h) h = (double)(1.0f / mHeight);
	if ((double)(1.0f / mDepth) > h) h = (double)(1.0f / mDepth);


	//	compute divergence
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
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
			for (int k = 1; k <= mDepth; k++) {
				vel.u[IX3(i, j, k)] -= 0.5f *(p[IX3(i + 1, j, k)] - p[IX3(i - 1, j, k)]) / h;
				vel.v[IX3(i, j, k)] -= 0.5f *(p[IX3(i, j + 1, k)] - p[IX3(i, j - 1, k)]) / h;
				vel.w[IX3(i, j, k)] -= 0.5f *(p[IX3(i, j, k + 1)] - p[IX3(i, j, k - 1)]) / h;
			}

	Set_boundary(1, vel.u);
	Set_boundary(2, vel.v);
	Set_boundary(3, vel.w);

	mLOG.Out();

}



/*
void GRIDY_SMOKE_3D::Project_with_jump_condition(GRIDY_VECTOR_FIELD& vel, double *div, double *p)
{

mLOG.In("Project with jump condition");

double h = (double)(1.0f / mWidth);
if((double)(1.0f/mHeight) > h) h = (double)(1.0f / mHeight);
if((double)(1.0f/mDepth) > h) h = (double)(1.0f / mDepth);



#pragma omp parallel for
for(int i=1; i<=mWidth; i++)
for(int j=1; j<=mHeight; j++)
for(int k=1; k<=mDepth; k++)
{
div[IX3(i,j,k)] = -0.5f*h*(vel.u[IX3(i+1,j,k)] - vel.u[IX3(i-1,j,k)] + vel.v[IX3(i,j+1,k)] - vel.v[IX3(i,j-1,k)] + vel.w[IX3(i,j,k+1)] - vel.w[IX3(i,j,k-1)] );
p[IX3(i,j,k)] = 0;
}

Set_boundary(0, div);
Set_boundary(0, p);


Gause_seidel_solver(0, p, div, 1, 6);


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
*/




void GRIDY_SMOKE_3D::Velocity_step(void)
{
	mLOG.In("Velocity step");


	Add_velocity(velocity, velocity0, mDt);		// Add velocity0 to velocity

	Add_bouyance(velocity, density, 0.1f, mDt);

	SWAP_value(velocity, velocity0);

	//	Diffuse velocity (use only if we need)
	mLOG.In("Diffuse velocity");
	Diffuse(1, velocity.u, velocity0.u, mDt);
	Diffuse(2, velocity.v, velocity0.v, mDt);
	Diffuse(3, velocity.w, velocity0.w, mDt);
	mLOG.Out();



	//Add_normal_force(velocity0, 5.0f);



	Add_vorticity_confinement(velocity, 4.0f, mDt);

	Project(velocity, velocity0.u, velocity0.v);
	SWAP_value(velocity, velocity0);


	Advect_velocity(velocity, velocity0, mDt);
	Project(velocity, velocity0.u, velocity0.v);

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Diffuse(int b, double *x, double *x0, double dt)
{
	double a = mDt * mDiffuse * mWidth * mHeight * mDepth;

	Gause_seidel_solver(b, x, x0, a, 1 + 6 * a);
}





void GRIDY_SMOKE_3D::Density_step(GRIDY_SCALAR_FIELD& x, GRIDY_SCALAR_FIELD& x0, GRIDY_VECTOR_FIELD& vel, double diff, double dt)
{
	mLOG.In("Density step");

	Add_source(x, x0, dt);
	SWAP_value(x0, x);
	Diffuse(0, x.value, x0.value, dt);
	SWAP_value(x0, x);
	Advect_density(x, x0, vel, dt);

	mLOG.Out();
}



void GRIDY_SMOKE_3D::Diffuse_for_ani(int b, double *x, double *x0, double dt)
{
	double a = 10.5;
	//cout << "a : " << a << endl;

	Gause_seidel_solver(b, x, x0, a, 1 + 6 * a);
}



void GRIDY_SMOKE_3D::Density_step_for_ani(GRIDY_SCALAR_FIELD& x, GRIDY_SCALAR_FIELD& x0, GRIDY_VECTOR_FIELD& vel, double diff, double dt)
{
	mLOG.In("Density step");

	Add_source(x, x0, dt);
	SWAP_value(x0, x);
	Diffuse_for_ani(0, x.value, x0.value, dt);
	SWAP_value(x0, x);
	Advect_density(x, x0, vel, dt);

	mLOG.Out();
}




void GRIDY_SMOKE_3D::Set_boundary(int b, double *x)
{
	bool is_top_opened = true;



#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++) {
		for (int j = 1; j <= mDepth; j++) {
			x[IX3(i, 0, j)] = (b == 2 ? -x[IX3(i, 1, j)] : x[IX3(i, 1, j)]);

			if (is_top_opened) x[IX3(i, mHeight + 1, j)] = x[IX3(i, mHeight, j)] = 0;
			else x[IX3(i, mHeight + 1, j)] = (b == 2 ? -x[IX3(i, mHeight, j)] : x[IX3(i, mHeight, j)]);
		}
	}

#pragma omp parallel for
	for (int i = 1; i <= mHeight; i++) {
		for (int j = 1; j <= mDepth; j++) {
			x[IX3(0, i, j)] = (b == 1 ? -x[IX3(1, i, j)] : x[IX3(1, i, j)]);
			x[IX3(mWidth + 1, i, j)] = (b == 1 ? -x[IX3(mWidth, i, j)] : x[IX3(mWidth, i, j)]);
		}
	}

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++) {
		for (int j = 1; j <= mHeight; j++) {
			x[IX3(i, j, 0)] = (b == 3 ? -x[IX3(i, j, 1)] : x[IX3(i, j, 1)]);
			x[IX3(i, j, mDepth + 1)] = (b == 3 ? -x[IX3(i, j, mDepth)] : x[IX3(i, j, mDepth)]);
		}
	}

	/*
	//	set obstacle (ver 2)
	#pragma omp parallel for
	for(int i=1; i<=mWidth; i++)
	for(int j=1; j<=mHeight; j++)
	for(int k=1; k<=mDepth; k++)
	{
	if(mCircle.Levelset(i,j,k) <= 0)
	{
	switch(b)
	{
	case 1:
	if(mCircle.Levelset(i-1,j,k) >= 0) x[IX3(i,j,k)] = -x[IX3(i-1,j,k)];
	if(mCircle.Levelset(i+1,j,k) >= 0) x[IX3(i,j,k)] = -x[IX3(i+1,j,k)];
	break;

	case 2:
	if(mCircle.Levelset(i,j-1,k) >= 0) x[IX3(i,j,k)] = -x[IX3(i,j-1,k)];
	if(mCircle.Levelset(i,j+1,k) >= 0) x[IX3(i,j,k)] = -x[IX3(i,j+1,k)];
	break;

	case 3:
	if(mCircle.Levelset(i,j,k-1) >= 0) x[IX3(i,j,k)] = -x[IX3(i,j,k-1)];
	if(mCircle.Levelset(i,j,k+1) >= 0) x[IX3(i,j,k)] = -x[IX3(i,j,k+1)];
	break;

	case 0:
	int count = 0;
	double value_sum = 0.0f;

	if(mCircle.Levelset(i-1,j,k) >= 0) { value_sum += x[IX3(i-1,j,k)]; count++; }
	if(mCircle.Levelset(i+1,j,k) >= 0) { value_sum += x[IX3(i+1,j,k)]; count++; }
	if(mCircle.Levelset(i,j-1,k) >= 0) { value_sum += x[IX3(i,j-1,k)]; count++; }
	if(mCircle.Levelset(i,j+1,k) >= 0) { value_sum += x[IX3(i,j+1,k)]; count++; }
	if(mCircle.Levelset(i,j,k-1) >= 0) { value_sum += x[IX3(i,j,k-1)]; count++; }
	if(mCircle.Levelset(i,j,k+1) >= 0) { value_sum += x[IX3(i,j,k+1)]; count++; }

	if(count == 0)
	{
	//if(x[IX3(i,j,k)] > 0.00001f) cout << x[IX3(i,j,k)] << endl;
	x[IX3(i,j,k)] = 0;
	}
	else x[IX3(i,j,k)] = value_sum / count;

	break;
	}
	}
	}
	*/

	x[IX3(0, 0, 0)] = (x[IX3(1, 0, 0)] + x[IX3(0, 1, 0)] + x[IX3(0, 0, 1)]) / 3.0f;
	x[IX3(mWidth + 1, 0, 0)] = (x[IX3(mWidth, 0, 0)] + x[IX3(mWidth + 1, 1, 0)] + x[IX3(mWidth + 1, 0, 1)]) / 3.0f;
	x[IX3(0, 0, mDepth + 1)] = (x[IX3(1, 0, mDepth + 1)] + x[IX3(0, 1, mDepth + 1)] + x[IX3(0, 0, mDepth)]) / 3.0f;
	x[IX3(mWidth + 1, 0, mDepth + 1)] = (x[IX3(mWidth, 0, mDepth + 1)] + x[IX3(mWidth + 1, 1, mDepth + 1)] + x[IX3(mWidth + 1, 0, mDepth)]) / 3.0f;

	if (is_top_opened) {
		x[IX3(0, mHeight + 1, 0)] = x[IX3(0, mHeight + 1, mDepth + 1)] = x[IX3(mWidth + 1, mHeight + 1, 0)] = x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = 0;
	} else {
		x[IX3(0, mHeight + 1, 0)] = (x[IX3(1, mHeight + 1, 0)] + x[IX3(0, mHeight, 0)] + x[IX3(0, mHeight + 1, 1)]) / 3.0f;
		x[IX3(0, mHeight + 1, mDepth + 1)] = (x[IX3(1, mHeight + 1, mDepth + 1)] + x[IX3(0, mHeight, mDepth + 1)] + x[IX3(0, mHeight + 1, mDepth)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, 0)] = (x[IX3(mWidth, mHeight + 1, 0)] + x[IX3(mWidth + 1, mHeight, 0)] + x[IX3(mWidth + 1, mHeight + 1, 1)]) / 3.0f;
		x[IX3(mWidth + 1, mHeight + 1, mDepth + 1)] = (x[IX3(mWidth, mHeight + 1, mDepth + 1)] + x[IX3(mWidth + 1, mHeight, mDepth + 1)] + x[IX3(mWidth + 1, mHeight + 1, mDepth)]) / 3.0f;
	}

}




void GRIDY_SMOKE_3D::Write_data(int frame)
{
	mLOG.In("Write simulated data");


	mLOG.In("Writing Density");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) mFILE_IO.Write_data_binary(density.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	//else if(mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_TEXT) mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	else mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, mSimulation_path, frame, "density");
	mLOG.Out();

	mLOG.In("Writing Velocity");
	if (mData_file_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) {
		mFILE_IO.Write_data_binary(velocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data_binary(velocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data_binary(velocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	} else {
		mFILE_IO.Write_data(velocity.u, mWidth, mHeight, mDepth, mSimulation_path, frame, "u");
		mFILE_IO.Write_data(velocity.v, mWidth, mHeight, mDepth, mSimulation_path, frame, "v");
		mFILE_IO.Write_data(velocity.w, mWidth, mHeight, mDepth, mSimulation_path, frame, "w");
	}
	mLOG.Out();

	//mLOG.In("Writing Sphere Position Data");
	//	mFILE_IO.Write_sphere_data((int)mCircle.center_x, (int)mCircle.center_y, (int)mCircle.center_z, mCircle.radius, mSimulation_path, frame, "sphere");
	//mLOG.Out();


	mLOG.Out();
}





void GRIDY_SMOKE_3D::Advect_velocity_with_obstacle(GRIDY_VECTOR_FIELD& vel, GRIDY_VECTOR_FIELD& vel0, GRIDY_SIGNED_DISTANCE_FIELD& sdf, double dt)
{
	mLOG.In("Advect velocity with obstacle");

	double dt0 = dt * mWidth;
	if (dt * mHeight > dt0) dt0 = dt * mHeight;
	if (dt * mDepth > dt0) dt0 = dt * mDepth;

#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
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

				
				// obstacle control
				if((sdf.sign[IX3(i,j,k)] == sdf.GRIDY_SDF_INTERFACE && mTargeting == 1) || (mTargeting == 2 && mBunny.sign[IX3(i,j,k)] == mBunny.GRIDY_SDF_INTERFACE))
				{
				double coef_secondary_jump = 0.8f;

				vel0.u[IX3(i-1,j,k)] += vel.u[IX3(i,j,k)] * coef_secondary_jump;
				vel0.u[IX3(i+1,j,k)] += vel.u[IX3(i,j,k)] * coef_secondary_jump;
				vel0.v[IX3(i,j-1,k)] += vel.v[IX3(i,j,k)] * coef_secondary_jump;
				vel0.v[IX3(i,j+1,k)] += vel.v[IX3(i,j,k)] * coef_secondary_jump * 3.0f;
				vel0.w[IX3(i,j,k-1)] += vel.w[IX3(i,j,k)] * coef_secondary_jump;
				vel0.w[IX3(i,j,k+1)] += vel.w[IX3(i,j,k)] * coef_secondary_jump;

				//vel.u[IX3(i,j,k)] = -vel.u[IX3(i,j,k)];
				//vel.v[IX3(i,j,k)] = -vel.v[IX3(i,j,k)];
				//vel.w[IX3(i,j,k)] = -vel.w[IX3(i,j,k)];
				}
				
			}

	Set_boundary(1, vel.u);
	Set_boundary(2, vel.v);
	Set_boundary(3, vel.w);

	mLOG.Out();
}





void GRIDY_SMOKE_3D::Add_normal_force(GRIDY_VECTOR_FIELD &vel0, double coef_normal_velocity, double dt)
{
	mLOG.In("Add normal force");


#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
				if (mSource.Levelset(i, j, k) == 0) {

					//double coef_normal_velocity = 20.0f;
					double circle_normal[3];

					mSource.Get_normal(circle_normal, i, j, k);

					// Add normal velocity from circle
					vel0.u[IX3(i - 1, j, k)] += coef_normal_velocity * circle_normal[0] / 2.0f * mDt;
					vel0.v[IX3(i, j - 1, k)] += coef_normal_velocity * circle_normal[1] / 2.0f * mDt;
					vel0.w[IX3(i, j, k - 1)] += coef_normal_velocity * circle_normal[2] / 2.0f * mDt;

					vel0.u[IX3(i + 1, j, k)] += coef_normal_velocity * circle_normal[0] / 2.0f * mDt;
					vel0.v[IX3(i, j + 1, k)] += coef_normal_velocity * circle_normal[1] / 2.0f * mDt;
					vel0.w[IX3(i, j, k + 1)] += coef_normal_velocity * circle_normal[2] / 2.0f * mDt;

					/*
					vel0.u[IX3(i,j,k)] += coef_normal_velocity * circle_normal[0] * mDt;
					vel0.v[IX3(i,j,k)] += coef_normal_velocity * circle_normal[1] * mDt;
					vel0.w[IX3(i,j,k)] += coef_normal_velocity * circle_normal[2] * mDt;
					*/
				}
			}



	mLOG.Out();
}




void GRIDY_SMOKE_3D::Add_vorticity_confinement(GRIDY_VECTOR_FIELD &vel, double epsilon, double dt)
{
	mLOG.In("Add vorticity confinement");

	GRIDY_VECTOR_FIELD vortex(vel.Get_width_size(), vel.Get_height_size(), vel.Get_width_size());

	//	compute vortex
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
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
			for (int k = 1; k <= mDepth; k++) {

				double DoDx = (vortex.Get_norm2(i + 1, j, k) - vortex.Get_norm2(i - 1, j, k)) * 0.5f;
				double DoDy = (vortex.Get_norm2(i, j + 1, k) - vortex.Get_norm2(i, j - 1, k)) * 0.5f;
				double DoDz = (vortex.Get_norm2(i, j, k + 1) - vortex.Get_norm2(i, j, k - 1)) * 0.5f;


				double len = sqrt(DoDx*DoDx + DoDy*DoDy + DoDz*DoDz);

				if (len != 0) {
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




//	Add bouyance must be execute before projection
void GRIDY_SMOKE_3D::Add_bouyance(GRIDY_VECTOR_FIELD &vel, GRIDY_SCALAR_FIELD &dens, double coef_bouyance, double dt)
{
#pragma omp parallel for
	for (int i = 1; i <= mWidth; i++)
		for (int j = 1; j <= mHeight; j++)
			for (int k = 1; k <= mDepth; k++) {
				// is in inside?
				double coef_sign = 1.0f;


				// test for collision
				if (mSource.Levelset(i, j, k) <= 0) {
					vel.v[IX3(i, j - 1, k)] += coef_sign * coef_bouyance * dt * 0.5f;
					vel.v[IX3(i, j + 1, k)] += coef_sign * coef_bouyance * dt * 0.5f;
				} else {
					vel.v[IX3(i, j - 1, k)] += 0.5f * coef_sign * coef_bouyance * dens.value[IX3(i, j, k)] * mDt;
					vel.v[IX3(i, j + 1, k)] += 0.5f * coef_sign * coef_bouyance * dens.value[IX3(i, j, k)] * mDt;
				}



			}

}








void GRIDY_SMOKE_3D::Simulate(int start_frame = 1, int end_frame = 300)
{
	mStart_frame = start_frame;
	mEnd_frame = end_frame;

	string prefix_frame = "FRAME #";
	char str_frame[10] = {};

	GRIDY_SIGNED_DISTANCE_FIELD mBunny, mAngel;


	mLOG.In("SIMULATE");


	mLOG.In("Write Info File");
	mFILE_IO.Write_simulation_info(mSimulation_path, mWidth, mHeight, mDepth, mStart_frame, mEnd_frame, mData_file_type);
	mLOG.Out();





	double mSource_x = (mWidth / 2.0f);




	for (int frame = mStart_frame; frame <= mEnd_frame; frame++) {
		_itoa_s(frame, str_frame, 10);
		mLOG.In(prefix_frame + str_frame);


		Clear_previous_step_data();
		double xx;
		double yy;
		double zz;

		if (frame <= 50) {
			xx = mWidth / 2.0 + (50 - frame) / 50.0 * mWidth / 3.0;
			yy = mHeight / 15.0 + (1 - pow((frame) / 50.0, 2)) * mHeight / 2.0;
			zz = mDepth / 2.0;
		} else if (frame <= 80) {
			xx = mWidth / 2.0 + (50 - frame) / 50.0 * mWidth / 3.0;
			yy = mHeight / 15.0 + (1 - pow((50 - frame + 15) / 15.0, 2)) * mHeight / 8.0;
			zz = mDepth / 2.0;
		} else if (frame <= 90) {
			xx = mWidth / 2.0 + (50 - frame) / 50.0 * mWidth / 3.0;
			yy = mHeight / 15.0 + (1 - pow((80 - frame + 5) / 5.0, 2)) * mHeight / 24.0;
			zz = mDepth / 2.0;
		} else {
			xx = mWidth / 2.0;
			yy = mHeight / 15.0;
			zz = mDepth / 2.0;

			xx -= mWidth * 4.0 / 15.0 * cos((frame - 90) / 20.0);
			zz += mWidth * 4.0 / 15.0 * sin((frame - 90) / 20.0);
		}


		xx = mWidth / 2.0;
		yy = mHeight / 15.0;
		zz = mDepth / 2.0;

		mSource.Init(xx, yy, zz, mWidth / 5.0f);


#pragma omp parallel for
		for (int i = 1; i <= mWidth; i++)
			for (int j = 1; j <= mHeight; j++)
				for (int k = 1; k <= mDepth; k++) {
					if (mSource.Levelset(i, j, k) < 0) {
						density0.value[IX3(i, j, k)] += (double)(mWidth / 32.0f);
						velocity0.v[IX3(i, j, k)] += (double)(mHeight / 128.0f);
					}
				}

		Velocity_step();

		Density_step(density, density0, velocity, mDiffuse, mDt);

		Write_data(frame);


		mLOG.Out();
	}

	mLOG.Out();
}










#endif