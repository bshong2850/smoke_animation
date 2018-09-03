/******************************************************************************
*	GRIDY_GL_VIEWER v0.6.0 (with Shader)
*
*												Designed by TaeHyeong, EunKi
*
*
******************************************************************************/



//#pragma comment(lib, "glew32.lib")
//#pragma comment(lib, "FreeImage.lib")

//#include "../LIBRARY/GRIDY_SHADER_FOR_GL_VIEWER.h"

//#define STB_IMAGE_IMPLEMENTATION
//#include "../../LIBRARY/stb_image.h"

//#include <GL/glut.h>
//#include <GL/gl.h>
//#include <GL/glu.h>



//#include "../../LIBRARY/GRIDY_COMMON.h"
#include "../LIBRARY/GRIDY_SCALAR_FIELD.h"
#include "../LIBRARY/GRIDY_VECTOR_FIELD.h"
#include "../LIBRARY/GRIDY_FILE_IO.h"
#include "../LIBRARY/GRIDY_IMPLICIT_OBJECTS.h"
//#include "../../LIBRARY/GRIDY_PBD.h"


#include "../LIBRARY/GRIDY_CONFIG.h"

#include <iostream>
//#include <FreeImage.h>
#include <GL/glut.h>

using namespace std;

//#define IX3(i,j,k)	((i)+(mWidth+2)*(j)+(mWidth+2)*(mHeight+2)*(k))

static int mWin_id;
static int mWin_width, mWin_height;
int mMouse_down[3];
int mMouse_x, mMouse_y, mMouse_x0, mMouse_y0;

GLfloat mZoom;
GLfloat mRotateX, mRotateY, mRotateZ, mAngle;
GLfloat mTransX, mTransY, mTransZ, mTransSize;
GLfloat mScaleX, mScaleY, mScaleZ;
GLint mLastX, mLastY;
GLfloat mLookAt;
GLdouble mFov, mAspect, mNear, mFar;

float mSlide;

int mDimension;
int mFrame_no;
int mFrame_start, mFrame_end;
int mWidth, mHeight, mDepth;
int mFile_type;

int mPsize;
int mFsize;

string full_file_path;

int option_show_velocity = 0;
int option_show_levelset = 0;
int option_show_temperature = 0;
int option_show_phi = 0;
int option_pause = 0;
int option_slide = 0;
int option_sphere = 0;
int option_convert = 0;
int option_pbd = 0;
bool option_show_particle = false;
bool option_show_particle_only = false;

string mView_mode;

//float option_range = 0.001f;
//float option_range = 0.225f; // 기본값
float option_range = 0.115f; // 기본값

GRIDY_SCALAR_FIELD density, temperature;
GRIDY_VECTOR_FIELD velocity;
GRIDY_IMPLICIT_OBJECT_SPHERES fuel_particle;
GRIDY_FILE_IO mFILE_IO;

//PBD pbd; 

int mSphere_x, mSphere_y, mSphere_z; float mSphere_r;

void Init_data(void);
void delete_data(void);

//GRIDY_SHADER_FOR_GL_VIEWER *shader;
static GLuint flake_texture;
/*
static void Flake_texture_init()
{
	glGenTextures(1, &flake_texture);
	glBindTexture(GL_TEXTURE_2D, flake_texture);

//	FILE *fp = fopen("F:\\GRIDY_RESULT\\MODEL\\fireflake_by_inmc.png", "rb");
//	FILE *fp = fopen("F:\\GRIDY_RESULT\\MODEL\\flower2.png", "rb");

	int w, h, c;
	auto rd = stbi_load_from_file(fp, &w, &h, &c, 4);


	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, rd);
	glBindTexture(GL_TEXTURE_2D, 0);

	free(rd);

	fclose(fp);
}
*/
//static void Shader_init()
//{
//	GLenum err = glewInit();
//
//	if (GLEW_OK != err)
//	{
//		cout << "hull\n" << endl;
//		cout << glewGetErrorString(err) << endl;
//	}
//	shader = new GRIDY_SHADER_FOR_GL_VIEWER();
//	shader->InitShader("./GRIDY_shader_v.glsl", "./GRIDY_shader_f.glsl");
//	glBindAttribLocation(shader->program, 1, "Radius");
//	glBindAttribLocation(shader->program, 7, "Angle");
//	glBindAttribLocation(shader->program, 8, "Temperature");
//	shader->Link();
//
//	//Flake_texture_init();
//
//	// create and set-up the vertex array object
//	glGenVertexArrays(1, &shader->vaoHandle);
//	glBindVertexArray(shader->vaoHandle);
//}
/*
----------------------------------------------------------------------
OpenGL specific drawing routines
----------------------------------------------------------------------
*/

void Pre_display(void)
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	glColor3f(0.8, 0.8, 0.8);
	glLoadIdentity();


	gluLookAt(0.0, mLookAt, 15.0, 0.0, 0.0, -10.0, 0.0, 1.0, 0.0);

	glTranslatef(0, 0, -mZoom);
	glTranslatef(mTransX, mTransY, mTransZ);
	glRotatef(mRotateX, 1.0, 0.0, 0.0);
	glRotatef(mRotateY, 0.0, 1.0, 0.0);
	glRotatef(mRotateZ, 0.0, 0.0, 1.0);
	glScalef(mScaleX, mScaleY, mScaleZ);


	glViewport(0, 0, mWin_width, mWin_height);
	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	gluPerspective(mAngle, mAspect, 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
}

//static int capture_frame = 0;
//static int doo = 900;
//
//
//void capture()
//{
//	if (doo > 0){ doo = 0; }
//	else{ doo++; return; }
//	if (capture_frame >= 250) exit(0);
//	// Make the BYTE array, factor of 3 because it's RBG.
//	BYTE* pixels = new BYTE[3 * mWin_width * mWin_height];
//
//	glReadPixels(0, 0, mWin_width, mWin_height, GL_BGR, GL_UNSIGNED_BYTE, pixels);
//
//	// Convert to FreeImage format & save to file
//	FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, mWin_width, mWin_height, 3 * mWin_width, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
//	char filename[200];
//	sprintf(filename, "E:\\GRIDY\\RESULT\\Captured\\%05d.png", capture_frame);
//	capture_frame++;
//	FreeImage_Save(FIF_PNG, image, filename, 0);
//
//	// Free resources
//	FreeImage_Unload(image);
//	delete[] pixels;
//}

static void Post_display(void)
{
	//capture();
	glutSwapBuffers();
} 

static void Draw_velocity_3D(void)
{
	int i, j, k;
	float x, y, z, h;

	int N = mDepth;
	if (mHeight > N) N = mHeight;
	else if (mWidth > N) N = mWidth;

	h = 10.0f / N;

	float color_map[5][3] = { { 0, 0, 1 }, { 0, 1, 1 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 } };
	int c[5] = { 0, 0, 0, 0, 0 };

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);


	for (i = 1; i <= mWidth; i++) {
		x = (i - mWidth / 2 - 0.5f)*h;
		for (j = 1; j <= mHeight; j++) {
			y = (j - mHeight / 2 - 0.5f)*h;
			for (k = 1; k <= mDepth; k++) {
				//z = (mDepth/2-k-0.5f)*h;
				z = (k - mDepth / 2 - 0.5f)*h;


				glColor3f(velocity.u[IX3(i, j, k)], velocity.v[IX3(i, j, k)], velocity.w[IX3(i, j, k)]);

				float temp_color = (abs(velocity.u[IX3(i, j, k)]) + abs(velocity.v[IX3(i, j, k)]) + abs(velocity.w[IX3(i, j, k)]));
				if (temp_color <= 0.0001f) continue;
				else if (temp_color > 10.0f) temp_color = 10.0f;
				temp_color *= 10.0f;

				int color_range_low = (int)floor(temp_color);
				int color_range_high = (int)ceil(temp_color);
				int color_range_value = temp_color - color_range_low;


				glColor3f((color_map[color_range_low][0])*(1 - color_range_value) + (color_map[color_range_high][0])*color_range_value,
					(color_map[color_range_low][1])*(1 - color_range_value) + (color_map[color_range_high][1])*color_range_value,
					(color_map[color_range_low][2])*(1 - color_range_value) + (color_map[color_range_high][2])*color_range_value);

				glVertex3f(x, y, z);
				glVertex3f(x + (velocity.u[IX3(i, j, k)])*h*5.0f, y + (velocity.v[IX3(i, j, k)])*h*5.0f, z + (velocity.w[IX3(i, j, k)])*h*5.0f);

				/*
				if(velocity.v[IX3(i,j,k)] < 0.05f && velocity.v[IX3(i,j,k)] > -0.05f) return;
				if(velocity.v[IX3(i,j,k)] > 0 )glColor3f(1,0,1);
				else glColor3f(0,1,1);

				glVertex3f ( x, y, z );
				glVertex3f ( x, y+(velocity.v[IX3(i,j,k)])*h*5.0f, z );
				*/
			}
		}
	}

	glEnd();
}

static void Draw_base_line(void)
{
	// Draw Bottom Grids and Initial items
	glBegin(GL_LINES);
	
	// Draw bottom grids
	glColor4f(0.2f, 0.2f, 0.2f, 0.3f);

	for (int i = -10; i <= 10; i++)
	{
		glVertex3f(-10, -5, i);
		glVertex3f(10, -5, i);

		glVertex3f(i, -5, -10);
		glVertex3f(i, -5, 10);
	}
	glEnd();


	// Draw 0,0 block (blue one)
	glBegin(GL_QUADS);
	glColor4f(0.2f, 0.2f, 0.4f, 0.3f);
	glVertex3f(-10, -5, 10); glVertex3f(-10, -5, 9); glVertex3f(-9, -5, 9); glVertex3f(-9, -5, 10);
	glEnd();

	/*
	// Bottom Grid
	glColor4f(0.95f, 0.95f, 0.95f, 0.5f);
	glBegin(GL_LINE_STRIP);
	glVertex3f(0, 0, 0); glVertex3f(2, 0, 0); glVertex3f(2, 0, 2); glVertex3f(0, 0, 2); glVertex3f(0, 0, 0);
	glEnd();
	*/

	int N = mDepth;
	if (mHeight > N) N = mHeight;
	else if (mWidth > N) N = mWidth;

	float x1, x2, y1, y2, z1, z2;
	float h = 10.0f / N;

	x1 = (1 - mWidth / 2 - 0.5f)*h;
	x2 = (mWidth - mWidth / 2 - 0.5f)*h;

	y1 = (1 - mHeight / 2 - 0.5f)*h;
	y2 = (mHeight - mHeight / 2 - 0.5f)*h;

	z1 = (mDepth / 2 - 1 - 0.5f)*h;
	z2 = (mDepth / 2 - mDepth - 0.5f)*h;

	// Draw GRID boundary (4x4x4)
	//glBegin(GL_QUADS);
	glBegin(GL_LINE_STRIP);
	glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
	glVertex3f(x1, y2, z2); glVertex3f(x2, y2, z2); glVertex3f(x2, y2, z1); glVertex3f(x1, y2, z1); glVertex3f(x1, y2, z2);
	glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
	glVertex3f(x1, y1, z2); glVertex3f(x2, y1, z2); glVertex3f(x2, y1, z1); glVertex3f(x1, y1, z1); glVertex3f(x1, y1, z2);
	glEnd();

}

static void Draw_grid_line(void)
{
	float x, y, z;

	int N = mDepth;
	if (mHeight > N) N = mHeight;
	else if (mWidth > N) N = mWidth;

	float h = 10.0f / N;

	for (int i = 1; i <= mWidth; i++)
	{
		x = (i - mWidth / 2 - 0.5f)*h;

		for (int j = 1; j <= mHeight; j++)
		{
			y = (j - mHeight / 2 - 0.5f)*h;

			for (int k = 1; k <= mDepth; k++)
			{
				glColor3f(0.3f, 0.3f, 0.3f);
				//	z = (k-mDepth/2-0.5f)*h;
				z = (mDepth / 2 - k + 0.5f)*h;

				glBegin(GL_LINE_STRIP);
				//	Front : 000-100-110-010
				glVertex3f(x, y, z);
				glVertex3f(x + h, y, z);
				glVertex3f(x + h, y + h, z);
				glVertex3f(x, y + h, z);
				glVertex3f(x, y, z);
				glEnd();

				glBegin(GL_LINE_STRIP);
				//	Back : 001-101-111-011
  				glVertex3f(x, y, z - h);
				glVertex3f(x + h, y, z - h);
				glVertex3f(x + h, y + h, z - h);
				glVertex3f(x, y + h, z - h);
				glVertex3f(x, y + h, z - h);
				glEnd();

				glBegin(GL_LINES);
				glVertex3f(x, y, z);		glVertex3f(x, y, z - h);
				glVertex3f(x + h, y, z);	glVertex3f(x + h, y, z - h);
				glVertex3f(x + h, y + h, z);	glVertex3f(x + h, y + h, z - h);
				glVertex3f(x, y + h, z);	glVertex3f(x, y + h, z - h);
				glEnd();
			}
		}
	}


}

static void Draw_density_3D(void)
{
	int i, j, k;
	float x, y, z, h; //d000, d100, d101, d001, d010, d110, d111, d011;

	int N = (float)mDepth;
	if (mHeight > N) N = mHeight;
	else if (mWidth > N) N = mWidth;

	h = 10.0f / N;
	//if(mWin_width > mWin_height) h = (mWin_height * 0.1f) / (float)N;
	//else h = (mWin_width * 0.1f) / (float)N;

	//glUseProgram(shader->program);
	glBegin(GL_QUADS);
	//for (i = 1; i <= mWidth; i++)
	//{
	//	x = (i - mWidth - 0.5f)*h;

	//	for (j = 1; j <= mHeight; j++)
	//	{
	//		y = (j - mHeight - 0.5f)*h;

	//		for (k = 1; k <= mDepth; k++)
	//		{
	//			z = (k - mDepth - 0.5f)*h;
	//			//	z = (mDepth/2-k-0.5f)*h;

	//			if (option_show_temperature)
	//			{
	//				//if (temperature.value[IX3(i, j, k)] > 0.001f)
	//				if (temperature.value[IX3(i, j, k)] > option_range)
	//				{
	//					glColor4f(2.299*temperature.value[IX3(i, j, k)], .587*temperature.value[IX3(i, j, k)], .114*temperature.value[IX3(i, j, k)], 1.0);
	//				}
	//				else continue;
	//			}
	//			else
	//			{
	//				//if (density.value[IX3(i, j, k)] > 0.001f)
	//				if (density.value[IX3(i, j, k)] > option_range)
	//					glColor4f(density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)]);
	//				else continue;
	//			}

	//			double c = density.value[IX3(i, j, k)] * 2;
	//			if (option_show_temperature)
	//				c = temperature.value[IX3(i, j, k)] * 2;
	//			//	glColor4f(c*c*c*c, c*c, c, c);
	//			glColor4f(c, c, c, pow(c, 2));
	//			//glColor4f(c, c*c, c*c*c*c, pow(c, 2));

	//			// Draw two-planes
	//			//glVertex3f(x, y, z);			glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z + h);		glVertex3f(x, y + h, z + h);
	//			//glVertex3f(x + h, y, z);		glVertex3f(x, y + h, z);		glVertex3f(x, y + h, z + h);			glVertex3f(x + h, y + h, z);
	//			//glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );			glVertex3f ( x, y+h, z+h );


	//			// Draw cube
	//			//	Front : 000-100-110-010
	//			glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y + h, z);			glVertex3f(x, y + h, z);
	//			//	Left : 000-001-011-010
	//			glVertex3f(x, y, z);			glVertex3f(x, y, z - h);		glVertex3f(x, y + h, z - h);			glVertex3f(x, y + h, z);
	//			//	Right : 100-101-111-110
	//			glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x + h, y + h, z);
	//			//	Back : 001-101-111-011
	//			glVertex3f(x, y, z - h);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
	//			//	Top : 010-110-111-011
	//			glVertex3f(x, y + h, z);		glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
	//			//	Bottom : 000-100-101-001
	//			glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);			glVertex3f(x, y, z - h);

	//		}
	//	}
	//}
	if (option_slide)
	{
		for (i = 1; i <= mWidth; i++)
		{
			x = (i - mWidth / 2 - 0.5f)*h;

			for (j = 1; j <= mHeight; j++)
			{
				y = (j - mHeight / 2 - 0.5f)*h;

				k = mSlide;

				z = (k - mDepth / 2 - 0.5f)*h;
				//	z = (mDepth/2-k-0.5f)*h;
				float c = density.value[IX3(i, j, k)];
				//if (density.value[IX3(i, j, k)] < 0.5) glColor4f(0, density.value[IX3(i, j, k)], 0, density.value[IX3(i, j, k)]);
				glColor4f(c*c*c*c, c*c, c, c);
				
				
				glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );		glVertex3f ( x, y+h, z+h );
				glVertex3f ( x+h, y, z );		glVertex3f ( x, y+h, z );		glVertex3f ( x, y+h, z+h );			glVertex3f ( x+h, y+h, z );
				glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );		glVertex3f ( x, y+h, z+h );
				

				/*
				//	Front : 000-100-110-010
				glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y + h, z);			glVertex3f(x, y + h, z);
				//	Left : 000-001-011-010
				glVertex3f(x, y, z);			glVertex3f(x, y, z - h);		glVertex3f(x, y + h, z - h);			glVertex3f(x, y + h, z);
				//	Right : 100-101-111-110
				glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x + h, y + h, z);
				//	Back : 001-101-111-011
				glVertex3f(x, y, z - h);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
				//	Top : 010-110-111-011
				glVertex3f(x, y + h, z);		glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
				//	Bottom : 000-100-101-001
				glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);			glVertex3f(x, y, z - h);
				*/
			}
		}
	}
	else
	{
		for (i = 1; i <= mWidth; i++)
		{
			x = (i - mWidth / 2 - 0.5f)*h;

			for (j = 1; j <= mHeight; j++)
			{
				y = (j - mHeight / 2 - 0.5f)*h;

				for (k = 1; k <= mDepth; k++)
				{
					z = (k - mDepth / 2 - 0.5f)*h;
					
					if (density.value[IX3(i, j, k)] > option_range)
						glColor4f(density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)]);
					else continue;
					

					double c = density.value[IX3(i, j, k)] * 2;
					if (option_show_temperature)
						c = temperature.value[IX3(i, j, k)] * 2;
				//	glColor4f(c*c*c*c, c*c, c, c);
					glColor4f(c, c, c, pow(c, 2));
					//glColor4f(c, c*c, c*c*c*c, pow(c, 2));
					
					// Draw two-planes
					/*glVertex3f(x, y, z);			glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z + h);		glVertex3f(x, y + h, z + h);
					glVertex3f(x + h, y, z);		glVertex3f(x, y + h, z);		glVertex3f(x, y + h, z + h);			glVertex3f(x + h, y + h, z);*/
					//glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );			glVertex3f ( x, y+h, z+h );
					
					
					// Draw cube
					//	Front : 000-100-110-010
					glVertex3f ( x, y, z );			glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y+h, z );			glVertex3f ( x, y+h, z );
					//	Left : 000-001-011-010
					glVertex3f ( x, y, z );			glVertex3f ( x, y, z-h );		glVertex3f ( x, y+h, z-h );			glVertex3f ( x, y+h, z );
					//	Right : 100-101-111-110
					glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y, z-h );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x+h, y+h, z );
					//	Back : 001-101-111-011
					glVertex3f ( x, y, z-h );		glVertex3f ( x+h, y, z-h );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x, y+h, z-h );
					//	Top : 010-110-111-011
					glVertex3f ( x, y+h, z );		glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x, y+h, z-h );
					//	Bottom : 000-100-101-001
					glVertex3f ( x, y, z );			glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y, z-h );			glVertex3f ( x, y, z-h );
					
				}
			}
		}
		//for (i = 1; i <= mWidth / 2; i++)
		//{
		//	x = (i - mWidth / 2 - 0.5f)*h;
		//	for (j = 1; j <= mHeight; j++)
		//	{
		//		y = (j - mHeight / 2 - 0.5f)*h;
		//		for (k = 1; k <= mDepth; k++)
		//		{
		//			z = (k - mDepth / 2 - 0.5f)*h;
		//			//	z = (mDepth/2-k-0.5f)*h;
		//			if (density.value[IX3(i, j, k)] > option_range)
		//				glColor4f(density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)]);
		//			else continue;
		//			float c = density.value[IX3(i, j, k)] * 2;
		//			if (option_show_temperature)
		//				c = temperature.value[IX3(i, j, k)] * 2;
		//			glColor4f(c, c, c, pow(c, 2)); 
		//			//glColor4f(c, c*c, c*c*c*c, pow(c, 2));
		//			// Draw two-planes
		//			/*glVertex3f(x, y, z);			glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z + h);		glVertex3f(x, y + h, z + h);
		//			glVertex3f(x + h, y, z);		glVertex3f(x, y + h, z);		glVertex3f(x, y + h, z + h);			glVertex3f(x + h, y + h, z);*/
		//			//glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );			glVertex3f ( x, y+h, z+h );
		//			
		//			// Draw cube
		//			//	Front : 000-100-110-010
		//			glVertex3f ( x, y, z );			glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y+h, z );			glVertex3f ( x, y+h, z );
		//			//	Left : 000-001-011-010
		//			glVertex3f ( x, y, z );			glVertex3f ( x, y, z-h );		glVertex3f ( x, y+h, z-h );			glVertex3f ( x, y+h, z );
		//			//	Right : 100-101-111-110
		//			glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y, z-h );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x+h, y+h, z );
		//			//	Back : 001-101-111-011
		//			glVertex3f ( x, y, z-h );		glVertex3f ( x+h, y, z-h );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x, y+h, z-h );
		//			//	Top : 010-110-111-011
		//			glVertex3f ( x, y+h, z );		glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z-h );		glVertex3f ( x, y+h, z-h );
		//			//	Bottom : 000-100-101-001
		//			glVertex3f ( x, y, z );			glVertex3f ( x+h, y, z );		glVertex3f ( x+h, y, z-h );			glVertex3f ( x, y, z-h );
		//			
		//		}
		//	}
		//}
		//for (i = mWidth; i > mWidth / 2; i--)
		//{
		//	x = (i - mWidth / 2 - 0.5f)*h;
		//	for (j = 1; j <= mHeight; j++)
		//	{
		//		y = (j - mHeight / 2 - 0.5f)*h;
		//		for (k = 1; k <= mDepth; k++)
		//		{
		//			z = (k - mDepth / 2 - 0.5f)*h;
		//			//	z = (mDepth/2-k-0.5f)*h;
		//			if (density.value[IX3(i, j, k)] > option_range)
		//				glColor4f(density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)], density.value[IX3(i, j, k)]);
		//			else continue;
		//			float c = density.value[IX3(i, j, k)] * 2;
		//			if (option_show_temperature)
		//				c = temperature.value[IX3(i, j, k)] * 2;
		//			glColor4f(c, c, c, pow(c, 2));
		//			//glColor4f(c, c*c, c*c*c*c, pow(c, 2));
		//			// Draw two-planes
		//			/*glVertex3f(x, y, z);			glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z + h);		glVertex3f(x, y + h, z + h);
		//			glVertex3f(x + h, y, z);		glVertex3f(x, y + h, z);		glVertex3f(x, y + h, z + h);			glVertex3f(x + h, y + h, z);*/
		//			//glVertex3f ( x, y, z );			glVertex3f ( x+h, y+h, z );		glVertex3f ( x+h, y+h, z+h );			glVertex3f ( x, y+h, z+h );

		//			// Draw cube
		//			//	Front : 000-100-110-010
		//			glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y + h, z);			glVertex3f(x, y + h, z);
		//			//	Left : 000-001-011-010
		//			glVertex3f(x, y, z);			glVertex3f(x, y, z - h);		glVertex3f(x, y + h, z - h);			glVertex3f(x, y + h, z);
		//			//	Right : 100-101-111-110
		//			glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x + h, y + h, z);
		//			//	Back : 001-101-111-011
		//			glVertex3f(x, y, z - h);		glVertex3f(x + h, y, z - h);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
		//			//	Top : 010-110-111-011
		//			glVertex3f(x, y + h, z);		glVertex3f(x + h, y + h, z);		glVertex3f(x + h, y + h, z - h);		glVertex3f(x, y + h, z - h);
		//			//	Bottom : 000-100-101-001
		//			glVertex3f(x, y, z);			glVertex3f(x + h, y, z);		glVertex3f(x + h, y, z - h);			glVertex3f(x, y, z - h);

		//		}
		//	}
		//}
	}
	glEnd();
	//glUseProgram(0);
}

static void Light_init()
{

	glEnable(GL_DEPTH_TEST);


	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);

	GLfloat LightPosition[4] = { 0, 100, 0, 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

	GLfloat LightAmbient[4] = { 0.2, 0.2, 0.2, 0.2 };
	GLfloat LightDiffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat LightSpecular[4] = { 1.0, 1.0, 1.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
}

void Set_viewer_attr(void)
{
	mZoom = 5.0f;
	mRotateX = 0, mRotateY = 0, mRotateZ = 0, mAngle = 45.0;
	mTransX = 0, mTransY = 0, mTransZ = 0, mTransSize = 0.5;
	mScaleX = 1.0, mScaleY = 1.0, mScaleZ = 1.0;
	mLastX = 0, mLastY = 0;
	mLookAt = 0.0;

	mFov = 45.0, mAspect = 1.0, mNear = -2.0, mFar = 100.0;

	mMouse_down[0] = 0;
	mMouse_down[1] = 0;
	mMouse_down[2] = 0;

	/*
	GLfloat lightZeroPosition[] = {-150.0, -150.0, 0.0, 50.0};
	GLfloat lightZeroColor[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lightOnePosition[] = {-50.0, -50.0, -50.0, 50.0};
	GLfloat lightOneColor[] = {1.0, 1.0, 1.0, 1.0};
	*/

	// printf("Viewer Attribute Set complete.\n");

	mView_mode = "";
}

static void Read_data(int frame_no)
{
	if (mFile_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) 	// 1: binary, 0: text
	{
		if (option_show_velocity)
		{
			mFILE_IO.Read_data_binary(velocity.u, full_file_path, frame_no, "u");
			mFILE_IO.Read_data_binary(velocity.v, full_file_path, frame_no, "v");
			mFILE_IO.Read_data_binary(velocity.w, full_file_path, frame_no, "w");
		}
		else if (option_show_temperature) mFILE_IO.Read_data_binary(temperature.value, full_file_path, frame_no, "temperature");
		else mFILE_IO.Read_data_binary(density.value, full_file_path, frame_no, "density");

		if (option_show_particle || option_show_particle_only)
		{
			mFILE_IO.Read_particle_data(fuel_particle, full_file_path, frame_no, "particle");
		}
	}
	else
	{
		if (option_show_velocity)
		{
			mFILE_IO.Read_data(velocity.u, full_file_path, frame_no, "u");
			mFILE_IO.Read_data(velocity.v, full_file_path, frame_no, "v");
			mFILE_IO.Read_data(velocity.w, full_file_path, frame_no, "w");
		}
		else if (option_show_temperature) mFILE_IO.Read_data_binary(temperature.value, full_file_path, frame_no, "temperature");
		else mFILE_IO.Read_data(density.value, full_file_path, frame_no, "density");
	}

	/*
	if(option_pbd)
	{
	mFILE_IO.Read_data_binary(pbd.obj->facePositionData, full_file_path, mFsize*12, frame_no, "f_pos");
	mFILE_IO.Read_data_binary(pbd.obj->faceNormalData, full_file_path, mFsize*12, frame_no, "f_norm");
	}

	if(option_sphere)
	{
	mFILE_IO.Read_sphere_data(mSphere_x, mSphere_y, mSphere_z, mSphere_r, full_file_path, frame_no, "sphere");
	}
	*/

	//printf("# %5d frame.\n", frame_no);
}

static void Convert_file_type(int frame_no)
{
	string converted_file_path = full_file_path + "\\GRIDY_FILE_CONVERTED";

	if (option_convert == 2)		// 현재 frame만 convert
	{
		if (mFile_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) 	// 1: binary, 0: text
		{
			
			mFILE_IO.Read_data_binary(density.value, full_file_path, frame_no, "density");
			mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, converted_file_path, frame_no, "density");
			/*
			mFILE_IO.Read_data_binary(density.value, full_file_path, frame_no, "temperature");
			mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, converted_file_path, frame_no, "temperature");*/
		}
		else
		{
			mFILE_IO.Read_data(density.value, full_file_path, frame_no, "density");
			mFILE_IO.Write_data_binary(density.value, mWidth, mHeight, mDepth, converted_file_path, frame_no, "density");
		}
		/*
		mFILE_IO.Read_data_binary(pbd.obj->facePositionData, full_file_path, mFsize*12, frame_no, "f_pos");
		mFILE_IO.Read_data_binary(pbd.obj->faceNormalData, full_file_path, mFsize*12, frame_no, "f_norm");
		printf("# %5d frame data converted.\n", frame_no);
		*/
	}
	else	// 모든 파일 convert
	{
		for (int frame = mFrame_start; frame <= mFrame_end; frame++)
		{
			if (mFile_type == mFILE_IO.GRIDY_DATA_FILE_TYPE_BINARY) 	// 1: binary, 0: text
			{
				/*
				mFILE_IO.Read_data_binary(velocity.u, full_file_path, frame, "u");
				mFILE_IO.Write_data(velocity.u, mWidth, mHeight, mDepth, converted_file_path, frame, "u");

				mFILE_IO.Read_data_binary(velocity.v, full_file_path, frame, "v");
				mFILE_IO.Write_data(velocity.v, mWidth, mHeight, mDepth, converted_file_path, frame, "v");

				mFILE_IO.Read_data_binary(velocity.w, full_file_path, frame, "w");
				mFILE_IO.Write_data(velocity.w, mWidth, mHeight, mDepth, converted_file_path, frame, "w");
				*/

				mFILE_IO.Read_data_binary(density.value, full_file_path, frame, "density");
				mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, converted_file_path, frame, "density");
				
				//mFILE_IO.Read_data_binary(density.value, full_file_path, frame, "temperature");
				//mFILE_IO.Write_data(density.value, mWidth, mHeight, mDepth, converted_file_path, frame, "temperature");
			}
			else
			{
				/*
				mFILE_IO.Read_data(velocity.u, full_file_path, frame, "u");
				mFILE_IO.Write_data_binary(velocity.u, mWidth, mHeight, mDepth, converted_file_path, frame, "u");
				mFILE_IO.Read_data(velocity.v, full_file_path, frame, "v");
				mFILE_IO.Write_data_binary(velocity.v, mWidth, mHeight, mDepth, converted_file_path, frame, "v");
				mFILE_IO.Read_data(velocity.w, full_file_path, frame, "w");
				mFILE_IO.Write_data_binary(velocity.w, mWidth, mHeight, mDepth, converted_file_path, frame, "w");
				mFILE_IO.Read_data(density.value, full_file_path, frame, "density");
				mFILE_IO.Write_data_binary(density.value, mWidth, mHeight, mDepth, converted_file_path, frame, "density");
				*/
			}
			/*
			mFILE_IO.Read_data_binary(pbd.obj->facePositionData, full_file_path, mFsize*12, frame, "f_pos");
			mFILE_IO.Read_data_binary(pbd.obj->faceNormalData, full_file_path, mFsize*12, frame, "f_norm");
			mFILE_IO.Write_data_binary(pbd.obj->facePositionData, mFsize*12, converted_file_path, frame,"f_pos");
			mFILE_IO.Write_data_binary(pbd.obj->faceNormalData, mFsize*12, converted_file_path, frame,"f_norm");
			*/
			printf("# %5d frame data converted.\n", frame);
		}
	}

	option_convert = 0;
}

/*
----------------------------------------------------------------------
GLUT callback routines
----------------------------------------------------------------------
*/
static void Key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'b':
	case 'B':


	case 'r':
	case 'R':
		Init_data();
		//delete_data();
		break;

	case 'c':
	case 'C':
		//clear_data ();
		break;

	case 'q':
	case 'Q':
	case 27:	// esc
		exit(0);
		break;

	case 'l':
	case 'L':
		option_show_levelset = !option_show_levelset;
		option_show_temperature = 0;
		option_show_phi = 0;
		break;

	case 't':
	case 'T':
		option_show_temperature = !option_show_temperature;
		option_show_levelset = 0;
		option_show_phi = 0;
		
		if(option_show_temperature)	mView_mode = "temperature";
		else mView_mode = "density";
		break;

	case 'd':
	case 'D':
		option_show_levelset = 0;
		option_show_temperature = 0;
		option_show_phi = 0;
		break;

	case 'p':
	case 'P':
		option_show_phi = !option_show_phi;
		option_show_levelset = 0;
		option_show_temperature = 0;
		break;


	case 'v':
	case 'V':
		option_show_velocity = !option_show_velocity;
		option_show_phi = 0;
		option_show_levelset = 0;
		option_show_temperature = 0;
		break;

	case 32:	// space bar
		option_pause = !option_pause;
		break;

		
	case '6':
		option_show_particle = !option_show_particle;
		option_show_particle_only = false;
		break;
		
	
	case '7':
		option_show_particle_only = !option_show_particle_only;
		option_show_particle = false;
		cout << "MODE : SHOW PARTICLE ONLY" << endl;
		break;

	case 'a':
	case 'A':
		option_slide = !option_slide;
		break;

	case 's':
	case 'S':
		option_sphere = !option_sphere;
		break;

	case ']':
		option_convert = 1;
		break;

	case '[':
		option_convert = 2;
		break;

	case '<':
	case ',':
		if (mFrame_no > mFrame_start)
		{
			mFrame_no--;
			Read_data(mFrame_no);;
			glutPostRedisplay();
		}
		break;

	case '>':
	case '.':
		if (mFrame_no < mFrame_end)
		{
			mFrame_no++;
			Read_data(mFrame_no);
			glutPostRedisplay();
		}
		break;

	case '+':
	case '=':
		if (option_range >= 1.0f) option_range = 1.0f;
		else option_range += 0.005f;

		printf("option_range : %f\n", option_range);
		break;

	case '-':
	case '_':
		if (option_range <= 0.001f) option_range = 0.001f;
		else option_range -= 0.005f;

		printf("option_range : %f\n", option_range);
		break;
		
	}
}

static void Mouse_func(int button, int state, int x, int y)
{
	mMouse_x0 = mMouse_x = x;
	mMouse_y0 = mMouse_y = y;

	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		mMouse_down[0] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	case GLUT_MIDDLE_BUTTON:
		mMouse_down[1] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	case GLUT_RIGHT_BUTTON:
		mMouse_down[2] = ((GLUT_DOWN == state) ? 1 : 0);
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

static void Motion_func(int x, int y)
{
	int diffx = x - mMouse_x;
	int diffy = y - mMouse_y;

	mMouse_x = x;
	mMouse_y = y;

	if (mMouse_down[0] && mMouse_down[2])
	{
		mZoom -= (float)0.05f * diffy;
	}
	else
	{
		if (option_slide && mMouse_down[0])
		{
			if (mSlide < 0) mSlide = mDepth;
			else if (mSlide > mDepth) mSlide = 0;
			else mSlide -= (float) 0.5f * diffy;

			//printf("%f %d\n", mSlide, mDepth);
		}
		else if (mMouse_down[0])
		{
			mRotateX += (float) 0.5f * diffy;
			mRotateY += (float) 0.5f * diffx;
		}
		else if (mMouse_down[2])
		{
			mTransX += (float) 0.01f * diffx;
			mTransY -= (float) 0.01f * diffy;
		}
	}

	glutPostRedisplay();
}

static void Reshape_func(int mWidth, int mHeight)
{
	glutSetWindow(mWin_id);
	glutReshapeWindow(mWin_width, mWin_height);
}

static void Idle_func(void)
{
	if (option_convert) Convert_file_type(mFrame_no - 1);
	if (mFrame_no <= mFrame_end && !option_pause)
	{
		if (mFrame_no % 3 == 2)
			Read_data(mFrame_no);
		mFrame_no++;
	}

	glutSetWindow(mWin_id);
	glutPostRedisplay();
}

static void Draw_particles(void)
{

	int N = mDepth;
	if (mHeight > N) N = mHeight;
	else if (mWidth > N) N = mWidth;

	float h = 10.0f / N;
	//glMatrixMode(GL_MODELVIEW);
	/*
//	printf("Draw fuel particles : %d\n", fuel_particle.mNo);
	glUseProgram(shader->program);
	glEnable(GL_TEXTURE_2D);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, flake_texture);
	glUniform1i(glGetUniformLocation(shader->program, "Flake_texture"), 1);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(7);
	glEnableVertexAttribArray(8);
	glColor3f(1.0f, 0.7f, 0.0f);	// orange
	//		glPointSize(0.0f);
	glBegin(GL_POINTS);
	printf("%d\n", fuel_particle.mNo);
	for (int i = 0; i < fuel_particle.mNo; i++)
	{
		glVertexAttrib1f(1, fuel_particle.mRadius[i]);
		glVertexAttrib1f(7, fuel_particle.mAngle[i]);
		glVertexAttrib1f(8, fuel_particle.mTemperature[i]);
		glNormal3f(fuel_particle.mNx[i], fuel_particle.mNy[i], fuel_particle.mNz[i]);
		glVertex3f((fuel_particle.mCx[i] - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] - (mDepth + 2) / 2)*h);
	}
	glEnd();
	glUseProgram(0);
	glBindTexture(GL_TEXTURE_2D, 0);
	*/

	//	포인트로 그리기
	/*
	glColor3f(1.0f, 0.7f, 0.0f);	// orange
	glBegin(GL_POINTS);
	//printf("%d\n", fuel_particle.mNo);
	for (int i = 0; i < fuel_particle.mNo; i++)
	{
		if (fuel_particle.mCx[i] <= mWidth && fuel_particle.mCy[i] <= mHeight && fuel_particle.mCz[i] <= mDepth)
			//glVertex3f(fuel_particle.mCx[i], fuel_particle.mCy[i], fuel_particle.mCz[i]);
			glVertex3f((fuel_particle.mCx[i] - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] - (mDepth + 2) / 2)*h);
			//glVertex3f((fuel_particle.mCx[i] - mWidth / 2)*h, (fuel_particle.mCy[i] - mHeight / 2)*h, (fuel_particle.mCz[i] - mDepth / 2)*h);
	}
	glEnd();
	*/

	//	velocity 고려하여 직선으로 그리기
	glColor3f(1.0f, 0.7f, 0.0f);	// orange
	glBegin(GL_LINES);
	float max_tail_length = 3.0f;
	//float min_tail_length = 0.0001f;

	//printf("%d\n", fuel_particle.mNo);
	for (int i = 0; i < fuel_particle.mNo; i++)
	{
		if (fuel_particle.mCx[i] <= mWidth && fuel_particle.mCy[i] <= mHeight && fuel_particle.mCz[i] <= mDepth)
		{
			glVertex3f((fuel_particle.mCx[i] - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] - (mDepth + 2) / 2)*h);

			float tail_x = fuel_particle.mNx[i] * 0.3f;
			float tail_y = fuel_particle.mNy[i] * 0.3f;
			float tail_z = fuel_particle.mNz[i] * 0.3f;

			if (abs(fuel_particle.mNx[i]) > max_tail_length) tail_x = (tail_x / tail_x) * max_tail_length;
			if (abs(fuel_particle.mNy[i]) > max_tail_length) tail_y = (tail_y / tail_y) * max_tail_length;
			if (abs(fuel_particle.mNz[i]) > max_tail_length) tail_z = (tail_z / tail_z) * max_tail_length;
			//if (abs(tail_y) > max_tail_length) tail_y = ((tail_y / tail_y) * max_tail_length - (mHeight + 2) / 2)*h;
			//if (abs(tail_z) > max_tail_length) tail_z = ((tail_z / tail_z) * max_tail_length - (mDepth + 2) / 2)*h;

			glVertex3f((fuel_particle.mCx[i] + tail_x - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] + tail_y - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] + tail_z - (mDepth + 2) / 2)*h);
			//glVertex3f((fuel_particle.mCx[i] + tail_x - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] + tail_y - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] + tail_z - (mDepth + 2) / 2)*h);

			//if (tail_y > max_tail_length) tail_y = max_tail_length;
			//if (tail_z > max_tail_length) tail_z = max_tail_length;

			//if (tail_x < min_tail_length) tail_x = min_tail_length;
			//if (tail_y < min_tail_length) tail_y = min_tail_length;
			//if (tail_z < min_tail_length) tail_z = min_tail_length;


			//glVertex3f(tail_x, tail_y, tail_z);
		}
	}
	glEnd();

	/*
	// 노말 그리는 부분. 
	glColor3f(0.7, 1.0, 1.0);
	glBegin(GL_LINES);
	for (int i = 0; i < fuel_particle.mNo; i++)
	{
	glVertex3f((fuel_particle.mCx[i] - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] - (mDepth + 2) / 2)*h);
	glVertex3f((fuel_particle.mCx[i] + fuel_particle.mNx[i] / 1.0 - (mWidth + 2) / 2)*h, (fuel_particle.mCy[i] + fuel_particle.mNy[i] / 1.0 - (mHeight + 2) / 2)*h, (fuel_particle.mCz[i] + fuel_particle.mNz[i] / 1.0 - (mDepth + 2) / 2)*h);
	}
	glEnd();
	*/
}

static void Display_func(void)
{
	Pre_display();

	Draw_base_line();

	string window_title;
	window_title = "GRIDY | GL_VIEWER v0.6.0, res(" + std::to_string(mWidth) + "x" + to_string(mHeight) + "x" + to_string(mDepth) + "), " + to_string(mFrame_no-1) + "/" + to_string(mFrame_end) + " " + mView_mode;
	glutSetWindowTitle(window_title.c_str());

	//draw_grid_line();
	//	if ( dvel ) draw_velocity ();
	//	else draw_density_3d();


	if (option_sphere)
	{
		int N = mDepth;
		if (mHeight > N) N = mHeight;
		else if (mWidth > N) N = mWidth;

		float h = 10.0f / N;
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef((mSphere_x - mWidth / 2)*h, (mSphere_y - mHeight / 2)*h, (mDepth / 2 - mSphere_z)*h);
		glutSolidSphere((mSphere_r + 0.5f)*h, 100, 100);
		glPopMatrix();

	}




//	if (option_pbd)
//		Draw_pbd_3d();
	//else

	if (option_show_particle)
	{
		Draw_particles();
	}
	if (!option_show_particle_only)
	{
		if (option_show_velocity) Draw_velocity_3D();
		else if (!option_show_particle_only) Draw_density_3D();
	}



	/*
	if(option_phi)
	{
	float x, y, z;
	float s = 0.2f;

	int N = mDepth;
	if(mHeight > N) N=mHeight;
	else if(mWidth > N) N=mWidth;

	float h = 10.0f/N;

	float normal[3];
	float ph0, ph1, ph2;


	glBegin(GL_LINES);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	for(int i=1; i<=mWidth ; i++)
	{
	x = (i-mWidth/2-0.5f)*h;

	for(int j=1 ; j<=mHeight ; j++)
	{
	y = (j-mHeight/2-0.5f)*h;

	for(int k=1; k<=mDepth; k++)
	{
	z = (mDepth/2-k-0.5f)*h;

	if(abs(phi[IX3(i,j,k)]) <= 1.732)
	{
	ph0 = phi[IX3(i+1,j,k)];
	ph1 = phi[IX3(i,j,k)];
	ph2 = phi[IX3(i-1,j,k)];
	normal[0] = ((ph1-ph0) + (ph2-ph1));

	ph0 = phi[IX3(i,j+1,k)];
	ph1 = phi[IX3(i,j,k)];
	ph2 = phi[IX3(i,j-1,k)];
	normal[1] = ((ph1-ph0) + (ph2-ph1));

	ph0 = phi[IX3(i,j,k+1)];
	ph1 = phi[IX3(i,j,k)];
	ph2 = phi[IX3(i,j,k-1)];
	normal[2] = ((ph1-ph0) + (ph2-ph1));

	// normalize
	float det = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
	det = sqrtf(det);

	if(det == 0.0f)	det = 1.0f;
	det = 1.0f/det;

	normal[0] *= det;
	normal[1] *= det;
	normal[2] *= det;
	glVertex3f(x, y, z);
	glVertex3f(x + normal[0]*s, y + normal[1]*s, z + normal[2]*s);

	//printf("%f %f %f -> %f %f %f\t", x, y, z, x + normal[0]*s, y + normal[1]*s, z + normal[2]*s);
	}
	}
	}
	}

	glEnd();
	}
	*/



	Post_display();
}

/*
----------------------------------------------------------------------
open_glut_window --- open a glut compatible window and set callbacks
----------------------------------------------------------------------
*/
static void Open_glut_window(void)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(mWin_width, mWin_height);

	mWin_id = glutCreateWindow("GRIDY | GL_VIEWER v0.6.0");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	Pre_display();

	glutKeyboardFunc(Key_func);
	glutMouseFunc(Mouse_func);
	glutMotionFunc(Motion_func);
	glutReshapeFunc(Reshape_func);
	glutIdleFunc(Idle_func);
	glutDisplayFunc(Display_func);
}


static void delete_data(void)
{
	density.~GRIDY_SCALAR_FIELD();
	velocity.~GRIDY_VECTOR_FIELD();
	temperature.~GRIDY_SCALAR_FIELD();
	fuel_particle.~GRIDY_IMPLICIT_OBJECT_SPHERES();
	mFILE_IO.Close_file_stream();
}

static void Init_data(void)
{
	printf("%s\n", full_file_path.c_str());

	mFILE_IO.Read_simulation_info(full_file_path, mDimension, mWidth, mHeight, mDepth, mFrame_start, mFrame_end, mFile_type);
	//mFILE_IO.Read_pbd_simulation_info(full_file_path, mPsize, mFsize, mFrame_start, mFrame_end);

	/*
	vector3f pbd_st,pbd_ed;
	vector3f smoke_st, smoke_ed;
	mFILE_IO.Read_pbd_couple_simulation_info(full_file_path, mPsize, mFsize, mWidth, mHeight, mDepth, pbd_st,pbd_ed,smoke_st,smoke_ed,mFrame_start, mFrame_end);
	*/
	printf("filepath : %s\n size : %d %d %d, (%d to %d)\n", full_file_path.c_str(), mWidth, mHeight, mDepth, mFrame_start, mFrame_end);

	int size = (mWidth + 2)*(mHeight + 2)*(mDepth + 2);

	density.Init(mWidth, mHeight, mDepth);
	velocity.Init(mWidth, mHeight, mDepth);
	temperature.Init(mWidth, mHeight, mDepth);


	/*
	pbd.obj = new PBDobject();
	pbd.obj->bind(mPsize, mFsize);
	*/
	//dens = new float[(mWidth+2)*(mHeight+2)*(mDepth+2)];

	mFrame_no = mFrame_start;
}

static void glInit(void)
{
	//	define windows size
	mWin_width = 1600;
	mWin_height = 980;

	Set_viewer_attr();
	Open_glut_window();
	//glEnable(GL_POINT_SPRITE);
	//glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glEnable(GL_DEPTH_TEST);

	//Light_init();


	/*
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glEnable(GL_LIGHTING);

	GLfloat ambientLight[] = {0.3f,0.3f,0.3f,1.0f};
	GLfloat diffuseLight[] = {0.7f,0.7f,0.7f,1.0f};
	GLfloat specular[]  = {1.0f,1.0f,1.0f,1.0f};
	GLfloat specref[]  = {1.0f,1.0f,1.0f,1.0f};
	GLfloat position[]  = {10.0f,10.0f,10.0f,1.0f};

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glEnable(GL_LIGHT0);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
	glMateriali(GL_FRONT, GL_SHININESS, 128);
	*/
}

/*
----------------------------------------------------------------------
main --- main routine
----------------------------------------------------------------------
*/
int main(int argc, char ** argv)
{

	glutInit(&argc, argv);

	glInit();
	//Shader_init();


	// TO DO : parsing하는 모듈이 들어가야 함
	if (argc == 2)
	{
		full_file_path = argv[1];
		//full_file_path = (char*)malloc(sizeof(char)*strlen(argv[1]));
		//full_file_path = argv[1];
		//strcpy_s(full_file_path, argv[1]);
	}
	else
	{
		full_file_path = GRIDY_RESULT_PATH;
//		full_file_path = "F:\\GRIDY_RESULT1";
		//full_file_path = "E:\\GRIDY\\RESULT";
//		full_file_path = "F:\\GRIDY_RESULT_ty2";

//		full_file_path = "F:\\GRIDY_RESULT_firehell_1";
//		full_file_path = "F:\\fire_tongtong_1";
//		full_file_path = "F:\\fire_obstacle";
//		full_file_path = "C:\\gridy_firehell";
//		full_file_path = "F:\\cFire_top_disspate(19618777)";
		/*
		cout << "USAGE : \nGRIDY_GL_VIEWER.exe <SIMULATION_DATA_DIRECTORY>" << endl;
		exit(0);
		*/
	}
	Init_data();
	// 임시 코드 끝

	option_show_velocity = 0;
	option_show_levelset = 0;
	option_show_temperature = 0;
	option_show_phi = 0;

	mSlide = mDepth;



	glutMainLoop();

	delete_data();
	exit(0);
}




 