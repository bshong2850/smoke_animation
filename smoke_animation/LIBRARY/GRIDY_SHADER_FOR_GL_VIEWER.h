#pragma warning(disable:4996)




//#include <GL/gl.h>
//#include <GL/glu.h>


#include <iostream>
#include <stdio.h>
#include <GL/glew.h>
//#include <GL/glut.h>





using namespace std;

class GRIDY_SHADER_FOR_GL_VIEWER
{
public:
	GLuint program;
	GLuint vaoHandle;

	GLuint positionBufferHandle;
	GLuint normalBufferHandle;
	GLuint colorBufferHandle;
	GLuint depthBufferHandle;

public:
	GRIDY_SHADER_FOR_GL_VIEWER(void);
	~GRIDY_SHADER_FOR_GL_VIEWER(void);

public:
	char* readShaderSource(const char* shaderFile);
	void InitShader(const char* vShaderFile, const char* fShaderFile);
	void Link();
};


GRIDY_SHADER_FOR_GL_VIEWER::GRIDY_SHADER_FOR_GL_VIEWER(void)
{

}


GRIDY_SHADER_FOR_GL_VIEWER::~GRIDY_SHADER_FOR_GL_VIEWER(void)
{

}


char* GRIDY_SHADER_FOR_GL_VIEWER::readShaderSource(const char* shaderFile)
{
    FILE* fp = fopen(shaderFile, "r");

    if ( fp == NULL ) { return NULL; }

    fseek(fp, 0L, SEEK_END);
    long size = ftell(fp);

    fseek(fp, 0L, SEEK_SET);
    char* buf = new char[size + 1];
	for(int i=0; i<size+1; i++) buf[i] = 0;

    fread(buf, 1, size, fp);

    buf[size] = '\0';
    fclose(fp);

    return buf;
}


void GRIDY_SHADER_FOR_GL_VIEWER::InitShader(const char* vShaderFile, const char* fShaderFile)
{
    struct Shader {
	const char*  filename;
	GLenum       type;
	GLchar*      source;
    }  shaders[2] = {
	{ vShaderFile, GL_VERTEX_SHADER, NULL },
	{ fShaderFile, GL_FRAGMENT_SHADER, NULL }
    };
	
    program = glCreateProgram();
    
	
	for (int i = 0; i < 2; ++i) {
		Shader& s = shaders[i];
		s.source = readShaderSource(s.filename);
		if (shaders[i].source == NULL) {
			std::cerr << "Failed to read " << s.filename << std::endl;
			exit(EXIT_FAILURE);
		}

		GLuint shader = glCreateShader(s.type);

		glShaderSource(shader, 1, (const GLchar**)&s.source, NULL);
		glCompileShader(shader);

		GLint  compiled;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
		if (!compiled) {
			std::cerr << s.filename << " failed to compile:" << std::endl;
			GLint  logSize;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logSize);
			char* logMsg = new char[logSize];
			glGetShaderInfoLog(shader, logSize, NULL, logMsg);
			std::cerr << logMsg << std::endl;
			delete[] logMsg;

			exit(EXIT_FAILURE);
		}

		delete[] s.source;

		glAttachShader(program, shader);
	}


    /* use program object */
    //glUseProgram(program);
}
void GRIDY_SHADER_FOR_GL_VIEWER::Link()
{

	/* link  and error check */
	glLinkProgram(program);

	GLint  linked;
	glGetProgramiv(program, GL_LINK_STATUS, &linked);
	if (!linked) {
		std::cerr << "Shader program failed to link" << std::endl;
		GLint  logSize;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logSize);
		char* logMsg = new char[logSize];
		glGetProgramInfoLog(program, logSize, NULL, logMsg);
		std::cerr << logMsg << std::endl;
		delete[] logMsg;

		exit(EXIT_FAILURE);
	}
}
/*
GRIDY_SHADER_FOR_GL_VIEWER *shader;


static void shader_init()
{
	GLenum err = glewInit();

	if(GLEW_OK != err)
	{
		cout << "hull\n" << endl;
		cout << glewGetErrorString(err) <<endl;
	}
	shader = new GRIDY_SHADER_FOR_GL_VIEWER();
	shader->InitShader("GRIDY_shader_v.glsl", "GRIDY_shader_f.glsl");

	// create and set-up the vertex array object
	glGenVertexArrays(1, &shader->vaoHandle);
	glBindVertexArray(shader->vaoHandle);
}
*/
/***********************SHADER DONE*********************************/