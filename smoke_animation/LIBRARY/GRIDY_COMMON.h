#ifndef _GRIDY_COMMON_H_
#define _GRIDY_COMMON_H_

#include <Windows.h>
#include <iostream>


using namespace std;


#define IX3(i,j,k)		((i) + (mWidth+2)*(j) + (mWidth+2)*(mHeight+2)*(k))
//#define IX3(i,j,k)	((i) + (width+2)*(j) + (width+2)*(height+2)*(k))
//#define SWAP(x0,x)	{float *tmp=x0; x0=x; x=tmp;}
#define MIN(i, j)		((i>j)? j : i)
#define MAX(i, j)		((i>j)? i : j)
//#define NORM2(i, j, k)	sqrt(i*i + j*j + k*k)


class GRIDY_COMMON
{
private:
//	errno_t err;

public:
	void Terminate_system(std::string msg)
	{
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
		
		cout << endl << endl;
		for(int i=0; i<80; i++) cout << "=";
		

		cout << "\t\t\t\t[ERROR REPORT]" << endl << msg << endl;
		
		for(int i=0; i<80; i++) cout << "="; cout << endl << endl;
		
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
		exit(0);
	}



};


#endif