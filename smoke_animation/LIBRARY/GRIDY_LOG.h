/*
 *	GRIDY_LOG
 *
 *	시뮬레이션 각 단계별 시간을 측정하고, 화면에 출력하며, 파일로 저장한다.
 *
*/

#ifndef _GRIDY_LOG_H_
#define _GRIDY_LOG_H_

#pragma warning(disable:4996)




#include <stack>
#include <iostream>
#include <string>
#include <Windows.h>
#include <time.h>




using namespace std;

class GRIDY_LOG
{
private :
	stack<string>	mWork_stack;
	stack<clock_t>	mTime_stack;

	clock_t	mTime;

	int mTab_size;							// 스텝-인 했을 때 띄우는 탭 사이즈. 기본 값은 4
	int mSpace_level;						// 작업이 끝났을 때 구별이 쉽도록 한 줄 더 띄우는 최소 단계. 기본 값은 2.
	string mPrev_work;						// 이전 작업 이름
	
	bool mIs_file_output = false;			// 파일로 출력할지 여부
	string mLog_file_path, mLog_file_name;	// 파일 경로(맨 마지막에 \ 붙지 않음)와 파일 이름
	FILE *mFP;								// 파일 포인터




	// tab_size만큼 빈 칸을 화면에 출력
	inline void Print_tab(int tab_size)
	{
		for(int i=0; i<tab_size; i++) 
			cout << " ";
	}

	// tab_size만큼 빈 칸을 파일에 출력
	inline void Print_tab_to_file(int tab_size)
	{
		for (int i = 0; i<tab_size; i++)
			fprintf(mFP, " ");
	}

	// 화면에 출력하는 문자의 색상을 결정
	inline void Set_text_color(int text_color)
	{
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), text_color);
	}

	// stack의 레벨(단계) 별로 다른 색상을 정의
	void Set_text_color_for_stack(int stack_size)
	{
		if(stack_size < 1) Set_text_color(11);
		else if(stack_size <= 1) Set_text_color(15);
		else if(stack_size <= 2) Set_text_color(7);
		else Set_text_color(8);
	}



public :
	GRIDY_LOG()
	{
		mTab_size = 4;
		mSpace_level = 2;
		mPrev_work = "";

		mIs_file_output = false;
	}



	~GRIDY_LOG() 
	{
		if (mFP) fclose(mFP);	// 파일이 열려있으면 포인터 해제
	}



	void GRIDY_LOG::Init(int tab_size, int space_level)
	{
		mTab_size = tab_size;
		mSpace_level = space_level;
		mPrev_work = "";

		mIs_file_output = false;
	}

	
	// 초기화
	void GRIDY_LOG::Init(string log_file_path, string log_file_name, int tab_size, int space_level)
	{
		mTab_size = tab_size;
		mSpace_level = space_level;
		mPrev_work = "";

		Set_file_out(log_file_path, log_file_name, tab_size, space_level);	// 파일 출력 기능은 아래의 setting을 호출하여 초기화
	}
	


	// 파일 출력을 위한 셋팅
	void GRIDY_LOG::Set_file_out(string log_file_path, string log_file_name, int tab_size, int space_level)
	{
		mTab_size = tab_size;
		mSpace_level = space_level;
		mPrev_work = "";

		mIs_file_output = true;
		mLog_file_path = log_file_path;
		mLog_file_name = log_file_name;
		
		// current date/time based on current system
		time_t now = time(0);
		tm *ltm = localtime(&now);
		
		// add time info to file name
		mLog_file_name += "_" + to_string(1900 + ltm->tm_year) + "-";	// year
		mLog_file_name += to_string(1 + ltm->tm_mon) + "-";				// month
		mLog_file_name += to_string(ltm->tm_mday) + "_";				// day
		mLog_file_name += to_string(ltm->tm_hour) + "-";
		mLog_file_name += to_string(ltm->tm_min);
		mLog_file_name += ".GRIDY_LOG";

		// init file pointer and open file to write logs
		errno_t err;
		if ((err = fopen_s(&mFP, (mLog_file_path + "\\" + mLog_file_name).c_str(), "w")) != 0)
		{
			printf("Cannot open %s\n", (mLog_file_path + "\\" + mLog_file_name).c_str());
			exit(0);
		}
		else
		{
			Message("Write Simulation Log to\n" + (mLog_file_path + "\\" + mLog_file_name) + "\n");
		}
	}





	// work stack에 새로운 작업을 넣음
	void In(string work)
	{
		clock_t in_time = clock();	// 들어온 시각

		Set_text_color_for_stack((int)mWork_stack.size());

		Print_tab(mTab_size * (int)mWork_stack.size());
		cout << "+ " << work << endl;

		if (mIs_file_output)	// 파일로 출력
		{
			Print_tab_to_file((int)mTab_size * (int)mWork_stack.size());
			fprintf(mFP, "+ %s\n", work.c_str());
		}
		
		Set_text_color(15);

		mWork_stack.push(work);		// 작업 스택에 작업 명을 넣음
		mTime_stack.push(in_time);	// 시간 스택에 현재 시각을 넣음

		mPrev_work = work;			// 이전 작업명 갱신
	}




	// work stack에 넣어둔 작업이 끝났을 때 호출
	void Out(void)
	{
		clock_t in_time, out_time;
		int print_tab_size;
		string work;


		if(mWork_stack.size() < 1 || mTime_stack.size() < 1) 
		{
			cout << "\n\nERROR : GRIDY_LOG stack is empty.\n\n";
			return;
		}


		work = mWork_stack.top();	// 작업 스택의 맨 위(가장 최근) 작업 명을 가져 옴



		mWork_stack.pop();
		Set_text_color_for_stack((int)mWork_stack.size());
		
		print_tab_size = mTab_size * (int)mWork_stack.size();
		Print_tab(mTab_size * (int)mWork_stack.size());
		if (mIs_file_output) Print_tab_to_file((mTab_size * (int)mWork_stack.size()));
		
		in_time = mTime_stack.top();
		mTime_stack.pop();
		out_time = clock();

		if(mPrev_work == work)	// 중간에 다른 작업이 없음. 즉, 서브 작업이 없음.
		{
			cout << "- ";
			cout.width(75 - (print_tab_size + 2));	// +- 연달아서 작업 명이 나오지 않도록

			if(mIs_file_output)
			{
				fprintf(mFP, "- ");
				Print_tab_to_file((75 - (print_tab_size + 2)));	// +- 연달아서 작업 명이 나오지 않도록
			}
		}
		else
		{
			cout << "- " << work;
			cout.width(75 - (print_tab_size + work.length() + 2));

			if (mIs_file_output)
			{
				fprintf(mFP, "- %s", work.c_str());
				Print_tab_to_file((75 - (print_tab_size + work.length() + 2)));
			}
		}

		cout << ((double)(out_time-in_time)/CLOCKS_PER_SEC) << "s." << endl;
		Set_text_color(15);

		if (mIs_file_output) fprintf(mFP, "%.3fs.\n", (double)(out_time - in_time) / CLOCKS_PER_SEC);

		if ((int)mWork_stack.size() < mSpace_level)
		{
			cout << endl;		// 이미 위에서 pop 되면서 level이 하나 줄었음
			if (mIs_file_output) fprintf(mFP, "\n");
		}
	}





	void Message(string message)
	{
		Set_text_color_for_stack(mWork_stack.size());
		
		int print_tab_size = mTab_size * mWork_stack.size();
		Print_tab(mTab_size * mWork_stack.size());

		cout << message << endl;
	}



	void Message(int message)
	{
		Set_text_color_for_stack(mWork_stack.size());
		
		int print_tab_size = mTab_size * mWork_stack.size();
		Print_tab(mTab_size * mWork_stack.size());

		cout << message << endl;
	}



	void Message(float message)
	{
		Set_text_color_for_stack(mWork_stack.size());
		
		int print_tab_size = mTab_size * mWork_stack.size();
		Print_tab(mTab_size * mWork_stack.size());

		cout << message << endl;
	}



	void Message(string message, int stack_size)
	{
		Set_text_color_for_stack(stack_size);
		
		int print_tab_size = mTab_size * stack_size;
		Print_tab(mTab_size * stack_size);

		cout << message << endl;
	}	


};





#endif	//_GRIDY_LOG_H_