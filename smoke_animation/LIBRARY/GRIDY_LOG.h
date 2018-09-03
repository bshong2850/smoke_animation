/*
 *	GRIDY_LOG
 *
 *	�ùķ��̼� �� �ܰ躰 �ð��� �����ϰ�, ȭ�鿡 ����ϸ�, ���Ϸ� �����Ѵ�.
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

	int mTab_size;							// ����-�� ���� �� ���� �� ������. �⺻ ���� 4
	int mSpace_level;						// �۾��� ������ �� ������ ������ �� �� �� ���� �ּ� �ܰ�. �⺻ ���� 2.
	string mPrev_work;						// ���� �۾� �̸�
	
	bool mIs_file_output = false;			// ���Ϸ� ������� ����
	string mLog_file_path, mLog_file_name;	// ���� ���(�� �������� \ ���� ����)�� ���� �̸�
	FILE *mFP;								// ���� ������




	// tab_size��ŭ �� ĭ�� ȭ�鿡 ���
	inline void Print_tab(int tab_size)
	{
		for(int i=0; i<tab_size; i++) 
			cout << " ";
	}

	// tab_size��ŭ �� ĭ�� ���Ͽ� ���
	inline void Print_tab_to_file(int tab_size)
	{
		for (int i = 0; i<tab_size; i++)
			fprintf(mFP, " ");
	}

	// ȭ�鿡 ����ϴ� ������ ������ ����
	inline void Set_text_color(int text_color)
	{
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), text_color);
	}

	// stack�� ����(�ܰ�) ���� �ٸ� ������ ����
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
		if (mFP) fclose(mFP);	// ������ ���������� ������ ����
	}



	void GRIDY_LOG::Init(int tab_size, int space_level)
	{
		mTab_size = tab_size;
		mSpace_level = space_level;
		mPrev_work = "";

		mIs_file_output = false;
	}

	
	// �ʱ�ȭ
	void GRIDY_LOG::Init(string log_file_path, string log_file_name, int tab_size, int space_level)
	{
		mTab_size = tab_size;
		mSpace_level = space_level;
		mPrev_work = "";

		Set_file_out(log_file_path, log_file_name, tab_size, space_level);	// ���� ��� ����� �Ʒ��� setting�� ȣ���Ͽ� �ʱ�ȭ
	}
	


	// ���� ����� ���� ����
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





	// work stack�� ���ο� �۾��� ����
	void In(string work)
	{
		clock_t in_time = clock();	// ���� �ð�

		Set_text_color_for_stack((int)mWork_stack.size());

		Print_tab(mTab_size * (int)mWork_stack.size());
		cout << "+ " << work << endl;

		if (mIs_file_output)	// ���Ϸ� ���
		{
			Print_tab_to_file((int)mTab_size * (int)mWork_stack.size());
			fprintf(mFP, "+ %s\n", work.c_str());
		}
		
		Set_text_color(15);

		mWork_stack.push(work);		// �۾� ���ÿ� �۾� ���� ����
		mTime_stack.push(in_time);	// �ð� ���ÿ� ���� �ð��� ����

		mPrev_work = work;			// ���� �۾��� ����
	}




	// work stack�� �־�� �۾��� ������ �� ȣ��
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


		work = mWork_stack.top();	// �۾� ������ �� ��(���� �ֱ�) �۾� ���� ���� ��



		mWork_stack.pop();
		Set_text_color_for_stack((int)mWork_stack.size());
		
		print_tab_size = mTab_size * (int)mWork_stack.size();
		Print_tab(mTab_size * (int)mWork_stack.size());
		if (mIs_file_output) Print_tab_to_file((mTab_size * (int)mWork_stack.size()));
		
		in_time = mTime_stack.top();
		mTime_stack.pop();
		out_time = clock();

		if(mPrev_work == work)	// �߰��� �ٸ� �۾��� ����. ��, ���� �۾��� ����.
		{
			cout << "- ";
			cout.width(75 - (print_tab_size + 2));	// +- ���޾Ƽ� �۾� ���� ������ �ʵ���

			if(mIs_file_output)
			{
				fprintf(mFP, "- ");
				Print_tab_to_file((75 - (print_tab_size + 2)));	// +- ���޾Ƽ� �۾� ���� ������ �ʵ���
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
			cout << endl;		// �̹� ������ pop �Ǹ鼭 level�� �ϳ� �پ���
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