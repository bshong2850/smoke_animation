//	
//	GRIDY_PARSER.h
//
//
//
//	This file is a part of GRIDY library.
//	Designed by TaeHyeong Kim (usemagic@gmail.com)
//



#pragma once
#ifndef _GRIDY_PARSER_H_
#define _GRIDY_PARSER_H_



#include "GRIDY_FILE_IO.h"
//#include "GRIDY_ADT.h"
#include <map>
#include <string>

#define buf_line 1024	// 읽어들일 한 줄의 최대 크기






//template <typename T>
class GRIDY_PARSER : GRIDY_FILE_IO
{
private:
	//FILE* fp;
	//GRIDY_FILE_IO mFILE_IO;
//	list<string> m_stack_id;
//	list<string> m_stack_value;
	map<string, string> m_map;
	map<string, string>::const_iterator find_iter;
	void parse_init_file(FILE *fp);



public:
	GRIDY_PARSER(void);
	~GRIDY_PARSER(void);

	void parse_setup(string full_file_path_name);
	void put(string id, string value);
	//void put(string id, int value);
	void put(string id, double value);		// 주의 : 실험 결과 소수점 5자리 이하는 무시되는 경향이 있음

	void get(string id, int* vari, int default);
	void get(string id, float* vari, float default);
	void get(string id, double* vari, double default);
	void get(string id, string* vari, string default);

	void show_all_items(void);



private:
	template<typename T, typename P>
	T remove_if(T beg, T end, P pred)		// by stackoverflow (http://stackoverflow.com/questions/83439/remove-spaces-from-stdstring-in-c)
	{
		T dest = beg;
		for (T itr = beg; itr != end; ++itr)
			if (!pred(*itr))
				*(dest++) = *itr;
		return dest;
	}	
};




GRIDY_PARSER::GRIDY_PARSER(void)
{

}




GRIDY_PARSER::~GRIDY_PARSER(void)
{
	m_map.clear();
	//Close_file_stream();
}





void GRIDY_PARSER::parse_setup(string full_path_file_name)
{
	char buffer[buf_line];

	Open_file_stream(full_path_file_name, "r");

	while (fgets(buffer, buf_line, fp) != NULL)
	{
		if (buffer[0] != '#')	// comment
		{
			string temp_str(buffer);
			string::size_type temp_pos = temp_str.find_first_of('=');	// '='의 위치를 찾는다

			if (temp_pos != string::npos)	// 해당 line이 id = data 구조라면
			{
				
				string temp_left = temp_str.substr(0, temp_pos);
				string temp_right = temp_str.substr(temp_pos + 1);

				temp_left.erase(remove_if(temp_left.begin(), temp_left.end(), isspace), temp_left.end());
				temp_right.erase(remove_if(temp_right.begin(), temp_right.end(), isspace), temp_right.end());

				put(temp_left, temp_right);
			}
		}
		
	}


	//show_all_items();



	Close_file_stream();
}





void GRIDY_PARSER::put(string id, string value)
{
	m_map.insert( map<string, string>::value_type(id, value) );
}


/*
void GRIDY_PARSER::put(string id, int value)
{
	m_map.insert(map<string, string>::value_type(id, to_string(value)));
}
*/


void GRIDY_PARSER::put(string id, double value)
{
	//cout << "insert : " << value << endl;
	m_map.insert( map<string, string>::value_type(id, to_string(value)));
}





void GRIDY_PARSER::get(string id, int* vari, int default = 0)
{
	find_iter = m_map.find(id);

	if (find_iter == m_map.end())
	{
//		cout << id << " is not found." << endl;
		*vari = default;
	}
	else
		*vari = atoi((find_iter->second).c_str());
}





void GRIDY_PARSER::get(string id, float* vari, float default = 0.0f)
{
	find_iter = m_map.find(id);

	if (find_iter == m_map.end())
	{
//		cout << id << " is not found." << endl;
		*vari = default;
	}
	else
		*vari = (float)atof( (find_iter->second).c_str() );
}





void GRIDY_PARSER::get(string id, double* vari, double default = 0.0)
{
	find_iter = m_map.find(id);

	if (find_iter == m_map.end())
	{
		//		cout << id << " is not found." << endl;
		*vari = default;
	}
	else
		*vari = atof((find_iter->second).c_str());
}




//string GRIDY_PARSER::get(string id)
void GRIDY_PARSER::get(string id, string* vari, string default = "")
{
	find_iter = m_map.find(id);

	if (find_iter == m_map.end())
	{
		//		cout << id << " is not found." << endl;
		*vari = default;
	}
	else
		*vari = (find_iter->second);
}






void GRIDY_PARSER::show_all_items(void)
{
	for(find_iter = m_map.begin(); find_iter != m_map.end(); ++find_iter)
	{
		cout << find_iter->first << " : " << find_iter->second << endl;
	}
	cout << endl;

}







#endif