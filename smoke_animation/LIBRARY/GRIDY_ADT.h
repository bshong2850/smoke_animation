#pragma once
#include <list>
#include <stack>
#include <iostream>
using namespace std;

template <typename T>
class GRIDY_LIST_STACK3
{
	list<T> m_stack;

public:
	GRIDY_LIST_STACK3(void) {}
	~GRIDY_LIST_STACK3(void) {}

	void push(T item1, T item2, T item3) 
	{
		m_stack.push_front(item1); 
		m_stack.push_front(item2); 
		m_stack.push_front(item3); 
	}

	T pop()
	{
		T r;

		r = m_stack.front();
		m_stack.pop_front();

		return r;
	}
	//T getTop() { return m_stack.front(); }
	int getSize() { return m_stack.size(); }
};




/*
template <typename T>
class GRIDY_LIST_STACK2
{
	list<T> m_stack;

public:
	GRIDY_LIST_STACK2(void) {}
	~GRIDY_LIST_STACK2(void) {}

	void push(T item1, T item2)
	{
		m_stack.push_front(item1);
		m_stack.push_front(item2);
	}

	T pop()
	{
		T result;

		result = m_stack.front();
		m_stack.pop_front();

		return result;
	}

	int getSize() { return m_stack.size(); }
};
*/


/*
template <typename T>
class GRIDY_STACK
{
	stack<int> m_stack;

public:
	GRIDY_STACK(void) {}
	~GRIDY_STACK(void) {}

	void push(T item1, T item2, T item3)
	{
		m_stack.push(item1);
		m_stack.push(item2);
		m_stack.push(item3);
	}
	
	T* top(void)
	{
		return {m_stack.
};
*/