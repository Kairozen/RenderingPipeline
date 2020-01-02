#include "../include/libmatrix.hpp"
#include <iostream>
#include <assert.h>

using namespace libmatrix;
using namespace std;

int main(int argc, char const *argv[])
{
	float epsilon = 0.000001;

	Matrix<4,3, int> m1;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m1[i][j] = i+j;
		}
	}
	
	assert(m1[0][1] == 1);
	assert(m1[1][1] == 2);
	assert(m1[2][1] == 3);
	assert(m1.at(0,0) == 0);
	assert(m1.at(3,2) == 5);
	
	bool exceptionThrown = false;
	try
	{
		m1.at(2,5);
	}
	catch (exception e)
	{
		exceptionThrown = true;
	}
	assert(exceptionThrown);
	
	exceptionThrown = false;
	try
	{
		m1.at(5,2);
	}
	catch (exception e)
	{
		exceptionThrown = true;
	}
	assert(exceptionThrown);

	Matrix<4,3,int> m2;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m2[i][j] = 1;
		}
	}

	m1 += m2;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			assert(m1.at(i,j) == i+j+1);
		}
	}

	Matrix<4,3,int> m3 = m1 + m2;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			assert(m3.at(i,j) == i+j+2);
		}
	}

	Matrix<3,3,float> test_inv;
	test_inv[0][0] = 2, test_inv[0][1] = 3, test_inv[0][2] = 8;
	test_inv[1][0] = 6, test_inv[1][1] = 0, test_inv[1][2] = -3;
	test_inv[2][0] = -1, test_inv[2][1] = 3, test_inv[2][2] = 2;
	//cout << test_inv << endl;
	Matrix<3,3,float> inv = test_inv.inverse();
	// cout << inv*test_inv << endl;
	// cout << inv << endl;
	// cout << Identity44r << endl;
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			if(i==j)
				assert(Identity44r.at(i,j) == 1);
			else
				assert(Identity44r.at(i,j) == 0);
		}		
	}

	cout << "All is right" << endl;
	return 0;
}