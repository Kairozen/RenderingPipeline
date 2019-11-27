#include "../include/libmatrix.hpp"
#include <iostream>
#include <assert.h>

using namespace libmatrix;
using namespace std;

int main(int argc, char const *argv[])
{
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
	cout << test_inv << endl;
	Matrix<3,3,float> inv = test_inv.inverse();
	cout << inv*test_inv << endl;
	cout << inv << endl;
	/*
	Vector<3,int> v2;
	v2[0] = 5;
	v2[1] = -4;
	v2[2] = 1;

	float dot = v.dot(v2);
	assert(dot == 1 * 5 + 2 * -4 + 3 * 0);

	Vector<2,int> o1;
	o1[0] = 1;
	Vector<2,int> o2;
	o2[1] = 1;
	assert(o1.at(0) == 1);
	assert(o1.at(1) == 0);
	assert(o2[0] == 0);
	assert(o2[1] == 1);

	assert(o1.is_ortho(o2) == true);
	assert(o2.is_ortho(o1) == true);

	assert(v.is_ortho(v2) == false);
	assert(v2.is_ortho(v) == false);

	Vector<3,float> vn;
	vn[0] = 5.5;
	vn[1] = 5.5;
	vn[2] = 5.5;
	assert(vn.is_null() == false);
	vn[1] = sqrt(-1);
	assert(vn.is_null() == true);

	Vector<2,float> norm;
	norm[0] = 3;
	norm[1] = 2;
	float epsilon = 0.0000001;
	assert(abs(norm.norm() - sqrt(3 * 3 + 2 * 2)) < epsilon);

	Vector<2,float> unit = norm.to_unit();
	assert(unit[0] == (1/norm.norm())*norm[0]);
	assert(unit[1] == (1/norm.norm())*norm[1]);
	assert(unit.is_unit());

	Vector<2,int> op1, op2;
	op1[0] = 2;
	op1[1] = 3;
	op2[0] = 4;
	op2[1] = 5;
	assert((op1+op2)[0] == 6);
	assert((op1+op2)[1] == 8);
	assert((op1-op2)[0] == -2);
	assert((op1-op2)[1] == -2);
	assert((-op1)[0] == -2);
	assert((-op1)[1] == -3);
	op1+=op2;
	assert(op1[0] == 6);
	assert(op1[1] == 8);
	*/
	cout << "All is right" << endl;
	return 0;
}