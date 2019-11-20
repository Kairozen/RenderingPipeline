#include "../include/libmatrix.hpp"
#include <iostream>

using namespace libmatrix;
using namespace std;

int main(int argc, char const *argv[])
{
	Matrix<2,3,int> m;
	m[0][1] = 1;
	m[1][1] = 2;
	cout << m << endl;
	cout << m.transpose() << endl;
	Matrix<2,3,int> m2;
	m2[1][1] = 1;
	cout << m + m2 << endl;
	return 0;
}