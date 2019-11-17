#include "../include/libmatrix.hpp"
#include <iostream>

using namespace libmatrix;
using namespace std;

int main(int argc, char const *argv[])
{
	Vector<3,int> v;
	v[0] = 1;
	v[1] = 2;
	v[2] = 3;

	Vector<3,int> v2;
	v2[0] = 5;
	v2[1] = -4;
	v2[2] = 1;
	cout << v << endl;
	cout << v.dot(v2) << endl;
	Vector<2,float> v3;
	v3[0] = 3;
	v3[1] = 2;
	cout << v3.norm() << endl;
	Vector<2,float> unit = v3.to_unit();
	cout << unit << endl;
	v3 *= 5;
	cout << v3 << endl;
	cout << -v3 << endl;
	cout << v3.is_ortho(unit) << endl;
	Vector<2,float> a; a[0] = 1;
	Vector<2,float> b; b[1] = 1;
	cout << a.is_ortho(b) << endl;
	
	return 0;
}