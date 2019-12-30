#include "../include/libmatrix.hpp"
#include "../include/libgeometry.hpp"
#include <iostream>
#include <assert.h>

using namespace libmatrix;
using namespace libgeometry;
using namespace std;

int main(int argc, char const *argv[])
{
	Quaternion<float> q;
	cout << q << endl;

	Direction<3,float> d;
	d[1] = 1;
	q = Quaternion<float>(30,d);
	cout << q << endl;
	// cout << d.fullDirection() << endl;
	Point<3,float> a,b;
	for (int i = 0; i < 3; i++)
	{
		a[i] = i * 2;
		b[i] = 1+i/2;
	}
	// cout << a << endl << b << endl;
	Transform t(q);
	cout << t << endl;
	Point<3,float> p({1,1,1});
	cout << t.apply(p) << endl;
	cout << t.to_quat() << endl;
	
	return 0;
}