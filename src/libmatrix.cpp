#ifndef LIBMATRIX_CPP
#define LIBMATRIX_CPP

#include "../include/libmatrix.h"
#include <iostream>
#include <exception>

using namespace libmatrix;
using namespace std;

/**
 * Addresses the i-th element of the vector
 * @param i the index of the element
 */
template<int n, class T>
T Vector<n,T>::at(int i)
{
	return vec[i];
}

/**
 * Cross product with another vector
 * @param v the other vector
 */
template<int n, class T>
float Vector<n,T>::cross(Vector<n,T> &v)
{
	return 0;
}

/**
 * Dot product with another vector
 * @param v the other vector
 */
template<int n, class T>
float Vector<n,T>::dot(Vector<n,T> &v)
{
	return 0;
}

/**
 * Returns true if the vector is orthogonal to another given as an argument, false otherwise
 * @param v the other vector
 */
template<int n, class T>
bool Vector<n,T>::is_ortho(Vector<n,T> &v)
{
	return true;
}

/**
 * Returns true if the vector contains an invalid value, false otherwise. Notably, if the vector contains nan as values
 */
template<int n, class T>
bool Vector<n,T>::is_null()
{
	return true;
}

/**
 * Returns the norm of the vector
 */
template<int n, class T>
float Vector<n,T>::norm()
{
	return 0;
}

/**
 * Returns a copy of the vector normalised.
 */
template<int n, class T>
Vector<n,T>& Vector<n,T>::to_unit()
{

}


template<int n, class T>
ostream& Vector<n,T>::operator<<(ostream &out)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator[](int i)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator+(Vector<n,T>& v)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator+=(Vector<n,T>& v)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator-(Vector<n,T>& v)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator-=(Vector<n,T>& v)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator*(Vector<n,T>& v)
{

}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator*=(Vector<n,T>& v)
{

}

int main(int argc, char const *argv[])
{
	return 0;
}

#endif