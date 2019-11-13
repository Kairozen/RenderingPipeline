#ifndef LIBMATRIX_CPP
#define LIBMATRIX_CPP

#include "../include/libmatrix.hpp"
#include <iostream>
#include <exception>
#include <cmath>

using namespace libmatrix;
using namespace std;

/**
 * Addresses the i-th element of the vector
 * @param i the index of the element
 */
template<int n, class T>
T Vector<n,T>::at(int i)
{
	if(i>=n)
		throw exception("Given index is out of range (" + to_string(i) + ").");
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
	float res = 0.f;
	for (int i = 0; i < n; ++i)
		res += vec[i]*v.vec[i];
	return res;
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
	return false;
}

/**
 * Returns the norm of the vector
 */
template<int n, class T>
float Vector<n,T>::norm()
{
	float norm = 0.f;
	for (int i = 0; i < n; ++i)
		norm += vec[i] * vec[i];
	return sqrt(norm);
}

/**
 * Returns a copy of the vector normalised.
 */
template<int n, class T>
Vector<n,T>& Vector<n,T>::to_unit()
{

}

template<int n, class T>
T& Vector<n,T>::operator[](int i)
{
	return vec[i];
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator+(const Vector<n,T>& v)
{
	Vector<n,T> res;
	for (int i = 0; i < n; ++i)
		res[i] = vec[i] + v.vec[i];
	return res;
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator+=(const Vector<n,T>& v)
{
	for (int i = 0; i < n; ++i)
		vec[i] += v.vec[i];
	return *this;
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator-(const Vector<n,T>& v)
{
	Vector<n,T> res;
	for (int i = 0; i < n; ++i)
		res[i] = vec[i] - v.vec[i];
	return res;
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator-=(const Vector<n,T>& v)
{
	for (int i = 0; i < n; ++i)
		vec[i] -= v.vec[i];
	return *this;
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator*(const Vector<n,T>& v)
{
	Vector<n,T> res;
	for (int i = 0; i < n; ++i)
		res[i] = vec[i] * v.vec[i];
	return res;
}

template<int n, class T>
Vector<n,T>& Vector<n,T>::operator*=(const Vector<n,T>& v)
{
	for (int i = 0; i < n; ++i)
		vec[i] *= v.vec[i];
	return *this;
}

#endif