#ifndef LIBMATRIX_HPP
#define LIBMATRIX_HPP

#include <iostream>

using namespace std;

namespace libmatrix
{
	template <int n, class T>
	class Vector
	{
	private:
		T vec[n];

	public:
		Vector<n,T>()
		{
			for (int i = 0; i < n; ++i)
				vec[i] = 0;
		}

		~Vector<n,T>() {}

		/**
		 * Addresses the i-th element of the vector
		 * @param i the index of the element
		 */
		T at(int i);

 		/**
		 * Cross product with another vector
		 * @param v the other vector
		 */
		float cross(Vector<n,T> &v);

 		/**
		 * Dot product with another vector
		 * @param v the other vector
		 */
		float dot(Vector<n,T> &v);

 		/**
		 * Returns true if the vector is orthogonal to another given as an argument, false otherwise
		 * @param v the other vector
		 */
		bool is_ortho(Vector<n,T> &v);

 		/**
		 * Returns true if the vector contains an invalid value, false otherwise. Notably, if the vector contains nan as values
		 */
		bool is_null();

 		/**
		 * Returns the norm of the vector
		 */
		float norm();

 		/**
		 * Returns a copy of the vector normalised.
		 */
		Vector<n,T>& to_unit();

		friend ostream& operator<<(ostream &out, const Vector<n,T> v)
		{
			out << "(";
			for(int i = 0; i < n; ++i)
			{
				if(i == 0)
					out << v.vec[i];
				else
					out << ", " << v.vec[i];
			}
			out << ")";
			return out;
		}

		T& operator[](int i);

		Vector<n,T>& operator+(const Vector<n,T>& v);

		Vector<n,T>& operator+=(const Vector<n,T>& v);

		Vector<n,T>& operator-(const Vector<n,T>& v);

		Vector<n,T>& operator-=(const Vector<n,T>& v);

		Vector<n,T>& operator*(const Vector<n,T>& v);

		Vector<n,T>& operator*=(const Vector<n,T>& v);

	};
}

#endif