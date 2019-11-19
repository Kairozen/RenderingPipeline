#ifndef LIBMATRIX_HPP
#define LIBMATRIX_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include <exception>
#include <typeinfo>

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
		T at(int i)
		{
			
			if(i>=n)
				throw out_of_range("Out of range exception : Vector::at(" + to_string(i) + ")");
			return vec[i];
		}

		/**
		 * Cross product with another vector
		 * @param v the other vector
		 */
		Vector<n,float> cross(Vector<n,T> &v)
		{
			if(n < 3)
				throw runtime_error("Vectors have less than 3 dimensions");
			Vector<n,float> res;
			res[0] = vec[1]*v.vec[2] - vec[2]*v.vec[1];
			res[1] = vec[2]*v.vec[0] - vec[0]*v.vec[2];
			res[2] = vec[0]*v.vec[1] - vec[1]*v.vec[0];
			return res;
		}

		/**
		 * Dot product with another vector
		 * @param v the other vector
		 */
		float dot(Vector<n,T> &v)
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
		bool is_ortho(Vector<n,T> &v)
		{
			return dot(v) == 0;
		}

		/**
		 * Returns true if the vector contains an invalid value, false otherwise. Notably, if the vector contains nan as values
		 */
		bool is_null()
		{
			return false;
		}

		/**
		 * Returns the norm of the vector
		 */
		float norm()
		{
			float norm = 0.f;
			for (int i = 0; i < n; ++i)
				norm += vec[i] * vec[i];
			return sqrt(norm);
		}

		/**
		 * Returns a copy of the vector normalised.
		 */
		Vector<n,float> to_unit()
		{
			Vector<n,float> res;
			float norm = this->norm();
			for (int i = 0; i < n; ++i)
				res[i] = (1/norm)*vec[i];
			return res;
		}

		T& operator[](int i)
		{
			return vec[i];
		}

		Vector<n,T> operator+(const Vector<n,T>& v)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] + v.vec[i];
			return res;
		}

		Vector<n,T>& operator+=(const Vector<n,T>& v)
		{
			for (int i = 0; i < n; ++i)
				vec[i] += v.vec[i];
			return *this;
		}

		Vector<n,T> operator-(const Vector<n,T>& v)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] - v.vec[i];
			return res;
		}

		Vector<n,T> operator-()
		{
			Vector<n,T> res;
			for (int i = 0; i < n; i++)
				res[i] = -vec[i];
			return res;
		}

		Vector<n,T>& operator-=(const Vector<n,T>& v)
		{
			for (int i = 0; i < n; ++i)
				vec[i] -= v.vec[i];
			return *this;
		}
		
		Vector<n,T> operator*(float f)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] * f;
			return res;
		}

		friend Vector<n,float> operator*(float f, const Vector<n,T> v)
		{
			Vector<n,float> res;
			for (int i = 0; i < n; ++i)
			{
				res[i] = v.vec[i] * f;
			}
			return res;
		}

		Vector<n,T>& operator*=(float f)
		{
			for (int i = 0; i < n; ++i)
				vec[i] *= f;
			return *this;
		}

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
	};
}

#endif