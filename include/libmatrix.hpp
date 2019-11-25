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
		Vector()
		{
			for (int i = 0; i < n; ++i)
				vec[i] = 0;
		}

		~Vector() {}

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
		Vector cross(Vector<n,T> &v)
		{
			if(n < 3)
				throw runtime_error("Vectors have less than 3 dimensions");
			Vector<n,T> res;
			res[0] = vec[1]*v.vec[2] - vec[2]*v.vec[1];
			res[1] = vec[2]*v.vec[0] - vec[0]*v.vec[2];
			res[2] = vec[0]*v.vec[1] - vec[1]*v.vec[0];
			return res;
		}

		/**
		 * Dot product with another vector
		 * @param v the other vector
		 */
		float dot(Vector &v)
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
		bool is_ortho(Vector &v)
		{
			return dot(v) == 0;
		}

		/**
		 * Returns true if the vector contains an invalid value, false otherwise. Notably, if the vector contains nan as values
		 */
		bool is_null()
		{
			for(int i = 0; i < n; ++i)
				if(isnan(vec[i]))
					return true;
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
		 * Returns true if the vector is a unit vector (its norm == 1)
		 */ 
		bool is_unit()
		{
			return norm() == 1.0f;
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

		Vector operator+(const Vector& v)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] + v.vec[i];
			return res;
		}

		Vector& operator+=(const Vector& v)
		{
			for (int i = 0; i < n; ++i)
				vec[i] += v.vec[i];
			return *this;
		}

		Vector operator-(const Vector& v)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] - v.vec[i];
			return res;
		}

		Vector operator-()
		{
			Vector<n,T> res;
			for (int i = 0; i < n; i++)
				res[i] = -vec[i];
			return res;
		}

		Vector& operator-=(const Vector& v)
		{
			for (int i = 0; i < n; ++i)
				vec[i] -= v.vec[i];
			return *this;
		}
		
		Vector operator*(float f)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] * f;
			return res;
		}

		friend Vector<n,float> operator*(float f, const Vector v)
		{
			Vector<n,float> res;
			for (int i = 0; i < n; ++i)
			{
				res[i] = v.vec[i] * f;
			}
			return res;
		}

		Vector& operator*=(float f)
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

	template <int n, int m, class T>
	class Matrix
	{
	private:
		T matrix[n][m];

	public:
		Matrix()
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] = 0;	
		}

		~Matrix() {}

		T at(int i, int j)
		{
			if(i >= n || j >= m)
				throw out_of_range("Out of range exception : Matrix::at(" + to_string(i) + "," + to_string(j) + ")");
			return matrix[i][j];
		}

		Matrix inverse()
		{

		}

		bool is_null()
		{
			for(int i = 0; i < n; ++i)
				for(int j = 0; j < m; ++j)
					if(isnan(matrix[i][j]))
						return true;
			return false;
		}

		bool is_ortho()
		{
			return false;
		}

		Matrix<m,n,T> transpose()
		{
			Matrix<m,n,T> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[j][i] = matrix[i][j];
			return res;
		}

		friend ostream& operator<<(ostream &out, Matrix mat)
		{
			out << "[";
			for (int i = 0; i < n; ++i)
			{
				out << "[";
				for (int j = 0; j < m; ++j)
				{
					if(j!=0)
						out << ", ";
					out << mat.matrix[i][j];
				}
				out << "]";
				if(i != n - 1)
					out << "," << endl;
			}
			out << "]";
			return out;
		}

		T* operator[](int i)
		{
			return matrix[i];
		}

		Matrix operator+(Matrix mat)
		{
			Matrix<n,m,T> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = mat.matrix[i][j] + matrix[i][j];
			return res;
		}

		Matrix operator+=(Matrix mat)
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] += mat.matrix[i][j];
			return *this;
		}

		Matrix operator*(float f)
		{
			Matrix<n,m,float> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = matrix[i][j] * f;
			return res;
		}

		friend Matrix operator*(float f, Matrix mat)
		{
			Matrix<n,m,float> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = f * mat.matrix[i][j];
			return res;
		}

		Matrix operator*(Vector<m,T> v)
		{
			Matrix<n,1,T> res;
			for (int i = 0; i < n; ++i)
			{
				float tmp = 0;
				for (int j = 0; j < m; ++j)
					tmp += matrix[i][j] * v[j];
				res[i][0] = tmp;
			}
			return res;
		}

		friend Matrix operator*(Vector<m,T> v, Matrix mat)
		{
			Matrix<n,1,T> res;
			for (int i = 0; i < n; ++i)
			{
				float tmp = 0;
				for (int j = 0; j < m; ++j)
					tmp += mat.matrix[i][j] * v[j];
				res[i][0] = tmp;
			}
			return res;
		}

		Matrix operator*=(Matrix mat)
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] *= mat.matrix[i][j];
			return *this;
		}
	};

	typedef Vector<2,int> Vec2i;
	typedef Vector<3,int> Vec3i;
	typedef Vector<4,int> Vec4i;
	typedef Vector<2,float> Vec2r;
	typedef Vector<3,float> Vec3r;
	typedef Vector<4,float> Vec4r;
	typedef Matrix<4,4,float> Mat44r;
}

#endif