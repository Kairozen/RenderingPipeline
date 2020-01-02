#ifndef LIBMATRIX_HPP
#define LIBMATRIX_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include <exception>
#include <typeinfo>
#include <initializer_list>
using namespace std;

namespace libmatrix
{
	const float epsilon = 0.000001;
	template <int n, class T>
	class Vector
	{
	protected:
		T vec[n];

	public:
		Vector()
		{
			for (int i = 0; i < n; ++i)
				vec[i] = 0;
		}

		Vector(const T array[n])
		{
			for (int i = 0; i < n; i++)
				vec[i] = array[i];
		}

		Vector(initializer_list<T> list)
		{
			int i = 0;
			for(typename initializer_list<T>::iterator it = list.begin(); it != list.end() && i < n; ++it)
				vec[i++] = *it;
		}

		~Vector() {}

		/**
		 * Addresses the i-th element of the vector
		 * @param i the index of the element
		 */
		T at(int i) const
		{
			
			if(i>=n)
				throw out_of_range("Out of range exception : Vector::at(" + to_string(i) + ")");
			return vec[i];
		}

		/**
		 * Cross product with another vector
		 * @param v the other vector
		 */
		Vector cross(const Vector<n,T> &v) const
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
		float dot(const Vector &v) const
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
		bool is_ortho(const Vector &v) const
		{
			return dot(v) == 0;
		}

		/**
		 * Returns true if the vector contains an invalid value, false otherwise. Notably, if the vector contains nan as values
		 */
		bool is_null() const
		{
			for(int i = 0; i < n; ++i)
				if(isnan(vec[i]))
					return true;
			return false;
		}

		/**
		 * Returns the norm of the vector
		 */
		float norm() const
		{
			float norm = 0.f;
			for (int i = 0; i < n; ++i)
				norm += vec[i] * vec[i];
			return sqrt(norm);
		}

		/**
		 * Returns true if the vector is a unit vector (its norm == 1)
		 */ 
		bool is_unit() const
		{
			return norm() == 1.0f;
		}

		/**
		 * Returns a copy of the vector normalised.
		 */
		Vector<n,float> to_unit() const
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

		Vector operator-(const Vector& v) const
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
				res[i] = vec[i] - v.vec[i];
			return res;
		}

		Vector operator-() const
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
		
		Vector operator*(float f) const
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

		Vector& operator*=(const float f)
		{
			for (int i = 0; i < n; ++i)
				vec[i] *= f;
			return *this;
		}

		friend ostream& operator<<(ostream &out, const Vector<n,T>& v)
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

	template <int n, class T>
	Vector<n,T> null_vector()
	{
		Vector<n,T> vec;
		for (int i = 0; i < n; i++)
			vec[i] = NAN;
		return vec;
	}

	template <int n, int m, class T>
	class Matrix
	{
	protected:
		T matrix[n][m];

		Matrix null_matrix()
		{
			Matrix<n,m,float> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = NAN;
			return res;
		}


	public:
		Matrix()
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] = 0;	
		}
		void swap_rows(int r1, int r2)
		{
			for (int i = 0; i < m; ++i)
			{
				T tmp = matrix[r1][i];
				matrix[r1][i] = matrix[r2][i];
				matrix[r2][i] = tmp; 
			}
		}

		~Matrix() {}

		T at (int i, int j) const
		{
			if(i >= n || j >= m)
				throw out_of_range("Out of range exception : Matrix::at(" + to_string(i) + "," + to_string(j) + ")");
			return matrix[i][j];
		}

		Matrix inverse()
		{
			if(n != m)
				return null_matrix();
			for (int i = 0; i < n; ++i)
			{
				bool col_of_zeros = true;
				for (int j = 0; j < n; ++j)
					if(matrix[i][j] != 0)
					{
						col_of_zeros = false;
						break;
					}
				if(col_of_zeros)
					return null_matrix();
			}
			for (int i = 0; i < n; ++i)
			{
				bool row_of_zeros = true;
				for (int j = 0; j < n; ++j)
					if(matrix[j][i] != 0)
					{
						row_of_zeros = false;
						break;
					}
				if(row_of_zeros)
					return null_matrix();
			}
			// Iniialisation of left and right matrix
			Matrix<n,n,T> left, right;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
				{
					left[i][j] = matrix[i][j];
					if(i==j)
						right[i][j] = 1;
				}
			// For each column of M'
			for (int col = 0; col < n; ++col)
			{
				// Find row with greatest non zero absolute value
				int row_max = col;
				for (int row = col; row < n; ++row)
				{
					if(abs(left.at(row,col)) > abs(left.at(row_max,col)) && left.at(row,col) != 0)
						row_max = row;
				}
				T max_value = left[row_max][col];
				// Swap rows i and j
				left.swap_rows(col,row_max);
				right.swap_rows(col,row_max);
				// Multiply row j by 1/M'ij
				for (int i = 0; i < n; ++i)
				{
					left[col][i] *= (1/max_value);
					right[col][i] *= (1/max_value); 
				}
				// For each row r != j, row r = row r + row j * -M'rj
				for (int r = 0; r < n; ++r)
				{
					T factor = - left[r][col];
					if(r != col)
					{
						for (int r_ind = 0; r_ind < n; ++r_ind)
						{
							left[r][r_ind] += left[col][r_ind] * factor;
							right[r][r_ind] += right[col][r_ind] * factor;
						}
					}
				}
			}
			return right;
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

		friend ostream& operator<<(ostream &out, const Matrix& mat)
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

		T* operator[](const int i)
		{
			return matrix[i];
		}

		Matrix operator+(const Matrix& mat) const
		{
			Matrix<n,m,T> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = mat.matrix[i][j] + matrix[i][j];
			return res;
		}

		Matrix& operator+=(const Matrix& mat)
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] += mat.matrix[i][j];
			return *this;
		}

		Matrix operator*(const float f) const
		{
			Matrix<n,m,float> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = matrix[i][j] * f;
			return res;
		}

		friend Matrix operator*(const float f, const Matrix& mat)
		{
			Matrix<n,m,float> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					res[i][j] = f * mat.matrix[i][j];
			return res;
		}

		Vector<n,T> operator*(const Vector<m,T>& v) const
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
			{
				float tmp = 0;
				for (int j = 0; j < m; ++j)
					tmp += matrix[i][j] * v.at(j);
				res[i] = tmp;
			}
			return res;
		}

		friend Vector<n,T> operator*(const Vector<m,T>& v, const Matrix& mat)
		{
			Vector<n,T> res;
			for (int i = 0; i < n; ++i)
			{
				float tmp = 0;
				for (int j = 0; j < m; ++j)
					tmp += mat.matrix[i][j] * v.at(j);
				res[i] = tmp;
			}
			return res;
		}

		template <int k>
		Matrix operator*(const Matrix<m,k,T>& mat) const
		{
			Matrix<n,k,T> res;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < k; ++j)
				{ 
					for (int z = 0; z < m; ++z)
						res[i][j] += matrix[i][z] * mat.matrix[z][j];
					if(abs(res[i][j]) < epsilon)
						res[i][j] = 0;
				}
			return res;
		}

		Matrix& operator*=(const Matrix& mat)
		{
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < m; ++j)
					matrix[i][j] *= mat.matrix[i][j];
			return *this;
		}
	};

	template<int n, class T>
	float dot(Vector<n,T> vec1, Vector<n,T> vec2)
	{
		return vec1.dot(vec2);
	}

	template<int n, class T>
	Vector<n,T> cross(Vector<n,T> vec1, Vector<n,T> vec2)
	{
		return vec1.cross(vec2);
	}


	typedef Vector<2,int> Vec2i;
	typedef Vector<3,int> Vec3i;
	typedef Vector<4,int> Vec4i;
	typedef Vector<2,float> Vec2r;
	typedef Vector<3,float> Vec3r;
	typedef Vector<4,float> Vec4r;
	typedef Matrix<4,4,float> Mat44r;

	const Vector<3,float> zerovector;
	const Vector<2,int> zerovec2i;
	const Vector<3,int> zerovec3i;
	const Vector<4,int> zerovec4i;
	const Vector<2,float> zerovec2r;
	const Vector<3,float> zerovec3r;
	const Vector<4,float> zerovec4r;

	template<int n, int m, class T>
	Matrix<n,m,T> identity()
	{
		Matrix<n,m,T> mat;
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < m; j++)
				if(i==j)
					mat[i][j] = 1;
		return mat;	
	}

	const Matrix<4,4,int> Identity44i = identity<4,4,int>();
	const Matrix<4,4,float> Identity44r = identity<4,4,float>();
}

#endif