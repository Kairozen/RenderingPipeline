#ifndef LIBGEOMETRY_HPP
#define LIBGEOMETRY_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include <exception>
#include <typeinfo>
#include <libmatrix.hpp>

using namespace std;
using namespace libmatrix;
namespace libgeometry
{
    const float epsilon = 0.000001;

    template <int N, class T>
    class Direction : public Vector<N,T>
    {
        public:
            Direction() {}
            Direction(Vector<N,T> dir) 
            {
                for (int i = 0; i < N; i++)
                    this->vec[i] = dir.at(i);
            }

            Direction<N+1,T> full_direction()
            {
                Direction<N+1,T> res;
                for (int i = 0; i < N; i++)
                    res[i] = this->vec[i];
                res[N] = 0;
                return res;
            }
    };

	template <class T>
	class Quaternion
	{
	private:
		float scalar;
		Vector<3,T> vector;
	public:
		Quaternion()
		{
			scalar = 0.f;
			vector = Vector<3,T>();
		}

		Quaternion(const float scalar, const Vector<3,T>& vec): scalar(scalar), vector(vec) {}

		Quaternion(const float angle, const Direction<3,float>& direction)
		{
            float angle_rad = (angle * M_PI) / 180.f;
            scalar = cos(angle_rad/2.f);

            Vector<3,float> directionVector = (direction.is_unit())?((Vector<3,float>)direction):(direction.to_unit());
            for (int i = 0; i < 3; i++)
                vector[i] = directionVector.at(i) * sin(angle_rad/2.f);
		}

		~Quaternion() {}

        float angle_rad() const
        {
            return (acos(scalar)*2.f);
        }

        float angle() const
        {
            return (angle_rad() * 180.f) / M_PI;
        }

		Quaternion conjugate() const
		{
            Quaternion<T> res;
            res.scalar = scalar;
            res.vector = -vector;
            return res;
		}

        float re() const
        {
            return scalar;
        }

        Vector<3,T> im() const
        {
            return vector;
        }

        Quaternion to_norm() const
        {
            Quaternion<T> res;
            float n = norm();
            if(n != 0)
            {
                n = 1/n;
                res.scalar = scalar * n;
                res.vector = vector * n;
                return res;
            }
            throw runtime_error("Division by zero while normalizing quaternion");
        }

		float norm() const
		{
			return sqrt(scalar*scalar + dot(vector,vector));
		}

        Quaternion inverse() const
        {
            Quaternion<T> res;
            float absv = norm();
            absv *= absv;
            absv = 1 / absv;
            Quaternion<T> conj = conjugate();
            res.scalar = conj.scalar * absv;
            res.vector = conj.vector * absv;
            return res;
        }

		friend Quaternion operator+(const float& s, const Quaternion& q)
		{
			return q + s;
		}

		Quaternion& operator+=(const float& s)
		{
			scalar += s;
			return *this;
		}

		Quaternion operator+(const float& s)
		{
			Quaternion<T> res;
			res.scalar = scalar + s;
			res.vector = vector;
			return res;
		}

		Quaternion& operator+=(const Quaternion& q)
		{
			scalar += q.scalar;
			vector += q.vector;
			return *this;
		}

		Quaternion operator+(const Quaternion& q)
		{
			Quaternion<T> res;
			res.scalar = scalar + q.scalar;
			res.vector = vector + q.vector;
			return res;
		}

        friend Quaternion operator-(const float& s, const Quaternion& q)
		{
			return q - s;
		}

		Quaternion& operator-=(const float& s)
		{
			scalar -= s;
			return *this;
		}

		Quaternion operator-(const float& s)
		{
			Quaternion<T> res;
			res.scalar = scalar - s;
			res.vector = vector;
			return res;
		}

		Quaternion& operator-=(const Quaternion& q)
		{
			scalar -= q.scalar;
			vector -= q.vector;
			return *this;
		}

		Quaternion operator-(const Quaternion& q)
		{
			Quaternion<T> res;
			res.scalar = scalar - q.scalar;
			res.vector = vector - q.vector;
			return res;
		}

        friend Quaternion operator*(const float& s, const Quaternion& q)
        {
            return q * s;
        }

        Quaternion& operator*(const float& s)
        {
            Quaternion<T> res;
            res.scalar = scalar * s;
            res.vector = vector * s;
            return res;
        }

        Quaternion& operator*=(const Quaternion& q)
        {
            Quaternion<T> tmp;
            tmp.scalar = scalar * q.scalar - vector.dot(q.vector);
            tmp.vector = q.vector * scalar + vector * q.scalar + vector.cross(q.vector);
            (*this) = tmp;
            return *this;
        }

        Quaternion operator*(const Quaternion& q)
        {
            Quaternion<T> res;
            res.scalar = scalar * q.scalar - vector.dot(q.vector);
            res.vector = q.vector * scalar + vector * q.scalar + vector.cross(q.vector);
            return res;
        }

		friend ostream& operator<<(ostream &out, const Quaternion<T> q)
		{
            out << q.scalar << " + " << q.vector;
			return out;
		}
    };

    class Plane;
    class Sphere;

    template <int N, class T>
    class Point : public Vector<N,T>
    {
        public:
            Point() {}

            Point(Vector<N,T> v)
            {
                for (int i = 0; i < N; i++)
                    this->vec[i] = v.at(i);                
            }

            bool behind(const Plane& plane) const;

            bool outside(const Sphere& sphere) const;

            Direction<N,T> length_to(const Point& p)
            {
                Direction<N,T> d = (Vec3r)p - (Vec3r)(*this);
                return d;
            }

            Point<N+1, T> full_point() const
            {
                Point<N+1,T> res;
                for (int i = 0; i < N; i++)
                    res[i] = this->vec[i];
                res[N] = 1;
                return res;
            }

            Point rotate(const Quaternion<T>& q)
            {
                Point res = 2.f * dot(q.im(),this->vec) * q.im() + (q.re() * q.re() - dot(q.im(), q.im())) * this->vec + 2.f * q.re() * cross(q.im(), this->vec);
                return res;
            }
    };
    
    class Plane : public Vec4r
    {
        public:
            Plane() {}

            Plane(const Vec4r& v)
            {
                for (int i = 0; i < 4; i++)
                    this->vec[i] = v.at(i);
            }
    };

    class Sphere
    {
        private:
            Point<3,float> center;
            float radius;

        public:
            Sphere(): center(Point<3,float>()), radius(0.f) {}
            Sphere(const Point<3,float>& center, float radius): center(center), radius(radius) {}

            Point<3,float> get_center() const
            {
                return center;
            }

            float get_radius() const
            {
                return radius;
            }

            bool is_null()
            {
                return center.is_null();
            }

            bool behind(const Plane& plane) const
            {
                //TODO
                return false;
            }

            friend ostream& operator<<(ostream& out, const Sphere& s)
            {
                out << "Sphere Center : " << s.center << " Radius : " << s.radius;
                return out;
            }
    };

    class Rectangle
    {
        private:
            Point<3,float> leftTop, rightBottom;
        public:
            Rectangle(Point<3,float> lt, Point<3,float> rb): leftTop(lt), rightBottom(rb) {}

            bool is_null()
            {
                return leftTop.is_null() || rightBottom.is_null();
            }

            friend ostream& operator<<(ostream& out, const Rectangle& r)
            {
                out << "Rectangle : LeftTop : " << r.leftTop << " RightBottom : " << r.rightBottom;
                return out;
            }

    };

    class Triangle
    {
        private:
            Point<3,float> a, b, c;
        public:
            Triangle(Point<3,float> a, Point<3,float> b, Point<3,float> c): a(a), b(b), c(c) {}

            Point<3,float> get_p0() const
            {
                return a;
            }
            Point<3,float> get_p1() const
            {
                return b;
            }
            Point<3,float> get_p2() const
            {
                return c;
            }

            float area()
            {
                Vec3r ab = b - a;
                Vec3r ac = c - a;
                return 0.5f * cross(ab,ac).norm();
            }

            bool is_null()
            {
                return a.is_null() || b.is_null() || c.is_null();
            }

            friend ostream& operator<<(ostream& out, const Triangle& t)
            {
                out << "Triangle : A : " << t.a << " B : " << t.b << " C : " << t.c;
                return out;
            }
        
    };

    class LineSegment
    {
        private:
            Point<3,float> begin, end;

        public:
            LineSegment(const Point<3,float>& a, const Point<3,float>& b)
            {
                begin = a;
                end = b;
            }

            Point<3,float> get_begin() const
            {
                return begin;
            }

            Point<3,float> get_end() const
            {
                return end;
            }

            bool is_null()
            {
                return begin.is_null() || end.is_null();
            }

            float inter_coef(const Plane& plane)
            {
                Vec4r segmentDirection = begin.length_to(end).full_direction();
                Vec4r fullBegin = begin.full_point();
                float num = dot(plane,fullBegin);
                float denum = dot(plane,segmentDirection);
                if(num == 0) 
                    return 0.f;
                if(denum == 0)
                    return NAN;
                return -(num/denum);
            }

            Point<3,float> inter(const Plane& plane)
            {
                Vec4r segmentDirection = begin.length_to(end).full_direction();
                Vec4r fullBegin = begin.full_point();
                float t = inter_coef(plane);
                if(t == 0)
                    return begin;
                if(isnan(t))
                {
                    return Point<3,float>(null_vector<3,float>());
                }
                Vec4r resVec = fullBegin + t * segmentDirection;
                Point<3,float> res;
                for (int i = 0; i < 3; i++)
                    res[i] = resVec[i];
                return res;
            }

            friend ostream& operator<<(ostream& out, const LineSegment& seg)
            {
                out << seg.begin << " - " << seg.end;
                return out;
            }
    };

    typedef enum {
        scaling, translation
    } transform_operation;

    class Transform
    {
        private:
            Mat44r transform_matrix;

        public:
            Transform(const Quaternion<float>& q)
            {
                float x,y,z,w;
                w = q.re();
                x = q.im().at(0);
                y = q.im().at(1);
                z = q.im().at(2);
                transform_matrix = Mat44r();
                transform_matrix[0][0] = 1 - (2 * y * y) - (2 * z * z);
                transform_matrix[0][1] = (2 * x * y) - (2 * w * z);
                transform_matrix[0][2] = (2 * x * z) + (2 * w * y);
                transform_matrix[1][0] = (2 * x * y) + (2 * w * z);
                transform_matrix[1][1] = 1 - (2 * x * x) - (2 * z * z);
                transform_matrix[1][2] = (2 * y * z) - (2 * w * x);
                transform_matrix[2][0] = (2 * x * z) - (2 * w * y);
                transform_matrix[2][1] = (2 * y * z) + (2 * w * x);
                transform_matrix[2][2] = 1 - (2 * x * x) - (2 * y * y);

                transform_matrix[3][3] = 1;
            }

            Transform(const float angle, const Direction<3,float>& axis)
            {
                Transform(Quaternion<float>(angle,axis));
            }

            Transform(const Vec3r& vect, transform_operation operation) 
            {
                if(operation == transform_operation::translation)
                {
                    transform_matrix = Identity44r;
                    for (int i = 0; i < 3; ++i)
                        transform_matrix[i][3] = vect.at(i);
                }
                else if(operation == transform_operation::scaling)
                {
                    transform_matrix = Identity44r;
                    for (int i = 0; i < 3; ++i)
                        transform_matrix[i][i] = vect.at(i);
                }
            }

            Transform(const Mat44r& m): transform_matrix(m) {}

            Mat44r get_matrix()
            {
                return transform_matrix;
            }

            Quaternion<float> to_quat() const
            {
                float scalar,x,y,z;
                float trace = 1 + transform_matrix.at(0,0) + transform_matrix.at(1,1) + transform_matrix.at(2,2);
                if(trace > epsilon) 
                {
                    float s = sqrt(trace) * 2;
                    x = (transform_matrix.at(2,1) - transform_matrix.at(1,2)) / s;
                    y = (transform_matrix.at(0,2) - transform_matrix.at(2,0)) / s;
                    z = (transform_matrix.at(1,0) - transform_matrix.at(0,1)) / s;
                    scalar = 0.25 * s;
                    return Quaternion<float>(scalar, Vec3r({x,y,z}));
                }
                else
                {
                    if(transform_matrix.at(0,0) > transform_matrix.at(1,1) && transform_matrix.at(0,0) > transform_matrix.at(2,2))
                    {
                        float s = sqrt(1 + transform_matrix.at(0,0) - transform_matrix.at(1,1) - transform_matrix.at(2,2)) * 2;
                        x = 0.25 * s;
                        y = (transform_matrix.at(1,0) + transform_matrix.at(0,1)) / s;
                        z = (transform_matrix.at(0,2) + transform_matrix.at(2,0)) / s;
                        scalar = (transform_matrix.at(2,1) + transform_matrix.at(1,2)) / s;
                        return Quaternion<float>(scalar, Vec3r({x,y,z}));
                    } 
                    else if(transform_matrix.at(1,1) > transform_matrix.at(2,2))
                    {
                        float s = sqrt(1 + transform_matrix.at(1,1) - transform_matrix.at(0,0) - transform_matrix.at(2,2)) * 2;
                        x = (transform_matrix.at(1,0) + transform_matrix.at(0,1)) / s;
                        y = 0.25 * s;
                        z = (transform_matrix.at(2,1) + transform_matrix.at(1,2)) / s;
                        scalar = (transform_matrix.at(0,2) + transform_matrix.at(2,0)) / s;
                        return Quaternion<float>(scalar, Vec3r({x,y,z}));
                    }
                    else
                    {
                        float s = sqrt(1 + transform_matrix.at(2,2) - transform_matrix.at(0,0) - transform_matrix.at(1,1)) * 2;
                        x = (transform_matrix.at(0,2) + transform_matrix.at(2,0)) / s;
                        y = (transform_matrix.at(2,1) + transform_matrix.at(1,2)) / s;
                        z = 0.25 * s;
                        scalar = (transform_matrix.at(1,0) + transform_matrix.at(0,1)) / s;
                        return Quaternion<float>(scalar, Vec3r({x,y,z}));
                    }   
                }
            }

            Transform concat(Transform t)
            {
                return Transform(t.transform_matrix * transform_matrix);
            }

            template <class T>
            Point<3,T> apply(const Point<3,T>& p)
            {
                Point<4,float> pt = apply(p.full_point());
                Point<3,float> res({pt.at(0), pt.at(1), pt.at(2)});
                return res;
            }

            template <class T>
            Direction<3,T> apply(const Direction<3,T>& d)
            {
                Direction<4,float> dir = apply(d.full_direction());
                Direction<3,float> res({dir.at(0), dir.at(1), dir.at(2)});
                return res;
            }

            template <class T>
            Point<4,T> apply(const Point<4,T>& p)
            {
                return Point<4,T>(transform_matrix * p); 
            }

            template <class T>
            Direction<4,T> apply(const Direction<4,T>& d)
            {
                return Point<4,T>(transform_matrix * d);
            }

            Triangle apply(const Triangle& t)
            {
                return Triangle(apply(t.get_p0()), apply(t.get_p1()), apply(t.get_p2()));
            }

            Sphere apply(const Sphere& s)
            {
                return Sphere(apply(s.get_center()), s.get_radius());
            }

            friend ostream& operator<<(ostream &out, const Transform t)
            {
                out << "Transform : " << endl << t.transform_matrix;
                return out;
            }
    };




    template <int N, class T>
    bool Point<N,T>::behind(const Plane& plane) const
    {
        if(N == 4)
            return dot(plane.to_unit(), Vec4r(this->vec)) < 0;
        else if(N == 3)
        {
            Vec4r point = Vec4r({this->vec[0],this->vec[1],this->vec[2], 1});
            return dot(plane.to_unit(), point) < 0;
        }
        return false;
    }

    template <int N, class T>
    bool Point<N,T>::outside(const Sphere& sphere) const
    {
        return length_to(sphere.get_center()).norm() > sphere.get_radius();
    }

    template <int N, class T>
    float distance(Point<N,T> a, Point<N,T> b)
    {
        float sum = 0.f;
        for(int i = 0; i < N; ++i)
            sum += (a.at(i) - b.at(i)) * (a.at(i) - b.at(i));
        return sqrt(sum);
    }

    template <class T>
    Point<3,T> point_to_3D_point(const Point<4,T>& point)
    {
        return Point<3,T>({point.at(0),point.at(1),point.at(2)});
    }
}
#endif