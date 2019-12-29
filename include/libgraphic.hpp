#ifndef LIBGRAPHIC_HPP
#define LIBGRAPHIC_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include <exception>
#include <typeinfo>
#include <libmatrix.hpp>
#include <libgeometry.hpp>
#include <scene_interface.h>
#include <vector>
#include <gui.h>
using namespace std;
using namespace libmatrix;
using namespace libgeometry;
using namespace gui;

namespace libgraphic
{
    class Frustum
    {
        private:
            Plane near, far, right, left, top, bottom;
            float focal_distance, horizontal_size, vertical_size, near_dist, far_dist;
            void calculate_planes()
            {
                float aspectRatio = horizontal_size/vertical_size;
                near = Plane({0,0,-1,-near_dist});
                far = Plane({0,0,1,far_dist});
                left = Plane({focal_distance,0,-1,0});
                right = Plane({-focal_distance, 0,-1,0});
                bottom = Plane({0,focal_distance,-aspectRatio,0});
                top = Plane({0,-focal_distance,-aspectRatio,0});
            }
        public:
            
            Frustum() {}

            Frustum(float focal_distance, float h, float v, float n, float f): focal_distance(focal_distance), horizontal_size(h), vertical_size(v), near_dist(n), far_dist(f)
            {
                calculate_planes();
            }

            float get_far() const
            {
                return far_dist;
            }

            float get_near() const
            {
                return near_dist;
            }

            bool oustide(const Point<4,float>& point) const
            {
                return point.behind(near) ||
                    point.behind(far) ||
                    point.behind(left) ||
                    point.behind(right) ||
                    point.behind(bottom) ||
                    point.behind(top);
            }

            bool outside(const Sphere& sphere) const
            {
                return sphere.behind(near) ||
                    sphere.behind(far) ||
                    sphere.behind(left) ||
                    sphere.behind(right) ||
                    sphere.behind(bottom) ||
                    sphere.behind(top);
            }

            LineSegment inter(const LineSegment& segment) const
            {
                //TODO
                return segment;
            }

            void update(float h, float v, float e)
            {
                horizontal_size = h;
                vertical_size = v;
                focal_distance = e;
                calculate_planes();
            }
    };

    class Object3D
    {
        private:
            vector<Triangle> faces;
        public:
            vector<Point<3,float>> vertices;
            string name;
            Vec3r position;

            void add_vertex(const Point<3,float>& p)
            {
                vertices.push_back(p);
            }

            Sphere bsphere() const
            {
                Point<3,float> resCenter;
                float resRadius = MAXFLOAT;
                for(Point<3,float> vertex : vertices)
                    resCenter += vertex;
                resCenter *= (1/vertices.size());
                for(Point<3,float> vertex : vertices)
                {
                    float dist = distance(resCenter, vertex);
                    if(dist > resRadius)
                        resRadius = dist;
                }
                return Sphere(resCenter,resRadius);
            }

            Triangle face(unsigned int x) const
            {
                return faces[x];
            }

            unsigned int num_faces() const
            {
                return faces.size();
            }

            void add_face(unsigned int v1, unsigned int v2, unsigned int v3)
            {
                faces.push_back(Triangle(vertices[v1-1],vertices[v2-1],vertices[v3-1]));
            }

            void remove_face(unsigned int x)
            {
                if(x < faces.size())
                    faces.erase(faces.begin() + x);
            }
    };

    class Camera
    {
        public:
            int height, width;
            float move_speed = 0.001f, orientation_speed = 0.001f, zoom_speed = 0.1f;

            Point<3,float> current_position;
            Quaternion<float> current_orientation;
            Vec3r current_speed;
            Vec3r current_orientation_speed;
            float current_zoom_speed;
            float near, far, focal_distance, fov_angle;
            Frustum frustum;

            Camera() {}

            Camera(float a, float h, float v, float n, float f)
            {
                current_position = Point<3,float>({0.f,0.f,-1.f});
                current_orientation = Quaternion<float>(1.f,Vec3r({0.f,0.f,0.f}));
                width = h;
                height = v;
                near = n;
                far = f;
                fov_angle = a;
                focal_distance = 1 / tan(((a*M_PI)/180.f)/2.f);
                frustum = Frustum(focal_distance,h,v,n,f);
            }

            void turn_up()
            {
                current_orientation_speed[0] = orientation_speed;
            }
            void turn_down()
            {
                current_orientation_speed[0] = -orientation_speed;
            }
            void turn_left()
            {
                current_orientation_speed[1] = -orientation_speed;
            }
            void turn_right()
            {
                current_orientation_speed[1] = orientation_speed;
            }
            void move_up()
            {
                current_speed[1] = move_speed;
            }
            void move_down()
            {
                current_speed[1] = -move_speed;
            }
            void move_forward()
            {
                current_speed[2] = move_speed;
            }
            void move_backward()
            {
                current_speed[2] = -move_speed;
            }
            void move_left()
            {
                current_speed[0] = -move_speed;
            }
            void move_right()
            {
                current_speed[0] = move_speed;
            }
            void zoom()
            {
                current_zoom_speed = zoom_speed;
            }
            void unzoom()
            {
                current_zoom_speed = -zoom_speed;
            }
            void stop_move_up()
            {
                current_speed[1] = 0.f;
            }
            void stop_move_down()
            {
                current_speed[1] = 0.f;
            }
            void stop_move_left()
            {
                current_speed[0] = 0.f;   
            }
            void stop_move_right()
            {
                current_speed[0] = 0.f;
            }
            void stop_move_forward()
            {
                current_speed[2] = 0.f;   
            }
            void stop_move_backward()
            {
                current_speed[2] = 0.f;
            }
            void stop_turn_up()
            {
                current_orientation_speed[0] = 0.f;
            }
            void stop_turn_down()
            {
                current_orientation_speed[0] = 0.f;
            }
            void stop_turn_left()
            {
                current_orientation_speed[1] = 0.f;
            }
            void stop_turn_right()
            {
                current_orientation_speed[1] = 0.f;
            }
            void stop_zoom()
            {
                current_zoom_speed = 0.f;
            }
            void stop_unzoom()
            {
                current_zoom_speed = 0.f;
            }

            Transform projection_matrix() const
            {
                Mat44r m;                
                float left = -near / focal_distance;
                float right = near / focal_distance;
                float top = near / focal_distance;
                float bottom = -near / focal_distance;
                m[0][0] = (2 * near) / (right - left);
                m[0][2] = (right + left) / (right - left);
                m[1][1] = (2 * near) / (top - bottom);
                m[1][2] = (top + bottom) / (top - bottom);
                m[2][2] = -((far + near) / (far - near));
                m[2][3] = -((2 * near * far) / (far - near));
                m[3][2] = -1.f;
                //cout << m << endl;
                return Transform(m);
            }

            Transform view_matrix() const
            {
                Transform rotation = Transform(current_orientation);
                Transform translation = Transform(current_position, transform_operation::translation);
                return Transform((translation.concat(rotation).get_matrix()).inverse());
            }

            Transform get_transform() const
            {
                return Transform(current_position, transform_operation::translation).concat(Transform(current_orientation));
            }

            bool outside_frustum(const Sphere& sphere) const
            {
                return frustum.outside(sphere);
            }

            bool sees(const Triangle& triangle) const
            {
                //TODO
                return true;
            }

            LineSegment visible_part(const LineSegment& segment) const
            {
                return frustum.inter(segment);
            }

            void update()
            {
                current_position += current_speed;
                current_orientation = current_orientation * Quaternion<float>(1.f,current_orientation_speed);
                fov_angle += current_zoom_speed;
                if(fov_angle > 179.f)
                    fov_angle = 179.f;
                focal_distance = 1 / tan(((fov_angle*M_PI)/180.f)/2.f);
                frustum = Frustum(focal_distance,width,height,near,far);
            }
    };
    class Scene : public SceneInterface
    {
        private:
            Gui *gui;
            Camera camera;
            vector<Object3D> objects;

        public:
     
            Scene()
            {
                gui = new Gui();
                camera = Camera(90,50,50,0.01f,10);
            }

            ~Scene()
            {
                gui->stop();
            }

            void start()
            {
                gui->start();
                gui->main_loop(this);
            }

            void add_object(const Object3D& obj)
            {
                objects.push_back(obj);
            }

            void draw_object(Object3D * const obj) const
            {
                //cout << endl << "############# OBJECT ###############" << endl;
                for (unsigned int i = 0; i < obj->num_faces(); i++)
                {
                    Triangle t = obj->face(i);
                    //cout << "TRIANGLE " << i << " " << t << endl;
                    Triangle t_view = camera.view_matrix().apply(t);
                    if(camera.sees(camera.projection_matrix().apply(t_view))) 
                        draw_wire_triangle(t_view);
                }
            }

            void draw_wire_triangle(const Triangle& triangle) const
            {
                //cout << endl;
                draw_edge(triangle.get_p0(), triangle.get_p1());
                draw_edge(triangle.get_p0(), triangle.get_p2());
                draw_edge(triangle.get_p1(), triangle.get_p2());
            }

            void draw_edge(const Point<3,float>& a, const Point<3,float>& b) const
            {
                LineSegment visible = camera.visible_part(LineSegment(a, b));
                if(!visible.is_null())
                {
                    //cout << "Before : " << a << " " << b << endl;
                    Point<3,float> a_persp = perspective_projection(visible.get_begin().full_point());
                    Point<3,float> b_persp = perspective_projection(visible.get_end().full_point());
                    Point<2,float> a_final = Point<2,float>({a_persp.at(0),a_persp.at(1)});
                    Point<2,float> b_final = Point<2,float>({b_persp.at(0),b_persp.at(1)});
                    //cout << "Drawing " << a_final << " " << b_final << endl;
                    gui->render_line(a_final, b_final, white);
                }
            }

            Point<3,float> perspective_projection(const Point<4,float>& p) const
            {
                Point<4,float> projected = camera.projection_matrix().apply(p);
                //cout << "Matrix" << endl;
                //cout << camera.projection_matrix().get_matrix() << endl;
                //cout << "After " << projected << endl;
                if(abs(projected.at(3)) < epsilon)
                    projected[3] = 1.f;
                return Point<3,float>({projected.at(0)/projected.at(3),projected.at(1)/projected.at(3),projected.at(2)/projected.at(3)});
            }

            void draw() const
            {
                for(Object3D obj : objects)
                {
                    Sphere obj_proj = camera.view_matrix().concat(camera.projection_matrix()).apply(obj.bsphere());
                    if(!camera.outside_frustum(obj_proj))
                        draw_object(&obj);
                }
            }

            void press_up()
            {
                camera.turn_up();
            }
            void press_down()
            {
                camera.turn_down();
            }
            void press_left()
            {
                camera.turn_left();
            }
            void press_right()
            {
                camera.turn_right();
            }
            void press_space()
            {

            }
            void press_w()
            {
                camera.move_up();
            }
            void press_s()
            {
                camera.move_down();
            }
            void press_a()
            {
                camera.move_left();
            }
            void press_d()
            {
                camera.move_right();
            }
            void press_q()
            {
                camera.unzoom();
            }
            void press_e()
            {
                camera.zoom();
            }
            void press_z()
            {
                camera.move_forward();
            }
            void press_x()
            {
                camera.move_backward();
            }

            void release_updown()
            {
                camera.stop_turn_up();
                camera.stop_turn_down();
            }
            void release_leftright()
            {
                camera.stop_turn_left();
                camera.stop_turn_right();
            }
            void release_space()
            {

            }
            void release_ws()
            {
                camera.stop_move_down();
                camera.stop_move_up();
            }
            void release_ad()
            {
                camera.stop_move_left();
                camera.stop_move_right();
            }
            void release_qe()
            {
                camera.stop_zoom();
                camera.stop_unzoom();
            }
            void release_zx()
            {
                camera.stop_move_forward();
                camera.stop_move_backward();
            }

            void update()
            {
                camera.update();
                draw();
            }
    };
}
#endif