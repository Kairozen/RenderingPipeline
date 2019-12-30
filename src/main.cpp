#include <libgraphic.hpp>
#include <iostream>
#include <fstream>
#include <vector>
using namespace libgraphic;
using namespace std;

void load_geo_file(const string filename, Scene& scene, int num_obj)
{
	ifstream file(filename);
	if(file)
	{
		Object3D obj;
		obj.name = filename;
		obj.position = {num_obj * 0.3f, 0.f, 0.f};

		int nb_points, nb_triangles;
		file >> nb_points;
		for (int i = 0; i < nb_points; i++)
		{
			float c1,c2,c3;
			file >> c1 >> c2 >> c3;
			obj.add_vertex(Point<3,float>({c1,c2,c3}));
		}
		file >> nb_triangles;
		for (int i = 0; i < nb_triangles; i++)
		{
			int i1,i2,i3;
			file >> i1 >> i2 >> i3;
			obj.add_face(i1,i2,i3);
		}
		scene.add_object(obj);
	}
	else
	{
		throw runtime_error("Impossible to read file " + filename);
	}
}

int main(int argc, char const *argv[])
{
	Scene scene;
	for (int i = 0; i < argc-1; i++)
	{
		load_geo_file(argv[i+1], scene, i);
	}
	scene.start();	
	return 0;
}