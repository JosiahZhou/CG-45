#ifndef BOX_H
#define BOX_H
#include "mesh.h"
#include <GL/glut.h>


class Box
{
public:
	Box(void);
	Box(const Mesh &mesh);
	Box(const Vec3Df min, const Vec3Df max, const Mesh &mesh);

	Vec3Df min;
	Vec3Df max;
	std::vector<Vertex> corners;
	std::vector<const Triangle*> triangles;

	bool contains(const Vertex &v);

	void highlightEdges();
	
	void trim(const Mesh &mesh);

	std::pair<Box, Box> split(const Mesh &mesh);

	std::pair<Box, Box> splitMiddle(const int & minTriangles, const Mesh &mesh);

	bool withinBoxFull(const Triangle t, const Mesh &MyMesh);
};

#endif