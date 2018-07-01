#ifndef BOX_H
#define BOX_H
#include "mesh.h"
#include <GL/glut.h>
#include <algorithm>


class Box
{
public:
	Box(void);
	Box(const Mesh &mesh);
	Box(const Vec3Df min, const Vec3Df max, std::vector<std::shared_ptr<const Triangle>>, Mesh &mesh);

	Vec3Df min;
	Vec3Df max;
	Vertex corners[8];
	std::vector<std::shared_ptr<const Triangle>> triangles;

	bool contains(const Vertex &v);

	bool contains(const Triangle &t, const Mesh &mesh);

	void highlightEdges();
	
	void trim(const Mesh &mesh);

	std::pair<Box, Box> splitAvg(Mesh &mesh);

	std::pair<Box, Box> splitMiddle(const int & minTriangles, Mesh &mesh);

	void clip(Box &leftBox, Box &rightBox, Mesh &mesh);
};

#endif