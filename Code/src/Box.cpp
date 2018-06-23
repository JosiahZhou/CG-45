#include "box.h"


Box::Box(void) {}

// Constructor for the Box class. The min and max of the mesh claculated.
// Triangles from mashes are added to the box which have all 
// of their vertices inside of the box.
Box::Box(const Mesh &mesh)
{
	Vec3Df newMin = Vec3Df(INT32_MAX, INT32_MAX, INT32_MAX);
	Vec3Df newMax = Vec3Df(INT32_MIN, INT32_MIN, INT32_MIN);

	for (int i = 0; i < mesh.vertices.size(); i++)
	{
		newMin[0] = mesh.vertices[i].p[0] < newMin[0] ? mesh.vertices[i].p[0] : newMin[0];
		newMin[1] = mesh.vertices[i].p[1] < newMin[1] ? mesh.vertices[i].p[1] : newMin[1];
		newMin[2] = mesh.vertices[i].p[2] < newMin[2] ? mesh.vertices[i].p[2] : newMin[2];

		newMax[0] = mesh.vertices[i].p[0] > newMax[0] ? mesh.vertices[i].p[0] : newMax[0];
		newMax[1] = mesh.vertices[i].p[1] > newMax[1] ? mesh.vertices[i].p[1] : newMax[1];
		newMax[2] = mesh.vertices[i].p[2] > newMax[2] ? mesh.vertices[i].p[2] : newMax[2];
	}

	this->min = newMin;
	this->max = newMax;

	//	+------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `+------+
	this->corners.push_back(Vertex(Vec3Df(min[0], min[1], min[2])));

	//	+------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], min[1], max[2])));

	//	X------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], max[1], min[2])));

	//	X------+
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], max[1], max[2])));

	//	X------+
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(max[0], min[1], min[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(max[0], max[1], min[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------X
	this->corners.push_back(Vertex(Vec3Df(max[0], min[1], max[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---X
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------X
	this->corners.push_back(Vertex(Vec3Df(max[0], max[1], max[2])));

	for (int i = 0; i < mesh.triangles.size(); ++i)
	{
		this->triangles.push_back(&mesh.triangles[i]);
	}
}

// Constructor for the Box class. The min and max vectors form a
// box. Triangles from mashes are added to the box which have all 
// of their vertices inside of the box.
Box::Box(const Vec3Df min, const Vec3Df max, const Mesh &mesh)
{
	this->min = min;
	this->max = max;

	//	+------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `+------+
	this->corners.push_back(Vertex(Vec3Df(min[0], min[1], min[2])));

	//	+------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], min[1], max[2])));

	//	X------+
	//  |`.    |`.
	//  |  `+--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], max[1], min[2])));

	//	X------+
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--+   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(min[0], max[1], max[2])));

	//	X------+
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(max[0], min[1], min[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------+
	this->corners.push_back(Vertex(Vec3Df(max[0], max[1], min[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---+
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------X
	this->corners.push_back(Vertex(Vec3Df(max[0], min[1], max[2])));

	//	X------X
	//  |`.    |`.
	//  |  `X--+---X
	//  |   |  |   |
	//  X---+--X   |
	//   `. |   `. |
	//     `X------X
	this->corners.push_back(Vertex(Vec3Df(max[0], max[1], max[2])));

	for (int i = 0; i < mesh.vertices.size(); ++i)
	{
		if (contains(mesh.vertices[mesh.triangles[i].v[0]])
			&& contains(mesh.vertices[mesh.triangles[i].v[1]])
			&& contains(mesh.vertices[mesh.triangles[i].v[2]]))
		{
			this->triangles.push_back(&mesh.triangles[i]);
		}
	}
}

// Returns true when all the x,y,z coordiantes are all within the corners
// of the bounding box.
bool Box::contains(const Vertex &v)
{
	if (v.p[0] >= corners[0].p[0] && v.p[0] <= corners[7].p[0]
		&& v.p[1] >= corners[0].p[1] && v.p[1] <= corners[7].p[1]
		&& v.p[2] >= corners[0].p[2] && v.p[2] <= corners[7].p[2])
	{
		return true;
	}
	return false;
}

// Highlights the bounding box by drawing green lines on the edges of the 
// bounding box.
void Box::highlightEdges()
{
	glBegin(GL_LINES);
	glColor3f(0, 1, 0);

	// 0------2
	// |`.    |`.
	// |  `4--+---5
	// |   |  |   |
	// 1---+--3   |
	//  `. |   `. |
	//    `6------7

	glVertex3f(corners[0].p[0], corners[0].p[1], corners[0].p[2]);
	glVertex3f(corners[1].p[0], corners[1].p[1], corners[1].p[2]);

	glVertex3f(corners[0].p[0], corners[0].p[1], corners[0].p[2]);
	glVertex3f(corners[2].p[0], corners[2].p[1], corners[2].p[2]);

	glVertex3f(corners[0].p[0], corners[0].p[1], corners[0].p[2]);
	glVertex3f(corners[4].p[0], corners[4].p[1], corners[4].p[2]);

	glVertex3f(corners[7].p[0], corners[7].p[1], corners[7].p[2]);
	glVertex3f(corners[3].p[0], corners[3].p[1], corners[3].p[2]);

	glVertex3f(corners[7].p[0], corners[7].p[1], corners[7].p[2]);
	glVertex3f(corners[5].p[0], corners[5].p[1], corners[5].p[2]);

	glVertex3f(corners[7].p[0], corners[7].p[1], corners[7].p[2]);
	glVertex3f(corners[6].p[0], corners[6].p[1], corners[6].p[2]);

	glVertex3f(corners[6].p[0], corners[6].p[1], corners[6].p[2]);
	glVertex3f(corners[4].p[0], corners[4].p[1], corners[4].p[2]);

	glVertex3f(corners[5].p[0], corners[5].p[1], corners[5].p[2]);
	glVertex3f(corners[4].p[0], corners[4].p[1], corners[4].p[2]);

	glVertex3f(corners[5].p[0], corners[5].p[1], corners[5].p[2]);
	glVertex3f(corners[2].p[0], corners[2].p[1], corners[2].p[2]);

	glVertex3f(corners[3].p[0], corners[3].p[1], corners[3].p[2]);
	glVertex3f(corners[2].p[0], corners[2].p[1], corners[2].p[2]);

	glVertex3f(corners[3].p[0], corners[3].p[1], corners[3].p[2]);
	glVertex3f(corners[1].p[0], corners[1].p[1], corners[1].p[2]);

	glVertex3f(corners[6].p[0], corners[6].p[1], corners[6].p[2]);
	glVertex3f(corners[1].p[0], corners[1].p[1], corners[1].p[2]);

	glEnd();
}

// Reduces the size of the bounding box. The new minimum vertex is equal to the
// smallest x,y,z coordinates respectively of all the vertices of all the
// triangles inside the bounding box.
void Box::trim(const Mesh &mesh)
{
	if (triangles.size() > 0)
	{
		Vec3Df newMin, newMax;
		for (int i = 0; i < triangles.size(); ++i)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int x = 0; x < 3; x++)
				{
					if (mesh.vertices[triangles[i]->v[y]].p[x] > newMax[x])
					{
						newMax[x] = mesh.vertices[triangles[i]->v[y]].p[x];
					}
					if (mesh.vertices[triangles[i]->v[y]].p[x] < newMin[x])
					{
						newMin[x] = mesh.vertices[triangles[i]->v[y]].p[x];
					}
				}
			}
		}
		this->min = newMin;
		this->max = newMax;
	}
}

// Returns two halves of a bounding box. The bounding box is first trimmed to
// into its smallest form. The average vertex position of all the triangles is
// calculated. The longest dimenention of the box is split in half and two
// box are made using the new min and max vertices.
std::pair<Box, Box> Box::split(const Mesh &mesh)
{
	trim(mesh);

	float edgeX = Vec3Df::squaredDistance(corners[4].p, corners[0].p);
	float edgeY = Vec3Df::squaredDistance(corners[2].p, corners[0].p);
	float edgeZ = Vec3Df::squaredDistance(corners[1].p, corners[0].p);

	Vec3Df avg;

	for (int i = 0; i < triangles.size(); ++i)
	{
		for (int m = 0; m < 3; ++m)
		{
			avg += mesh.vertices[triangles[i]->v[m]].p;
		}
	}

	avg /= ((float)triangles.size() * 3.0f);

	int edge = 2;
	if (edgeX > edgeY && edgeX > edgeZ)
	{
		edge = 0;
	}
	else if (edgeY > edgeX && edgeY > edgeZ)
	{
		edge = 1;
	}

	Vec3Df newMin = min;
	Vec3Df newMax = max;
	newMin[edge] = avg[edge];
	newMax[edge] = avg[edge];

	return std::pair<Box, Box>(Box(min, newMax, mesh), Box(newMin, max, mesh));
}

// Returns two halves of a bounding box. The bounding box is first trimmed to
// into its smallest form. The longest dimenention of the box is split in half and two
// box are made using the middle vertices.
std::pair<Box, Box> Box::splitMiddle(const int &minTriangles, const Mesh &mesh)
{
	trim(mesh);

	if (triangles.size() < minTriangles)
	{
		return std::pair<Box, Box>();
	}

	float edgeX = Vec3Df::squaredDistance(corners[4].p, corners[0].p);
	float edgeY = Vec3Df::squaredDistance(corners[2].p, corners[0].p);
	float edgeZ = Vec3Df::squaredDistance(corners[1].p, corners[0].p);

	Vec3Df newMin, newMax;;

	if (edgeX > edgeY && edgeX > edgeZ)
	{
		//	+------+
		//  |`.    |`.
		//  |  `3--X---7
		//  |   |  |   |
		//  0---+--+   |
		//   `. |   `. |
		//     `+------+
		newMin = (corners[3].p + corners[7].p) / 2.0f;

		//	+------+
		//  |`.    |`.
		//  |  `+--+---7
		//  |   |  |   |
		//  0---X--4   |
		//   `. |   `. |
		//     `+------+
		newMax = (corners[0].p + corners[4].p) / 2.0f;
	}
	else if (edgeY > edgeX && edgeY > edgeZ)
	{
		//	+------+
		//  |`.    |`.
		//  |  `+------7
		//  |   |  |   |
		//  0---+--+   X
		//   `. |   `. |
		//     `+------6
		newMin = (corners[6].p + corners[7].p) / 2.0f;

		//	2------+
		//  |`.    |`.
		//  X  `+------7
		//  |   |  |   |
		//  0---+--+   |
		//   `. |   `. |
		//     `+------+
		newMax = (corners[0].p + corners[2].p) / 2.0f;
	}
	else
	{
		//	+------5
		//  |`.    |`X
		//  |  `+------7
		//  |   |  |   |
		//  0---+--+   |
		//   `. |   `. |
		//     `+------+
		newMin = (corners[5].p + corners[7].p) / 2.0f;

		//	+------+
		//  |`.    |`.
		//  |  `+--+---7
		//  |   |  |   |
		//  0---|--+   |
		//   `X |   `. |
		//     1+------+
		newMax = (corners[0].p + corners[1].p) / 2.0f;
	}

	return std::pair<Box, Box>(Box(corners[0].p, newMax, mesh), Box(newMin, corners[7].p, mesh));
}
