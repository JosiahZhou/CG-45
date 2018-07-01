#include "box.h"
#include "raytracing.h"


Box::Box(void) {}

// Constructor for the Box class. The min and max of the mesh claculated.
// Triangles from mashes are added to the box which have all 
// of their vertices inside of the box.
Box::Box(const Mesh &mesh)
{
	Vec3Df newMin = Vec3Df(FLT_MAX, FLT_MAX, FLT_MAX);
	Vec3Df newMax = Vec3Df(FLT_MIN, FLT_MIN, FLT_MIN);

	for (int i = 0; i < mesh.vertices.size(); ++i)
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
		this->triangles.push_back(mesh.triangles[i]);
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

	for (int i = 0; i < mesh.triangles.size(); ++i)
	{
		if (contains(mesh.vertices[mesh.triangles[i].v[0]])
			&& contains(mesh.vertices[mesh.triangles[i].v[1]])
			&& contains(mesh.vertices[mesh.triangles[i].v[2]]))
		{
			this->triangles.push_back(mesh.triangles[i]);
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
		Vec3Df newMin = mesh.vertices[triangles[0].v[0]].p;
		Vec3Df newMax = mesh.vertices[triangles[0].v[0]].p;
		for (int i = 0; i < triangles.size(); ++i)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int x = 0; x < 3; x++)
				{
					if (withinBoxFull(triangles[i], mesh))
					{
						if (mesh.vertices[triangles[i].v[y]].p[x] > newMax[x])
						{
							newMax[x] = mesh.vertices[triangles[i].v[y]].p[x];
						}
						if (mesh.vertices[triangles[i].v[y]].p[x] < newMin[x])
						{
							newMin[x] = mesh.vertices[triangles[i].v[y]].p[x];
						}
					}
				}
			}
		}
		this->min = newMin;
		this->max = newMax;
	}
}

bool Box::withinBoxFull(const Triangle t, const Mesh &MyMesh)
{
	for (int i = 0; i < 3; i++)
	{
		// if vertex is not in box then we consider the triangle outside of the box
		if (!(MyMesh.vertices[t.v[i]].p[0] >= corners[0].p[0] && MyMesh.vertices[t.v[i]].p[0] <= corners[7].p[0] &&
			MyMesh.vertices[t.v[i]].p[1] >= corners[0].p[1] && MyMesh.vertices[t.v[i]].p[1] <= corners[7].p[1] &&
			MyMesh.vertices[t.v[i]].p[2] >= corners[0].p[2] && MyMesh.vertices[t.v[i]].p[2] <= corners[7].p[2]))
		{
			return false;
		}
	}

	// all vertices are in the box
	return true;
}

bool Box::withinBox(const Triangle t, const Mesh &MyMesh)
{
	for (int i = 0; i < 3; i++)
	{
		// if atleast one vertex is in box then we consider the triangle inside of the box
		if (!(MyMesh.vertices[t.v[i]].p[0] >= corners[0].p[0] && MyMesh.vertices[t.v[i]].p[0] <= corners[7].p[0] &&
			MyMesh.vertices[t.v[i]].p[1] >= corners[0].p[1] && MyMesh.vertices[t.v[i]].p[1] <= corners[7].p[1] &&
			MyMesh.vertices[t.v[i]].p[2] >= corners[0].p[2] && MyMesh.vertices[t.v[i]].p[2] <= corners[7].p[2]))
		{
			return true;
		}
	}

	// all vertices are outside the box
	return false;
}

bool Box::contains(const Triangle &t, const Mesh &mesh)
{
	if (contains(mesh.vertices[t.v[0]])
		&& contains(mesh.vertices[t.v[1]])
		&& contains(mesh.vertices[t.v[2]]))
	{
		return true;
	}
	return false;
}


// Returns two halves of a bounding box. The bounding box is first trimmed to
// into its smallest form. The average vertex position of all the triangles is
// calculated. The longest dimenention of the box is split in half and two
// box are made using the new min and max vertices.
std::pair<Box, Box> Box::split(Mesh &mesh)
{
	//float edgeX = Vec3Df::squaredDistance(corners[4].p, corners[0].p);
	//float edgeY = Vec3Df::squaredDistance(corners[2].p, corners[0].p);
	//float edgeZ = Vec3Df::squaredDistance(corners[1].p, corners[0].p);

	//Vec3Df avg = Vec3Df(0.0f, 0.0f, 0.0f);

	//for (int i = 0; i < triangles.size(); ++i)
	//{
	//	for (int m = 0; m < 3; ++m)
	//	{
	//		avg += mesh.vertices[triangles[i].v[m]].p;
	//	}
	//}

	//avg /= ((float)triangles.size() * 3.0f);

	//int edge = 2;
	//if (edgeX > edgeY && edgeX > edgeZ)
	//{
	//	edge = 0;
	//}
	//else if (edgeY > edgeX && edgeY > edgeZ)
	//{
	//	edge = 1;
	//}

	//Vec3Df oldMin = Vec3Df(min[0], min[1], min[2]);
	//Vec3Df oldMax = Vec3Df(max[0], max[1], max[2]);

	//Vec3Df newMin = Vec3Df(min[0], min[1], min[2]);
	//Vec3Df newMax = Vec3Df(max[0], max[1], max[2]);
	//newMin[edge] = avg[edge];
	//newMax[edge] = avg[edge];

	//Box leftNode = Box(oldMin, newMax, mesh);
	//Box rightNode = Box(newMin, oldMax, mesh);

	//for (int i = 0; i < triangles.size(); ++i)
	//{
	//	if (leftNode.contains(mesh.vertices[triangles[i].v[0]])
	//		|| leftNode.contains(mesh.vertices[triangles[i].v[1]])
	//		|| leftNode.contains(mesh.vertices[triangles[i].v[2]]))
	//	{
	//		if (rightNode.contains(mesh.vertices[triangles[i].v[0]])
	//			|| rightNode.contains(mesh.vertices[triangles[i].v[1]])
	//			|| rightNode.contains(mesh.vertices[triangles[i].v[2]]))
	//		{
	//			Triangle t = triangles[i];

	//			Vertex inver, outver1, outver2;
	//			int inverIdx, outver1Idx, outver2Idx;

	//			bool leftIsInver = false;

	//			if (leftNode.contains(mesh.vertices[t.v[0]]) && leftNode.contains(mesh.vertices[t.v[1]]))
	//			{
	//				inverIdx = t.v[2];
	//				outver1Idx = t.v[0];
	//				outver2Idx = t.v[1];
	//			}
	//			else if (leftNode.contains(mesh.vertices[t.v[0]]) && leftNode.contains(mesh.vertices[t.v[2]]))
	//			{
	//				inverIdx = t.v[1];
	//				outver1Idx = t.v[0];
	//				outver2Idx = t.v[2];
	//			}
	//			else if (leftNode.contains(mesh.vertices[t.v[1]]) && leftNode.contains(mesh.vertices[t.v[2]]))
	//			{
	//				inverIdx = t.v[0];
	//				outver1Idx = t.v[1];
	//				outver2Idx = t.v[2];
	//			}
	//			else if (leftNode.contains(mesh.vertices[t.v[0]]))
	//			{
	//				leftIsInver = true;
	//				inverIdx = t.v[0];
	//				outver1Idx = t.v[1];
	//				outver2Idx = t.v[2];
	//			}
	//			else if (leftNode.contains(mesh.vertices[t.v[1]]))
	//			{
	//				leftIsInver = true;
	//				inverIdx = t.v[1];
	//				outver1Idx = t.v[0];
	//				outver2Idx = t.v[2];
	//			}
	//			else
	//			{
	//				leftIsInver = true;
	//				inverIdx = t.v[2];
	//				outver1Idx = t.v[0];
	//				outver2Idx = t.v[1];
	//			}

	//			inver = mesh.vertices[inverIdx];
	//			outver1 = mesh.vertices[outver1Idx];
	//			outver2 = mesh.vertices[outver2Idx];

	//			Vertex pinver = Vertex(inver);
	//			Vertex poutver = Vertex(inver);

	//			Vec3Df pin, pout;
	//			Ray r;
	//			r.direction = inver.p - outver1.p;
	//			r.origin = inver.p;
	//			rayIntersectionPointBox(r, leftNode, pin, pout);

	//			Vec3Df pin2, pout2;
	//			r.direction = inver.p - outver2.p;
	//			r.origin = inver.p;
	//			rayIntersectionPointBox(r, leftNode, pin2, pout2);

	//			int size = mesh.vertices.size();

	//			pinver.p = pin;
	//			poutver.p = pin2;

	//			mesh.vertices.push_back(pinver);
	//			mesh.vertices.push_back(poutver);

	//			Triangle tri1 = Triangle(outver1Idx, t.t[0], size, t.t[1], outver2Idx, t.t[2]);
	//			Triangle tri2 = Triangle(outver2Idx, t.t[0], size, t.t[1], size + 1, t.t[2]);
	//			Triangle tri3 = Triangle(inverIdx, t.t[0], size, t.t[1], size + 1, t.t[2]);
	//			
	//			//error when pushed
	//			//mesh.triangles.push_back(tri1);
	//			//mesh.triangles.push_back(tri2);
	//			//mesh.triangles.push_back(tri3);

	//			if (leftIsInver)
	//			{
	//				leftNode.triangles.push_back(tri3);
	//				rightNode.triangles.push_back(tri1);
	//				rightNode.triangles.push_back(tri2);
	//			}
	//			else {
	//				rightNode.triangles.push_back(tri3);
	//				leftNode.triangles.push_back(tri1);
	//				leftNode.triangles.push_back(tri2);
	//			}
	//		}
	//	}
	//}

	//return std::pair<Box, Box>(leftNode, rightNode);
	float edgeX = Vec3Df::squaredDistance(corners[4].p, corners[0].p);
	float edgeY = Vec3Df::squaredDistance(corners[2].p, corners[0].p);
	float edgeZ = Vec3Df::squaredDistance(corners[1].p, corners[0].p);

	Vec3Df avg = Vec3Df(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < triangles.size(); ++i)
	{
		for (int m = 0; m < 3; ++m)
		{
			avg += mesh.vertices[triangles[i].v[m]].p;
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
	std::cout << "tr: " << triangles.size() << std::endl;

	Box leftBox = Box(min, newMax, mesh);
	Box rightBox = Box(newMin, max, mesh);

	for (int i = 0; i < triangles.size(); i++)
	{
		if (!leftBox.contains(triangles[i], mesh) && !rightBox.contains(triangles[i], mesh))
		{
			int inVer, outVer1, outVer2;
			if (leftBox.contains(mesh.vertices[triangles[i].v[0]]))
			{
				if (leftBox.contains(mesh.vertices[triangles[i].v[1]]))
				{
					inVer = triangles[i].v[2];
					outVer1 = triangles[i].v[0];
					outVer2 = triangles[i].v[1];
				}
				else
				{
					if (leftBox.contains(mesh.vertices[triangles[i].v[2]]))
					{
						inVer = triangles[i].v[1];
						outVer1 = triangles[i].v[0];
						outVer2 = triangles[i].v[2];
					}
					else
					{
						inVer = triangles[i].v[0];
						outVer1 = triangles[i].v[1];
						outVer2 = triangles[i].v[2];
					}
				}
			}
			else
			{
				if (leftBox.contains(mesh.vertices[triangles[i].v[1]]))
				{
					if (leftBox.contains(mesh.vertices[triangles[i].v[2]]))
					{
						inVer = triangles[i].v[0];
						outVer1 = triangles[i].v[1];
						outVer2 = triangles[i].v[2];
					}
					else
					{
						inVer = triangles[i].v[1];
						outVer1 = triangles[i].v[0];
						outVer2 = triangles[i].v[2];
					}
				}
				else
				{
					inVer = triangles[i].v[2];
					outVer1 = triangles[i].v[0];
					outVer2 = triangles[i].v[1];
				}
			}
			Vec3Df IntersectionPoint1 = mesh.vertices[inVer].p + (mesh.vertices[inVer].p - mesh.vertices[outVer1].p) * ((newMax[edge] - mesh.vertices[inVer].p[edge]) / (mesh.vertices[inVer].p[edge] - mesh.vertices[outVer1].p[edge]));
			Vec3Df IntersectionPoint2 = mesh.vertices[inVer].p + (mesh.vertices[inVer].p - mesh.vertices[outVer2].p) * ((newMax[edge] - mesh.vertices[inVer].p[edge]) / (mesh.vertices[inVer].p[edge] - mesh.vertices[outVer2].p[edge]));
			std::cout << "Inter1: " << IntersectionPoint1 << std::endl;
			std::cout << "Inter2: " << IntersectionPoint2 << std::endl;
			if (IntersectionPoint1[0] < -20)
			{
				std::cout << "stop" << std::endl;
			}
			Vertex IntersectionVertex1 = mesh.vertices[inVer];
			Vertex IntersectionVertex2 = mesh.vertices[inVer];

			IntersectionVertex1.p = IntersectionPoint1;
			IntersectionVertex2.p = IntersectionPoint2;

			Triangle inTri = Triangle(inVer, triangles[i].t[0], mesh.vertices.size(), triangles[i].t[1], mesh.vertices.size() + 1, triangles[i].t[2]);
			Triangle outTri1 = Triangle(outVer1, triangles[i].t[0], outVer2, triangles[i].t[1], mesh.vertices.size(), triangles[i].t[2]);
			Triangle outTri2 = Triangle(outVer2, triangles[i].t[0], mesh.vertices.size(), triangles[i].t[1], mesh.vertices.size() + 1, triangles[i].t[2]);

			mesh.vertices.push_back(IntersectionVertex1);
			mesh.vertices.push_back(IntersectionVertex2);

			mesh.triangles.push_back(inTri);
			mesh.triangles.push_back(outTri1);
			mesh.triangles.push_back(outTri2);

			mesh.triangleMaterials.push_back(mesh.triangleMaterials[i]);
			mesh.triangleMaterials.push_back(mesh.triangleMaterials[i]);
			mesh.triangleMaterials.push_back(mesh.triangleMaterials[i]);

			if (leftBox.contains(inTri, mesh))
			{
				leftBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 3]);
			}
			else {
				rightBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 3]);
			}
			if (leftBox.contains(outTri1, mesh))
			{
				leftBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 2]);
			}
			else {
				rightBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 2]);
			}
			if (leftBox.contains(outTri2, mesh))
			{
				leftBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 1]);
			}
			else {
				rightBox.triangles.push_back(mesh.triangles[mesh.triangles.size() - 1]);
			}
		}
	}
	return std::pair<Box, Box>(leftBox, rightBox);

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