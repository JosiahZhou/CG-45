#include <stdio.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#endif
#ifdef WIN32
#include <windows.h>
#endif
#include <algorithm> 
#include "raytracing.h"


//temporary variables
//these are only used to illustrate 
//a simple debug drawing. A ray 
Vec3Df testRayOrigin;
Vec3Df testRayDestination;

std::vector<AABB> boxes;

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should, if they 
	//are exported as WavefrontOBJ.
	//PLEASE ADAPT THE LINE BELOW TO THE FULL PATH OF THE dodgeColorTest.obj
	//model, e.g., "C:/temp/myData/GraphicsIsFun/dodgeColorTest.obj", 
	//otherwise the application will not load properly
	MyMesh.loadMesh("dodgeColorTest.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);
}


//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	Vec3Df direction = dest - origin;
	for (int b = 0; b < boxes.size(); b++)
	{
		AABB box = boxes[b];
		Vec3Df pin, pout;
		if (rayIntersectionPointBox(origin, direction, box, pin, pout))
		{
			for (int i = 0; i < MyMesh.triangles.size(); i++)
			{
				Vec3Df pointOfIntersection;
				float distanceRay;
				Triangle triangle = MyMesh.triangles[i];
				if (rayIntersectionPointTriangle(origin, direction, triangle, pointOfIntersection, distanceRay))
				{
					return Vec3Df(dest[0], dest[1], dest[2]);
				}
			}
		}
	}
	return Vec3Df(0, 0, 0);
	//return Vec3Df(dest[0], dest[1], dest[2]);
}

/**
 * Moller and Trumbore algorithnm: point(u,v) = (1-u-v)*p0 + u*p1 + v*p2
 * returns true if the ray intersects with the triangle.
*/
bool rayIntersectionPointTriangle(Vec3Df rayOrigin, Vec3Df rayDirection, Triangle triangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection)
{
	Vec3Df edges[2];
	float t;
	float u;
	float v;

	edges[0] = MyMesh.vertices[triangle.v[1]].p - MyMesh.vertices[triangle.v[0]].p;
	edges[1] = MyMesh.vertices[triangle.v[2]].p - MyMesh.vertices[triangle.v[0]].p;
	Vec3Df n = Vec3Df::crossProduct(edges[0], edges[1]);
	n.normalize();

	Vec3Df T = rayOrigin - MyMesh.vertices[triangle.v[0]].p;

	Vec3Df P = Vec3Df::crossProduct(rayDirection, edges[1]);
	Vec3Df Q = Vec3Df::crossProduct(T, edges[0]);

	float invDet = 1.0 / Vec3Df::dotProduct(P, edges[0]); // 1 / P . E1

	t = invDet * Vec3Df::dotProduct(Q, edges[1]);
	u = invDet * Vec3Df::dotProduct(P, T);
	v = invDet * Vec3Df::dotProduct(Q, rayDirection);

	// check whether there is an intersection
	if (u < 0.0 || v < 0.0 || u + v > 1.0) return false; // u >= 0, v >= 0, u + v <= 1

														 //	pointOfIntersection = MyMesh.vertices[triangle.v[0]].p + u * edges[0] + v * edges[1];
	pointOfIntersection = (1 - u - v) * MyMesh.vertices[triangle.v[0]].p + u * MyMesh.vertices[triangle.v[1]].p + v * MyMesh.vertices[triangle.v[2]].p;
	distanceLightToIntersection = t;

	return true;
}

/**
* pin - point of intersection going *in* the box
* pout - point of intersection going *out* the box
* returns true if the ray intersects with the box
* Reference: https://www.youtube.com/watch?v=USjbg5QXk3g (Math for Game Developers - Bullet Collision (Vector/AABB Intersection))
*/
bool rayIntersectionPointBox(Vec3Df rayOrigin, Vec3Df rayDirection, AABB box, Vec3Df& pin, Vec3Df& pout)
{
	Vec3Df min = box.minmax_.first;
	Vec3Df max = box.minmax_.second;

	float tminX, tminY, tminZ, tmaxX, tmaxY, tmaxZ, tinX, tinY, tinZ, toutX, toutY, toutZ, tin, tout;

	tminX = (min.p[0] - rayOrigin.p[0]) / rayDirection.p[0]; // (x - Ox) / Dx
	tminY = (min.p[1] - rayOrigin.p[1]) / rayDirection.p[1]; // (y - Oy) / Dy
	tminZ = (min.p[2] - rayOrigin.p[2]) / rayDirection.p[2]; // (z - Oz) / Dz

	tmaxX = (max.p[0] - rayOrigin.p[0]) / rayDirection.p[0]; // (x - Ox) / Dx
	tmaxY = (max.p[1] - rayOrigin.p[1]) / rayDirection.p[1]; // (y - Oy) / Dy
	tmaxZ = (max.p[2] - rayOrigin.p[2]) / rayDirection.p[2]; // (z - Oz) / Dz

	tinX = std::min(tminX, tmaxX);
	toutX = std::max(tminX, tmaxX);
	tinY = std::min(tminY, tmaxY);
	toutY = std::max(tminY, tmaxY);
	tinZ = std::min(tminZ, tmaxZ);
	toutZ = std::max(tminZ, tmaxZ);

	tin = std::max(tinX, std::max(tinY, tinZ));
	tout = std::min(toutX, std::min(toutY, toutZ));

	if (tin > tout || tout < 0)
	{
		return false;
	}

	pin = rayOrigin + rayDirection * tin;
	pout = rayOrigin + rayDirection * tout;
	return true;
}

std::pair<Vec3Df, Vec3Df> getMinAndMaxVertex()
{
	Vec3Df min = Vec3Df(INT32_MAX, INT32_MAX, INT32_MAX);
	Vec3Df max = Vec3Df(INT32_MIN, INT32_MIN, INT32_MIN);

	for (int i = 0; i < MyMesh.vertices.size(); ++i)
	{
		min[0] = MyMesh.vertices[i].p[0] < min[0] ? MyMesh.vertices[i].p[0] : min[0];
		min[1] = MyMesh.vertices[i].p[1] < min[1] ? MyMesh.vertices[i].p[1] : min[1];
		min[2] = MyMesh.vertices[i].p[2] < min[2] ? MyMesh.vertices[i].p[2] : min[2];

		max[0] = MyMesh.vertices[i].p[0] > max[0] ? MyMesh.vertices[i].p[0] : max[0];
		max[1] = MyMesh.vertices[i].p[1] > max[1] ? MyMesh.vertices[i].p[1] : max[1];
		max[2] = MyMesh.vertices[i].p[2] > max[2] ? MyMesh.vertices[i].p[2] : max[2];
	}

	return std::pair<Vec3Df, Vec3Df>(min, max);
}

// https://www.opengl.org/wiki/Calculating_a_Surface_Normal
Vec3Df calculateSurfaceNormal(Triangle triangle)
{
	Vec3Df trianglev0 = MyMesh.vertices[triangle.v[0]].p;
	Vec3Df trianglev1 = MyMesh.vertices[triangle.v[1]].p;
	Vec3Df trianglev2 = MyMesh.vertices[triangle.v[2]].p;

	Vec3Df U = trianglev1 - trianglev0;
	Vec3Df V = trianglev2 - trianglev0;

	float Ux = U.p[0];
	float Uy = U.p[1];
	float Uz = U.p[2];

	float Vx = V.p[0];
	float Vy = V.p[1];
	float Vz = V.p[2];

	float normalX = Uy * Vz - Uz * Vy;
	float normalY = Uz * Vx - Ux * Vz;
	float normalZ = Ux * Vy - Uy * Vx;

	return Vec3Df(normalX, normalY, normalZ); //return the normal of the triangle.
}

Vec3Df calculateCentroid(const Triangle t)
{
	Vec3Df trianglev0 = MyMesh.vertices[t.v[0]].p;
	Vec3Df trianglev1 = MyMesh.vertices[t.v[1]].p;
	Vec3Df trianglev2 = MyMesh.vertices[t.v[2]].p;

	return (trianglev0 + trianglev1 + trianglev2) / 3.0f;
}

void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame

	//let's draw the mesh
	MyMesh.draw();

	//let's draw the lights in the scene as points
	glPushAttrib(GL_ALL_ATTRIB_BITS); //store all GL attributes
	glDisable(GL_LIGHTING);
	glColor3f(1, 1, 1);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i = 0; i < MyLightPositions.size(); ++i)
		glVertex3fv(MyLightPositions[i].pointer());
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state. 
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous 
	//state after the pop.


	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(0, 1, 1);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glColor3f(0, 0, 1);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();


	// Traverse the tree
	//tree.getValue().highlightBoxEdges();

	for (AABB box : boxes)
	{
		box.highlightBoxEdges();
	}


	//draw whatever else you want...
	////glutSolidSphere(1,10,10);
	////allows you to draw a sphere at the origin.
	////using a glTranslate, it can be shifted to whereever you want
	////if you produce a sphere renderer, this 
	////triangulated sphere is nice for the preview
}

// https://stackoverflow.com/questions/13484943/print-a-binary-tree-in-a-pretty-way
int rec[1000006];
// Prints the BoxTree in directory-format
void printTree(struct BoxTree* curr, int depth)
{
	int i;
	if (curr == NULL)return;
	printf("\t");
	for (i = 0; i < depth; i++)
		if (i == depth - 1)
			printf("%s---", rec[depth - 1] ? "|" : "|");
		else
			printf("%s   ", rec[i] ? "|" : "  ");
	printf("%d\n", curr->data.triangles.size());
	rec[depth] = 1;
	printTree(curr->left, depth + 1);
	rec[depth] = 0;
	printTree(curr->right, depth + 1);
}

// Recursively adds the nodes of the BoxTree to "boxes" whom elements will be drawn on the screen.
void addBoxes(struct BoxTree* curr)
{
	if (curr == NULL)return;
	boxes.push_back(curr->data);
	addBoxes(curr->left);
	addBoxes(curr->right);
}

void addLeavesOnly(struct BoxTree* curr)
{
	if (curr == NULL)return;
	if (curr->left == NULL && curr->right == NULL) boxes.push_back(curr->data);
	addLeavesOnly(curr->left);
	addLeavesOnly(curr->right);
}

//yourKeyboardFunc is used to deal with keyboard input.
//t is the character that was pressed
//x,y is the mouse position in pixels
//rayOrigin, rayDestination is the ray that is going in the view direction UNDERNEATH your mouse position.
//
//A few keys are already reserved: 
//'L' adds a light positioned at the camera location to the MyLightPositions vector
//'l' modifies the last added light to the current 
//    camera position (by default, there is only one light, so move it with l)
//    ATTENTION These lights do NOT affect the real-time rendering. 
//    You should use them for the raytracing.
//'r' calls the function performRaytracing on EVERY pixel, using the correct associated ray. 
//    It then stores the result in an image "result.bmp".
//    Initially, this function is fast (performRaytracing simply returns 
//    the target of the ray - see the code above), but once you replaced 
//    this function and raytracing is in place, it might take a 
//    while to complete...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination)
{

	//here, as an example, I use the ray to fill in the values for my upper global ray variable
	//I use these variables in the debugDraw function to draw the corresponding ray.
	//try it: Press a key, move the camera, see the ray that was launched as a line.
	testRayOrigin = rayOrigin;
	testRayDestination = rayDestination;

	Vec3Df normRayDirection = rayDestination - rayOrigin;
	normRayDirection.normalize();

	// do here, whatever you want with the keyboard input t.
	switch (t) {
	case 'i':  // check intersection with triangle
	{
		//here, as an example, I use the ray to fill in the values for my upper global ray variable
		//I use these variables in the debugDraw function to draw the corresponding ray.
		//try it: Press a key, move the camera, see the ray that was launched as a line.

		for (int i = 0; i < MyMesh.triangles.size(); i++)
		{
			Triangle triangle = MyMesh.triangles[i];

			Vec3Df pointOfIntersection;
			float distanceRay;
			Vec3Df color;
			Vec3Df hit;
			//				if (rayIntersectionPointPlane(rayOrigin, normRayDirection, triangle, pointOfIntersection))
			//				{
			//					std::cout << pointOfIntersection << std::endl;
			//				}

			if (rayIntersectionPointTriangle(rayOrigin, normRayDirection, triangle, pointOfIntersection, distanceRay))
			{
				std::cout << "Ray InterSects Triangle: " << pointOfIntersection << std::endl;
				//					Material mat = MyMesh.materials[MyMesh.triangleMaterials[i]];
				//					mat.set_Kd(0, 0, 0);
				//					std::cout << (mat.Kd()) << std::endl;
			}
		}
	}
	break;
	case 'd': // constructs axis aligned bounding boxes
	{
		std::pair<Vec3Df, Vec3Df> minMax = getMinAndMaxVertex();
		AABB aabb = AABB(minMax.first, minMax.second);
		Vec3Df pin, pout;

		BoxTree root = BoxTree(aabb);
		//root.splitMiddle(4000);
		root.splitAvg(4000);

		//addBoxes(&root);
		addLeavesOnly(&root);
		printTree(&root, 0);

		/*if (rayIntersectionPointBox(rayOrigin, normRayDirection, boxes[0], pin, pout))
		{
			std::cout << "Ray InterSects Box: \n" << pin << "\n" << pout << std::endl;
		}*/
	}
	break;
	default:
		break;
	}


	std::cout << t << " pressed! The mouse was in location " << x << "," << y << "!" << std::endl;
}



/**********************************************************************************************
**Axis-Aligned BoundingBox class
***********************************************************************************************/

AABB::AABB() 
{
	minmax_ = std::pair<Vec3Df, Vec3Df>(Vec3Df(), Vec3Df());
	vertices_ = std::vector<Vertex>();
	sides_ = std::vector<std::pair<Vec3Df, Vec3Df>>();
	triangles = std::vector<Triangle>();
}

AABB::AABB(const Vec3Df min, const Vec3Df max)
{
	minmax_ = std::pair<Vec3Df, Vec3Df>(min, max);
	vertices_ = std::vector<Vertex>();
	sides_ = std::vector<std::pair<Vec3Df, Vec3Df>>();
	triangles = std::vector<Triangle>();

	//	+------+      
	//  |`.    |`.    
	//  |  `+--+---+  
	//  |   |  |   |  
	//  X---+--+   |  
	//   `. |   `. |  
	//     `+------+ 
	vertices_.push_back(Vertex(Vec3Df(min[0], min[1], min[2]))); // 0

	//	+------+      
	//  |`.    |`.    
	//  |  `+--+---+  
	//  |   |  |   |  
	//  X---+--+   |  
	//   `. |   `. |  
	//     `X------+ 
	vertices_.push_back(Vertex(Vec3Df(min[0], min[1], max[2])));

	//	X------+      
	//  |`.    |`.    
	//  |  `+--+---+  
	//  |   |  |   |  
	//  X---+--+   |  
	//   `. |   `. |  
	//     `X------+ 
	vertices_.push_back(Vertex(Vec3Df(min[0], max[1], min[2])));

	//	X------+      
	//  |`.    |`.    
	//  |  `X--+---+  
	//  |   |  |   |  
	//  X---+--+   |  
	//   `. |   `. |  
	//     `X------+ 
	vertices_.push_back(Vertex(Vec3Df(min[0], max[1], max[2])));

	//	X------+      
	//  |`.    |`.    
	//  |  `X--+---+  
	//  |   |  |   |  
	//  X---+--X   |  
	//   `. |   `. |  
	//     `X------+ 
	vertices_.push_back(Vertex(Vec3Df(max[0], min[1], min[2])));

	//	X------X      
	//  |`.    |`.    
	//  |  `X--+---+  
	//  |   |  |   |  
	//  X---+--X   |  
	//   `. |   `. |  
	//     `X------+ 
	vertices_.push_back(Vertex(Vec3Df(max[0], max[1], min[2])));

	//	X------X      
	//  |`.    |`.    
	//  |  `X--+---+  
	//  |   |  |   |  
	//  X---+--X   |  
	//   `. |   `. |  
	//     `X------X 
	vertices_.push_back(Vertex(Vec3Df(max[0], min[1], max[2])));

	//	X------X      
	//  |`.    |`.    
	//  |  `X--+---X  
	//  |   |  |   |  
	//  X---+--X   |  
	//   `. |   `. |  
	//     `X------X 
	vertices_.push_back(Vertex(Vec3Df(max[0], max[1], max[2]))); // 7

	// calculating the sides
	// -X  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(min[0], min[1], min[2]), Vec3Df(min[0], max[1], max[2])));

	// +X  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(max[0], min[1], min[2]), Vec3Df(max[0], max[1], max[2])));

	// -Y  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(min[0], min[1], min[2]), Vec3Df(max[0], min[1], max[2])));

	// +Y  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(min[0], max[1], min[2]), Vec3Df(max[0], max[1], max[2])));

	// -Z  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(min[0], min[1], min[2]), Vec3Df(max[0], max[1], min[2])));

	// +Z  
	sides_.push_back(std::pair<Vec3Df, Vec3Df>(Vec3Df(min[0], min[1], max[2]), Vec3Df(max[0], max[1], max[2])));

	// initialize triangles
	for (int i = 0; i < MyMesh.triangles.size(); i++)
	{
		// WithinBox/WithinBoxFull: You Decide...
		if (withinBox(MyMesh.triangles[i]))
		{
			triangles.push_back(MyMesh.triangles[i]);
		}
	}
}

bool AABB::withinBox(const Triangle t)
{
	for (int i = 0; i < 3; i++)
	{
		// if atleast one vertex is in box then we consider the triangle inside of the box
		if (MyMesh.vertices[t.v[i]].p[0] >= vertices_[0].p[0] && MyMesh.vertices[t.v[i]].p[0] <= vertices_[7].p[0] &&
			MyMesh.vertices[t.v[i]].p[1] >= vertices_[0].p[1] && MyMesh.vertices[t.v[i]].p[1] <= vertices_[7].p[1] &&
			MyMesh.vertices[t.v[i]].p[2] >= vertices_[0].p[2] && MyMesh.vertices[t.v[i]].p[2] <= vertices_[7].p[2])
		{
			return true;
		}
	}

	// all vertices are outside the box
	return false;
}

bool AABB::withinBoxFull(const Triangle t)
{
	for (int i = 0; i < 3; i++)
	{
		// if vertex is not in box then we consider the triangle outside of the box
		if (!(MyMesh.vertices[t.v[i]].p[0] >= vertices_[0].p[0] && MyMesh.vertices[t.v[i]].p[0] <= vertices_[7].p[0] &&
			MyMesh.vertices[t.v[i]].p[1] >= vertices_[0].p[1] && MyMesh.vertices[t.v[i]].p[1] <= vertices_[7].p[1] &&
			MyMesh.vertices[t.v[i]].p[2] >= vertices_[0].p[2] && MyMesh.vertices[t.v[i]].p[2] <= vertices_[7].p[2]))
		{
			return false;
		}
	}

	// all vertices are in the box
	return true;
}

void AABB::highlightBoxEdges()
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

	glVertex3f(vertices_[0].p[0], vertices_[0].p[1], vertices_[0].p[2]);
	glVertex3f(vertices_[1].p[0], vertices_[1].p[1], vertices_[1].p[2]);

	glVertex3f(vertices_[0].p[0], vertices_[0].p[1], vertices_[0].p[2]);
	glVertex3f(vertices_[2].p[0], vertices_[2].p[1], vertices_[2].p[2]);

	glVertex3f(vertices_[0].p[0], vertices_[0].p[1], vertices_[0].p[2]);
	glVertex3f(vertices_[4].p[0], vertices_[4].p[1], vertices_[4].p[2]);

	glVertex3f(vertices_[7].p[0], vertices_[7].p[1], vertices_[7].p[2]);
	glVertex3f(vertices_[3].p[0], vertices_[3].p[1], vertices_[3].p[2]);

	glVertex3f(vertices_[7].p[0], vertices_[7].p[1], vertices_[7].p[2]);
	glVertex3f(vertices_[5].p[0], vertices_[5].p[1], vertices_[5].p[2]);

	glVertex3f(vertices_[7].p[0], vertices_[7].p[1], vertices_[7].p[2]);
	glVertex3f(vertices_[6].p[0], vertices_[6].p[1], vertices_[6].p[2]);

	glVertex3f(vertices_[6].p[0], vertices_[6].p[1], vertices_[6].p[2]);
	glVertex3f(vertices_[4].p[0], vertices_[4].p[1], vertices_[4].p[2]);

	glVertex3f(vertices_[5].p[0], vertices_[5].p[1], vertices_[5].p[2]);
	glVertex3f(vertices_[4].p[0], vertices_[4].p[1], vertices_[4].p[2]);

	glVertex3f(vertices_[5].p[0], vertices_[5].p[1], vertices_[5].p[2]);
	glVertex3f(vertices_[2].p[0], vertices_[2].p[1], vertices_[2].p[2]);

	glVertex3f(vertices_[3].p[0], vertices_[3].p[1], vertices_[3].p[2]);
	glVertex3f(vertices_[2].p[0], vertices_[2].p[1], vertices_[2].p[2]);

	glVertex3f(vertices_[3].p[0], vertices_[3].p[1], vertices_[3].p[2]);
	glVertex3f(vertices_[1].p[0], vertices_[1].p[1], vertices_[1].p[2]);

	glVertex3f(vertices_[6].p[0], vertices_[6].p[1], vertices_[6].p[2]);
	glVertex3f(vertices_[1].p[0], vertices_[1].p[1], vertices_[1].p[2]);

	glEnd();
}

/**********************************************************************************************
**Axis-Aligned BoundingBox Tree Class
***********************************************************************************************/
BoxTree::BoxTree(const AABB data) 
{
	this->data = data;
	this->left = NULL;
	this->right = NULL;
}

BoxTree::BoxTree(const AABB data, BoxTree *left, BoxTree *right) 
{
	this->data = data;
	this->left = left;
	this->right = right;
}

void BoxTree::splitMiddle(int minTriangles)
{
	// reduces the boxsize to 'fit' the object (i.e. reduce the size of the boundingbox to the minimum required)
	if (data.triangles.size() > 0) {
		Vec3Df min = MyMesh.vertices[data.triangles[0].v[0]].p;
		Vec3Df max = MyMesh.vertices[data.triangles[0].v[0]].p;
		for (int z = 0; z < data.triangles.size(); z++)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int x = 0; x < 3; x++)
				{
					if (MyMesh.vertices[data.triangles[z].v[y]].p[x] > max[x])
					{
						max[x] = MyMesh.vertices[data.triangles[z].v[y]].p[x];
					}
					if (MyMesh.vertices[data.triangles[z].v[y]].p[x] < min[x])
					{
						min[x] = MyMesh.vertices[data.triangles[z].v[y]].p[x];
					}
				}
			}
		}

		data = AABB(min, max);
	}

	if (data.triangles.size() < minTriangles)
	{
		return;
	}

	float edgeX = Vec3Df::squaredDistance(data.vertices_[4].p, data.vertices_[0].p);
	float edgeY = Vec3Df::squaredDistance(data.vertices_[2].p, data.vertices_[0].p);
	float edgeZ = Vec3Df::squaredDistance(data.vertices_[1].p, data.vertices_[0].p);

	if (edgeX > edgeY && edgeX > edgeZ)
	{
		//	+------+      
		//  |`.    |`.    
		//  |  `3--X---7  
		//  |   |  |   |  
		//  0---+--+   |  
		//   `. |   `. |  
		//     `+------+ 
		Vec3Df midPoint = (data.vertices_[3].p + data.vertices_[7].p) / 2.0f;
		AABB leftNode = AABB(data.vertices_[0].p, midPoint);
		left = new BoxTree(leftNode); // Beware: Usage of "new"

		//	+------+      
		//  |`.    |`.    
		//  |  `+--+---7  
		//  |   |  |   |  
		//  0---X--4   |  
		//   `. |   `. |  
		//     `+------+ 
		midPoint = (data.vertices_[0].p + data.vertices_[4].p) / 2.0f;
		AABB rightNode = AABB(midPoint, data.vertices_[7].p);
		right = new BoxTree(rightNode); // Beware: Usage of "new"
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
		Vec3Df midPoint = (data.vertices_[6].p + data.vertices_[7].p) / 2.0f;
		AABB leftNode = AABB(data.vertices_[0].p, midPoint);
		left = new BoxTree(leftNode); // Beware: Usage of "new"

		//	2------+      
		//  |`.    |`.    
		//  X  `+------7  
		//  |   |  |   |  
		//  0---+--+   |  
		//   `. |   `. |  
		//     `+------+ 
		midPoint = (data.vertices_[0].p + data.vertices_[2].p) / 2.0f;
		AABB rightNode = AABB(midPoint, data.vertices_[7].p);
		right = new BoxTree(rightNode); // Beware: Usage of "new"

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
		Vec3Df midPoint = (data.vertices_[5].p + data.vertices_[7].p) / 2.0f;
		AABB leftNode = AABB(data.vertices_[0].p, midPoint);
		left = new BoxTree(leftNode); // Beware: Usage of "new"

		//	+------+      
		//  |`.    |`.    
		//  |  `+--+---7  
		//  |   |  |   |  
		//  0---|--+   |  
		//   `X |   `. |  
		//     1+------+ 
		midPoint = (data.vertices_[0].p + data.vertices_[1].p) / 2.0f;
		AABB rightNode = AABB(midPoint, data.vertices_[7].p);
		right = new BoxTree(rightNode); // Beware: Usage of "new"

	}

	left->splitMiddle(minTriangles);
	right->splitMiddle(minTriangles);
}

void BoxTree::splitAvg(int minTriangles)
{
	// save min and max
	Vec3Df oldMin = Vec3Df(data.minmax_.first[0], data.minmax_.first[1], data.minmax_.first[2]);
	Vec3Df oldMax = Vec3Df(data.minmax_.second[0], data.minmax_.second[1], data.minmax_.second[2]);

	Vec3Df newMin = Vec3Df(data.minmax_.first[0], data.minmax_.first[1], data.minmax_.first[2]);
	Vec3Df newMax = Vec3Df(data.minmax_.second[0], data.minmax_.second[1], data.minmax_.second[2]);

	// reduce empty space of bounding box
	if (data.triangles.size() > 0)
	{
		Vec3Df min = MyMesh.vertices[data.triangles[0].v[0]].p;
		Vec3Df max = MyMesh.vertices[data.triangles[0].v[0]].p;
		for (int z = 0; z < data.triangles.size(); z++)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int x = 0; x < 3; x++)
				{
					if (data.withinBoxFull(data.triangles[z]))
					{
						if (MyMesh.vertices[data.triangles[z].v[y]].p[x] > max[x])
						{
							max[x] = MyMesh.vertices[data.triangles[z].v[y]].p[x];
						}
						if (MyMesh.vertices[data.triangles[z].v[y]].p[x] < min[x])
						{
							min[x] = MyMesh.vertices[data.triangles[z].v[y]].p[x];
						}
					}
				}
			}
		}
		data = AABB(min, max);
	}

	if (data.triangles.size() < minTriangles)
	{
		return;
	}

	float edgeX = Vec3Df::squaredDistance(data.vertices_[4].p, data.vertices_[0].p);
	float edgeY = Vec3Df::squaredDistance(data.vertices_[2].p, data.vertices_[0].p);
	float edgeZ = Vec3Df::squaredDistance(data.vertices_[1].p, data.vertices_[0].p);

	Vec3Df avg = Vec3Df(0.0f, 0.0f, 0.0f);

	//add all the vertices inside the boundingbox to avg
	for (int l = 0; l < data.triangles.size(); ++l)
	{
		for (int m = 0; m < 3; ++m)
		{
			avg += MyMesh.vertices[data.triangles[l].v[m]].p;
		}
	}

	avg /= ((float)data.triangles.size() * 3.0f);

	if (edgeX > edgeY && edgeX > edgeZ)
	{
		newMin[0] = avg[0];
		newMax[0] = avg[0];
	}
	else if (edgeY > edgeX && edgeY > edgeZ)
	{
		newMin[1] = avg[1];
		newMax[1] = avg[1];
	}
	else
	{
		newMin[2] = avg[2];
		newMax[2] = avg[2];
	}

	AABB leftNode = AABB(oldMin, newMax);
	AABB rightNode = AABB(newMin, oldMax);
	left = new BoxTree(leftNode); // Beware: Usage of "new"
	right = new BoxTree(rightNode); // Beware: Usage of "new"

	left->splitAvg(minTriangles);
	right->splitAvg(minTriangles);
}