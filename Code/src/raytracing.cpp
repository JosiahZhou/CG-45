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
#include <random>
#include <chrono>

int light_speed_sphere = 12391;

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
    //MyLightPositions.push_back(MyCameraPosition);
    createLightPointer();
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

/**
 * Creates an area of lights (in the form of a sphere) around the light points.
 * This is necessary for soft shadows. Soft shadows is basically the same principle as hard shadows but then with multiple light sources.
 */
void setupMySphereLightPositions() {
    
    // Clear old light positions of the sphere.
    MySphereLightPositions.clear();
    
    // Loop through all the light centers.
    for (int i = 0; i < MyLightPositions.size(); i++) {
        
        // We use a seed so that every scene will be the same for all lights,
        // even though we are using random points.
        std::mt19937 seed(light_speed_sphere);
        std::uniform_real_distribution<double> rndFloat(0.0, 1.0);
        
        // Retrieve the values of the current light.
        Vec3Df lightPosition = MyLightPositions[i];
        float lightSphereWidth = MyLightPositionRadius[i];
        int lightSphereAmount = MyLightPositionAmount[i];
        
        // Create the list of points.
        std::vector<Vec3Df> currentLightSphere;
        
        // We only calculate the lightSphere if it is actually needed, else we just use the MyLightPositions.
        if (MyLightPositionAmount[i] > 1 && MyLightPositionRadius[i] > 0) {
            
            // Calculate position for every surface light.
            for (int i = 0; i < lightSphereAmount; i++) {
                double theta = 2 * M_PI * rndFloat(seed);
                double phi = acos(1 - 2 * rndFloat(seed));
                double x = lightPosition[0] + sin(phi) * cos(theta) * lightSphereWidth;
                double y = lightPosition[1] + sin(phi) * sin(theta) * lightSphereWidth;
                double z = lightPosition[2] + cos(phi) * lightSphereWidth;
                Vec3Df offset = Vec3Df(x, y, z);        
                currentLightSphere.push_back(offset);
            }
        } else {
            // We just add the normal light position.
            currentLightSphere.push_back(lightPosition);
        }
        
        // We add list of points around the sphere into the list.
        MySphereLightPositions.push_back(currentLightSphere);
    }
}

/**
 * Returns the intensity of the light.
 *
 * @param distance the distance between the object and the light.
 * @param power the power of the light.
 * @param minimum the minimum intensity of the light.
 * @return the intensity of the light between the object and the light.
 */
double intensityOfLight(const float &distance, const float &power, const float &minimum) {
    double intensity = 1 / (4 * M_PI * distance * distance * (1 / power) + 1);
    if (intensity > minimum) {
        return intensity;
    } else {
        return minimum;
    }
}


struct BoxTree {
	AABB data;
	BoxTree *left;
	BoxTree *right;
};

// TODO: complete and bugfix this method
// Pre-Condition, initial tree.data needs to be set.
void split(int minTriangles, BoxTree& tree)
{

	// if this.value has not been set, then triangles.size == 0
	if (tree.data.triangles.size() < minTriangles)
	{
		return;
	}

	/*float edgeX = Vec3Df::squaredDistance(tree.data.vertices_[4].p, tree.data.vertices_[0].p);
	float edgeY = Vec3Df::squaredDistance(tree.data.vertices_[2].p, tree.data.vertices_[0].p);
	float edgeZ = Vec3Df::squaredDistance(tree.data.vertices_[1].p, tree.data.vertices_[0].p);*/

	//BoxTree *l = (struct BoxTree *) malloc(sizeof(struct BoxTree));
	//BoxTree *r = (struct BoxTree *) malloc(sizeof(struct BoxTree));

	//if (edgeX > edgeY && edgeX > edgeZ)
	//{

	//	+------+      
	//  |`.    |`.    
	//  |  `3--X---7  
	//  |   |  |   |  
	//  0---+--+   |  
	//   `. |   `. |  
	//     `+------+ 
	Vec3Df midPoint = (tree.data.vertices_[3].p + tree.data.vertices_[7].p) / 2.0f;
	AABB leftNode = AABB(tree.data.vertices_[0].p, midPoint);
	//l->data = leftNode;

	//l->data = 0; //AABB(value.vertices_[0].p, midPoint);

	//	+------+      
	//  |`.    |`.    
	//  |  `+--+---7  
	//  |   |  |   |  
	//  0---X--4   |  
	//   `. |   `. |  
	//     `+------+ 
	midPoint = (tree.data.vertices_[0].p + tree.data.vertices_[4].p) / 2.0f;
	AABB rightNode = AABB(midPoint, tree.data.vertices_[7].p);
	//r->data = rightNode;

	/*}
	else if (edgeY > edgeX && edgeY > edgeZ)
	{

	}
	else
	{

	}*/

	std::cout << "SPLIT: " << tree.data.triangles.size() << std::endl;
	std::cout << "SPLIT: " << leftNode.triangles.size() << std::endl;
	std::cout << "SPLIT: " << rightNode.triangles.size() << std::endl;

	BoxTree l = { leftNode, NULL, NULL };
	tree.left = &l;

	std::cout << "SPLIT: " << tree.left->data.triangles.size() << std::endl;
	//tree.left = l;
	//tree.right = r;

	
	//std::cout << "SPLIT: " << l->data.triangles.size() << std::endl;
	//std::cout << "SPLIT: " << tree.left->data.triangles.size() << std::endl;

	//std::cout << edgeX << ", " << edgeY << ", " << edgeZ << "   " << tree.data.triangles.size() << std::endl;

	//split(minTriangles, l);
	//split(minTriangles, r);
}

void traverse(BoxTree tree)
{
	std::cout << tree.data.triangles.size() << std::endl;

	if (tree.left != NULL)
	{
		traverse(*tree.left);
	}

	if (tree.right != NULL)
	{
		traverse(*tree.right);
	}
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

		//std::pair<Vec3Df, Vec3Df> minMax = getMinAndMaxVertex();
		//AABB aabb = AABB(minMax.first, minMax.second);
		//BoxTree root;
		//root.data = aabb;
		//split(4000, root);
		//std::cout << "KEYB: " << root.left->data.triangles.size() << std::endl;
		//traverse(root);
		//std::cout << root.data.triangles.size() << std::endl;
		//std::cout << (*root.left).data.triangles.size() << std::endl;
		//std::cout << (*root.right).data.triangles.size() << std::endl;

		std::pair<Vec3Df, Vec3Df> minMax = getMinAndMaxVertex();
		std::cout << MyMesh.triangles.size() << std::endl;
		AABB aabb = AABB(minMax.first, minMax.second);
		boxes.push_back(aabb);
		std::cout << aabb.triangles.size() << std::endl;
		Vec3Df pin, pout;

		if (rayIntersectionPointBox(rayOrigin, normRayDirection, boxes[0], pin, pout))
		{
			std::cout << "Ray InterSects Box: \n" << pin << "\n" << pout << std::endl;
		}
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

AABB::AABB() {
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
