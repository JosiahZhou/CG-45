#ifndef RAYTRACING_H_1
#define RAYTRACING_H_1
#include <vector>
#include "mesh.h"

//Welcome to your MAIN PROJECT...
//THIS IS THE MOST relevant code for you!
//this is an important file, raytracing.cpp is what you need to fill out
//In principle, you can do the entire project ONLY by working in these two files

class AABB;
extern Mesh MyMesh; //Main mesh
extern std::vector<Vec3Df> MyLightPositions;
extern Vec3Df MyCameraPosition; //currCamera
extern unsigned int WindowSize_X;//window resolution width
extern unsigned int WindowSize_Y;//window resolution height
extern unsigned int RayTracingResolutionX;  // largeur fenetre
extern unsigned int RayTracingResolutionY;  // largeur fenetre

//use this function for any preprocessing of the mesh.
void init();

//you can use this function to transform a click to an origin and destination
//the last two values will be changed. There is no need to define this function.
//it is defined elsewhere
void produceRay(int x_I, int y_I, Vec3Df & origin, Vec3Df & dest);


//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest);

// Debug Ray
Vec3Df DebugRay(const Vec3Df & origin, const Vec3Df & dest, Triangle t);
//Shadow test
bool isInShadow(Vec3Df & intersection, Triangle & triangle);
//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

//want keyboard interaction? Here it is...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);

//intersection ray with [Blank]
bool rayIntersectionPointTriangle(Vec3Df rayOrigin, Vec3Df rayDirection, Triangle triangle, Triangle ignoreTriangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection);
bool rayIntersectionPointBox(Vec3Df rayOrigin, Vec3Df rayDirection, AABB box, Vec3Df& pin, Vec3Df& pout);

Vec3Df calculateSurfaceNormal(Triangle triangle);
Vec3Df calculateCentroid(const Triangle t);

/**********************************************************************************************
**Axis-Aligned BoundingBox class
***********************************************************************************************/
class AABB {
public:
	AABB();
	AABB(const Vec3Df min, const Vec3Df max);

	//Returns true if a triangle is partially inside the boundingbox
	bool withinBox(const Triangle t);

	//Returns true if a triangle is entirely inside the boundingbox
	bool withinBoxFull(const Triangle t);

	//highlights the box edges
	void highlightBoxEdges();

	// min and max value
	std::pair<Vec3Df, Vec3Df> minmax_;

	//vertices of the bounding box
	std::vector<Vertex> vertices_;

	//the sides of the box
	std::vector<std::pair<Vec3Df, Vec3Df>> sides_;

	//triangles residing inside the bounding box
	std::vector<Triangle> triangles;
};


/**********************************************************************************************
**Axis-Aligned BoundingBox Tree class
***********************************************************************************************/
class BoxTree {
public:
	BoxTree(const AABB data);
	BoxTree(const AABB data, BoxTree *left, BoxTree *right);

	// splits the box ("data") recursively into smaller parts, and adds them to the tree, until it the amount of triangles within the box is smaller than "minTriangles"
	void splitMiddle(int minTriangles);
	void splitAvg(int minTriangles);

	AABB data;
	BoxTree *left = NULL;
	BoxTree *right = NULL;
};
#endif
