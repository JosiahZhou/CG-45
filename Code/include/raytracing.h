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
extern std::vector<int> MyLightPositionAmount;
extern std::vector<float> MyLightPositionRadius;
extern std::vector<float> MyLightPositionPower;
extern std::vector<std::vector<Vec3Df>> MySphereLightPositions;
extern Vec3Df MyCameraPosition; //currCamera
extern unsigned int WindowSize_X;//window resolution width
extern unsigned int WindowSize_Y;//window resolution height
extern unsigned int RayTracingResolutionX;  // largeur fenetre
extern unsigned int RayTracingResolutionY;  // largeur fenetre
extern int MyLightPositionsPointer;
extern void createLightPointer();

// Ray structure
struct Ray {
	Vec3Df origin;
	Vec3Df direction;
};

struct Intersection {
	Vec3Df point;
	Triangle triangle;
	Material material;
	float distance;
	float schlickCosTheta;
};

//use this function for any preprocessing of the mesh.
void init();

//you can use this function to transform a click to an origin and destination
//the last two values will be changed. There is no need to define this function.
//it is defined elsewhere
void produceRay(int x_I, int y_I, Vec3Df & origin, Vec3Df & dest);


//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest);

// Debug Ray
Vec3Df DebugRay(const Vec3Df & origin, const Vec3Df & dest, Triangle & t);
//Shadow test
bool isInShadow(Vec3Df & intersection, int & shadowpoints);
//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

// function to setip a sphere around a point (light position).
void setupMySphereLightPositions();

//want keyboard interaction? Here it is...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);

//intersection ray with [Blank]
bool rayIntersectionPointTriangle(Ray r, Triangle triangle, Triangle ignoreTriangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection);
bool rayIntersectionPointBox(Ray r, AABB box, Vec3Df& pin, Vec3Df& pout);

// Functions for recursive raytracing of reflection and refraction
void ComputeDirectLight(Intersection intersect, Vec3Df& directColor);
bool ComputeReflectedRay(Ray origRay, Vec3Df pointOfIntersection, Triangle triangleOfIntersection, Ray& reflectedRay);
bool ComputeRefractedRay(Ray origRay, Intersection intersect, Ray& refractedRay);
void Trace(unsigned int level, Ray ray, Vec3Df& color, Triangle ignoreTriangle);
void BounceLight(Ray ray, Vec3Df& luminance, Triangle ignoreTriangle);
bool Intersect(unsigned int level, const Ray ray, Vec3Df& pointOfIntersection, Triangle& triangleOfIntersection, Triangle ignoreTriangle, float& distance);
Vec3Df specularFunction(const Vec3Df &vertexPosition, Vec3Df &normal, Material material);
// Prints the BoxTree in directory-format
void printTree(struct BoxTree* curr, int depth);

// Recursively adds the nodes of the BoxTree to "boxes" whom elements will be drawn on the screen.
void showBoxes(struct BoxTree* curr); // adds all boxes
void showLeavesOnly(struct BoxTree* curr); // adds only leafs
void showIntersectionBoxOnly(Ray r, struct BoxTree* curr); // adds all intersected boxes
void showIntersectionLeafOnly(Ray r, struct BoxTree* curr); // adds only intersected leafs

// gets the minimum and maximum vertex of all triangles, used to create a bounding box
std::pair<Vec3Df, Vec3Df> getMinAndMaxVertex();

/**
* Gets the first box int the tree "curr", closest to the camera, that intersects with the ray.
* This Function should only be called within these pre-conditions:
* 1 - BoxTree *curr is not NULL
* 2 - The ray, atleast, intersects the overlapping boundingbox, or in other words the initial "curr->data". (So it should be wrapped inside an if statement).
**/
AABB getFirstIntersectedBox(Ray r, BoxTree* curr, Vec3Df& pin, Vec3Df& pout);
BoxTree* getFirstIntersectedBoxFast(Ray r, BoxTree* curr, Vec3Df& pin, Vec3Df& pout);

// gets all the intersected boxes, same 2 points for this class
void getAllIntersectedLeafs(Ray r, BoxTree* curr, Vec3Df& pin, Vec3Df& pout, std::vector<AABB> &intersections);

// calculate the intensity of light
float intensityOfLight(const float &distance, const float &power, const float &minimum);

Vec3Df calculateSurfaceNormal(Triangle triangle);
Vec3Df calculateCentroid(const Triangle t);

/**********************************************************************************************
**Axis-Aligned BoundingBox class
***********************************************************************************************/
class AABB {
public:
	AABB();
	//AABB(const Vec3Df min, const Vec3Df max);
	AABB(const Vec3Df min, const Vec3Df max);

	//Returns true if a triangle is partially inside the boundingbox
	bool withinBox(const Triangle t);

	//Returns true if a triangle is entirely inside the boundingbox
	bool withinBoxFull(const Triangle t);

	//highlights the box edges
	void highlightBoxEdges();

	//Trims the box to match the boundaries of the triangles
	AABB trim();

	// min and max value
	std::pair<Vec3Df, Vec3Df> minmax_;

	//vertices of the bounding box
	std::vector<Vertex> vertices_;

	//the sides of the box
	std::vector<std::pair<Vec3Df, Vec3Df>> sides_;

	//triangles residing inside the bounding box
	std::vector<Triangle> triangles;

	std::vector<unsigned int> materials;
};

/**********************************************************************************************
**Axis-Aligned BoundingBox Tree class
***********************************************************************************************/
class BoxTree {
public:
	BoxTree(const AABB data, int level);
	BoxTree(const AABB data, BoxTree *left, BoxTree *right, int level);

	// splits the box ("data") recursively into smaller parts, and adds them to the tree, until it the amount of triangles within the box is smaller than "minTriangles"
	void splitMiddle(int minTriangles, int maxLevel);
	void splitAvg(int minTriangles, int maxLevel, int currentEdge);

	AABB data;
	BoxTree *parent = NULL;
	BoxTree *left = NULL;
	BoxTree *right = NULL;
	int level = 0;
};
#endif
