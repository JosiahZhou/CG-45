#ifndef RAYTRACING_H_1
#define RAYTRACING_H_1
#include <vector>
#include "mesh.h"
#include "tree.h"


//Welcome to your MAIN PROJECT...
//THIS IS THE MOST relevant code for you!
//this is an important file, raytracing.cpp is what you need to fill out
//In principle, you can do the entire project ONLY by working in these two files

extern Mesh MyMesh; //Main mesh
extern std::vector<Vec3Df> MyLightPositions;
extern std::vector<int> MyLightPositionAmount;
extern std::vector<float> MyLightPositionRadius;
extern std::vector<int> MyLightPositionPower;
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

void resetRecurseTestRays(void);

//use this function for any preprocessing of the mesh.
void init();

//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest);

// Debug Ray
Vec3Df DebugRay(const Vec3Df & origin, const Vec3Df & dest, Triangle t, Tree<Box>* root);

//Shadow test
bool isInShadow(Vec3Df & intersection, Triangle & intersectionTriangle, Tree<Box>* root);

// Functions for recursive raytracing of reflection and refraction
bool Intersect(unsigned int level, const Ray ray, Intersection& intersect, Triangle ignoreTriangle);
void Shade(unsigned int level, Ray origRay, Intersection intersect, Vec3Df& color);
void BounceLight(Ray ray, Vec3Df& luminance, Triangle ignoreTriangle);
void Trace(unsigned int level, Ray ray, Vec3Df& color, Triangle ignoreTriangle);
void ComputeDirectLight(Intersection intersect, Vec3Df& directColor);
bool ComputeRefractedRay(Ray origRay, Intersection intersect, Ray& refractedRay);
bool ComputeReflectedRay(Ray origRay, Vec3Df pointOfIntersection, Triangle triangleOfIntersection, Ray& reflectedRay);

//intersection ray with [Blank]
bool rayIntersectionPointTriangle(Ray r, const Triangle & triangle, const Triangle & ignoreTriangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection);
bool rayIntersectionPointBox(Ray r, Box box, Vec3Df& pin, Vec3Df& pout);

Vec3Df calculateSurfaceNormal(Triangle triangle);

void drawRay(Vec3Df rayOrig, Vec3Df rayDest, Vec3Df color);

Vec3Df calculateCentroid(const Triangle t);

//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

// function to setip a sphere around a point (light position).
void setupMySphereLightPositions();

// calculate the intensity of light
double intensityOfLight(const float &distance, const float &power, const float &minimum);

Vec3Df diffuseFunction(const Vec3Df &vertexPos, Vec3Df &normal, Material *material, Vec3Df lightPos);
Vec3Df specularFunction(const Vec3Df &vertexPosition, Vec3Df &normal, Material *material, Vec3Df lightPosition);
Vec3Df calculateShading(const Vec3Df &vertexPosition, Vec3Df &normal, Material *material);

std::vector<Box>* getIntersectionBoxes(Ray r, Tree<Box>* curr);
std::vector<Box>* getIntersectionLeaves(Ray r, Tree<Box>* curr);

/**
* Gets the first box int the tree "curr", closest to the camera, that intersects with the ray.
* This Function should only be called within these pre-conditions:
* 1 - BoxTree *curr is not NULL
* 2 - The ray, atleast, intersects the overlapping boundingbox, or in other words the initial "curr->data". (So it should be wrapped inside an if statement).
**/
Box getFirstIntersectedBox(Ray r, Tree<Box>* curr, Vec3Df& pin, Vec3Df& pout);
const Tree<Box>* getFirstIntersectedBoxFast(Ray r, const Tree<Box>* curr, Vec3Df& pin, Vec3Df& pout);

void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);
#endif