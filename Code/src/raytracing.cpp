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
#define _USE_MATH_DEFINES
#include <math.h>
#include <stack>

int light_speed_sphere = 12391;
bool diffuse = true;
bool specular = true;
bool ambient = true;
bool shading = true;

//temporary variables
//these are only used to illustrate
//a simple debug drawing. A ray
Vec3Df testRayOrigin;
Vec3Df testRayDestination;

Vec3Df recurseTestRayOrigins[50];
Vec3Df recurseTestRayDestinations[50];
bool drawRecurseRays;

unsigned int recurseTestRayCount;
Tree<Box>* root;


// std::map<std::string, unsigned int> materialIndex;

unsigned int maxRecursionLevel;
unsigned int shadowCounter;

void resetRecurseTestRays(void) {
	for (int i = 0; i < 50; i++) {
		recurseTestRayOrigins[i] = Vec3Df(0, 0, 0);
		recurseTestRayDestinations[i] = Vec3Df(0, 0, 0);
	}

	recurseTestRayCount = 0;
}

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
	MyMesh.loadMesh("test_scene_1.obj", true);
	MyMesh.computeVertexNormals();

	root = new Tree<Box>(Box(MyMesh));
	root->split(MyMesh.triangles.size() / 4.0, MyMesh);
	root->print(0);

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);

	maxRecursionLevel = 5;
	recurseTestRayCount = 0;

	resetRecurseTestRays();
	drawRecurseRays = false;
}


//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	Vec3Df direction = dest - origin;
	float minDist = INFINITY;
	Vec3Df foundIntersection;
	Triangle t;
	bool intersect;
	Ray r = { origin, direction };

	Vec3Df pin, pout;
	if (rayIntersectionPointBox(r, root->data, pin, pout))
	{
		//std::vector<AABB> intersections;
		//getAllIntersectedLeafs(r, &tree, pin, pout, intersections);
		//std::sort(intersections.begin(), intersections.end(), less_than());

		std::stack<Tree<Box>> s;
		s.push(*getFirstIntersectedBoxFast(r, root, pin, pout));
		while (!s.empty())
		{
			intersect = false;
			Tree<Box> Btree = s.top();
			Box box = Btree.data;
			s.pop();
			for (int i = 0; i < box.triangles.size(); i++)
			{
				Vec3Df pointOfIntersection;
				float distanceRay;
				const Triangle* triangle = box.triangles[i];
				if (rayIntersectionPointTriangle(r, *triangle, Triangle(), pointOfIntersection, distanceRay))
				{
					intersect = true;
					if (minDist > distanceRay && distanceRay > 0) {
						t = *triangle;
						foundIntersection = pointOfIntersection;
						minDist = distanceRay;
					}
				}
			}

			if (intersect) {
				if (isInShadow(foundIntersection, t, root)) {
					Vec3Df color = Vec3Df(1, 1, 1);
					int light = MyLightPositions.size() - shadowCounter;
					double weight = (double)shadowCounter / (double)MyLightPositions.size(); // percentage in shadow
					color = (1.0 - weight) * color;
					//std::cout << "[IsInShadow] : color: " << color << std::endl;
					return color; // shadow == black
				}
				else {
					// color and other stuff here as well...
					return Vec3Df(1, 1, 1); // light == white

				}
			}
			//else
			//{
			//	if (Btree.parent != NULL)
			//		s.push(*Btree.parent);
			//}
		}
	}
	//caclulate shadows --> only for the minimum distance ( closestIntersectionPoint)
	// color and other stuff here as well...
	return Vec3Df(1, 0, 0);
	//return Vec3Df(dest[0], dest[1], dest[2]);
}

Vec3Df DebugRay(const Vec3Df & origin, const Vec3Df & dest, Triangle t, Tree<Box>* root) {
	Vec3Df direction = dest - origin;
	float minDist = INFINITY;
	Vec3Df foundIntersection;
	std::stack<Tree<Box>*> stack;
	stack.push(root);
	Ray r = { origin, direction };

	while(!stack.empty())
	{
		Tree<Box>* curr = stack.top();
		Box box = stack.top()->data;
		stack.pop();
		Vec3Df pin, pout;
		if (rayIntersectionPointBox(r, box, pin, pout))
		{
			stack.push(curr->left);
			stack.push(curr->right);
			for (int i = 0; i < box.triangles.size(); ++i)
			{
				Vec3Df pointOfIntersection;
				float distanceRay;
				const Triangle* triangle = box.triangles[i];
				if (rayIntersectionPointTriangle(r, *triangle, Triangle(), pointOfIntersection, distanceRay))
				{
					if (minDist > distanceRay && distanceRay > 0) {
						t = *triangle;
						foundIntersection = pointOfIntersection;
						minDist = distanceRay;
					}
				}
			}
		}
	}
	return foundIntersection;
}

bool isInShadow(Vec3Df & intersection, Triangle & intersectionTriangle, Tree<Box>* root) {
	int counter = 0;
	bool shadow = false;
	int x = 0;
	for (x = 0; x < MyLightPositions.size(); x++) {
		Vec3Df origin = intersection;;
		Vec3Df dest = MyLightPositions[x];
		Triangle foundTriangle;
		float minDist = INFINITY;
		/*********************************************************/
		// Copied code from performRayTracing
		/*********************************************************/
		Vec3Df direction = dest - origin;
		origin = origin + 0.001f * direction;
		Ray ray;
		ray.origin = origin;
		ray.direction = direction;
		Vec3Df pin, pout;
		if (rayIntersectionPointBox(ray, root->data, pin, pout))
		{
			Box box = root->data;
			for (int i = 0; i < box.triangles.size(); i++) {
				Vec3Df intersect;
				float distanceRay;
				const Triangle* triangle = box.triangles[i];
				// if an intersection gets found, put the resulting point and triangle in the result vars
				if (rayIntersectionPointTriangle(ray, *triangle, Triangle(), intersect, distanceRay)) {
					//intersectBool = true;
					if (distanceRay > 0) {
						minDist = distanceRay;
						counter = counter + 1;
						//std::cout << "[IsInShadow] : counter: " << counter << std::endl;
						shadow = true;
						goto nextsource;
					}
				}
			}
		}

	nextsource:
		float random;

		/*if (hitSomething) {
			return hitSomething;
		}*/
	}
	shadowCounter = counter;
	return shadow;
}

// returns whether the ray hit something or not
bool Intersect(unsigned int level, const Ray ray, Intersection& intersect, Triangle ignoreTriangle) {
	// if (level > maxRecursionLevel) return false;

	intersect.distance = INFINITY;

	// trace the ray through all the triangles, optimization structure needs to go here still
	for (int i = 0; i < MyMesh.triangles.size(); i++) {
		Vec3Df intersectionPoint;
		float distance;
		Triangle triangle = MyMesh.triangles[i];
		// if an intersection gets found, put the resulting point and triangle in the result vars
		if (rayIntersectionPointTriangle(ray, triangle, ignoreTriangle, intersectionPoint, distance)) {
			// check if distance is smaller than previous result and larger than zero
			if (distance < intersect.distance && distance > 0) {
				intersect.point = intersectionPoint;
				intersect.triangle = triangle;
				intersect.distance = distance;
				intersect.material = MyMesh.materials[MyMesh.triangleMaterials[i]];

				Vec3Df n = calculateSurfaceNormal(triangle);
				n.normalize();

				intersect.schlickCosTheta = fabs(Vec3Df::dotProduct(n, ray.direction));
			}
		}
	}

	if (intersect.distance < 1000000) return true;
	else return false;
}

void Shade(unsigned int level, Ray origRay, Intersection intersect, Vec3Df& color) {
	Vec3Df directColor, reflectedColor, refractedColor, specularLuminance;
	directColor = Vec3Df(0, 0, 0);
	reflectedColor = Vec3Df(0, 0, 0);
	refractedColor = Vec3Df(0, 0, 0);

	Ray reflectedRay, refractedRay;

	bool computeDirect, computeReflect, computeRefract, computeSpecular;
	computeDirect = computeReflect = computeRefract = computeSpecular = false;
	Vec3Df mirrorReflectance = Vec3Df(0, 0, 0);
	Vec3Df glassRefractance = Vec3Df(0, 0, 0);
	Vec3Df diffuseContribution = Vec3Df(1, 1, 1);

	// determine which contributions need to be computed for which materials
	switch (intersect.material.illum()) {
		// color on, ambient on (basically only diffuse)
	case 1:
		computeDirect = true;
		break;
		// highlights on (diffuse + specular highlights)
	case 2:
		computeDirect = true;
		computeSpecular = true;
		break;
		// pure mirror
	case 3:
		computeDirect = true;
		computeReflect = true;
		computeSpecular = true;
		// TODO: tune values
		mirrorReflectance = intersect.material.Ka();
		diffuseContribution = Vec3Df(1, 1, 1) - mirrorReflectance;
		break;
		// // realistic glass (both reflection and refraction)
		// case 4:
		// 	computeReflect = true;
		// 	computeRefract = true;
		// 	computeSpecular = true;
		// 	break;
		// general glossy material (reflection/refraction, ray trace on, fresnel off)
		// TODO: look at partially opaque materials (frosted glass etc)
	case 6:
	{
		computeReflect = true;
		computeRefract = true;
		computeSpecular = true;
		diffuseContribution = Vec3Df(0, 0, 0);
		mirrorReflectance = intersect.material.Ka();

		float R0 = (intersect.material.Ni() - 1.0f) / (intersect.material.Ni() + 1.0f);
		R0 *= R0;
		float sc = (1 - intersect.schlickCosTheta);
		float schlickReflectance = R0 + (1 - R0)*sc*sc*sc*sc*sc;

		glassRefractance = Vec3Df(1, 1, 1) - mirrorReflectance;

		mirrorReflectance += glassRefractance * schlickReflectance;
		glassRefractance *= (1 - schlickReflectance);
		break;
	}
	// pure refractive material (no reflections)
	case 9:
		computeRefract = true;
		computeSpecular = true;
		diffuseContribution = Vec3Df(0, 0, 0);
		break;
		// by default compute nothing and give a warning in the terminal
	default:
		printf("Warning: unknown material type %d, please check...", intersect.material.illum());
		break;
	}

	if (computeDirect) ComputeDirectLight(intersect, directColor);

	if ((computeReflect || computeSpecular) && level + 1 <= maxRecursionLevel) {
		if (ComputeReflectedRay(origRay, intersect.point, intersect.triangle, reflectedRay)) {
			if (computeReflect) Trace(level + 1, reflectedRay, reflectedColor, intersect.triangle);
			if (computeSpecular) BounceLight(reflectedRay, specularLuminance, intersect.triangle);
		}
	}

	if (computeRefract && level + 1 <= maxRecursionLevel) {
		if (ComputeRefractedRay(origRay, intersect, refractedRay)) Trace(level + 1, refractedRay, refractedColor, intersect.triangle);
	}

	// TODO: figure out proper mirrorReflectance/specular usage, both should be handled differently
	color = diffuseContribution * intersect.material.Kd()*directColor + specularLuminance * intersect.material.Ks() + mirrorReflectance * reflectedColor + glassRefractance * refractedColor;

	if (drawRecurseRays && level == 0) std::cout << "Got color " << color << " from level " << level << std::endl;

	return;
}

void BounceLight(Ray ray, Vec3Df& luminance, Triangle ignoreTriangle) {
	// TODO
	luminance = Vec3Df(0, 0, 0);
}

void Trace(unsigned int level, Ray ray, Vec3Df& color, Triangle ignoreTriangle) {
	Intersection intersect;

	if (Intersect(level, ray, intersect, ignoreTriangle)) {
		if (intersect.distance < 0) {
			Vec3Df newRayDir = ray.origin - intersect.point;
			newRayDir.normalize();
			intersect.point = ray.origin + 5 * newRayDir;
		}

		if (drawRecurseRays) {
			recurseTestRayOrigins[recurseTestRayCount] = ray.origin;
			recurseTestRayDestinations[recurseTestRayCount] = intersect.point;
			recurseTestRayCount++;

			std::cout << "  Traced a ray on level " << level << " from " << recurseTestRayOrigins[recurseTestRayCount - 1] << " to " << recurseTestRayDestinations[recurseTestRayCount - 1] << ". Travelled " << intersect.distance << std::endl;
		}

		if (intersect.distance >= 0) Shade(level, ray, intersect, color);
	}
	else {
		color = Vec3Df(0, 0, 0);
	}

	return;
}

void ComputeDirectLight(Intersection intersect, Vec3Df& directColor) {
	// if (isInShadow(pointOfIntersection, Triangle())) {
	directColor = Vec3Df(1, 1, 1); // Hard shadow, item is in the shadow thus color is black
// }

/*
Better if this method is removed or the return should be a boolean, it is not necessary  to calculate colors, reflections if there's a shadow cast on the triangle/intersection.
*/

	return;
}

bool ComputeReflectedRay(Ray origRay, Vec3Df pointOfIntersection, Triangle triangleOfIntersection, Ray& reflectedRay) {
	Vec3Df n = calculateSurfaceNormal(triangleOfIntersection);
	n.normalize();

	reflectedRay.origin = pointOfIntersection;
	// calculate reflected direction vector
	reflectedRay.direction = origRay.direction - 2 * n*(Vec3Df::dotProduct(origRay.direction, n));
	reflectedRay.direction.normalize();

	return true;
}

bool ComputeRefractedRay(Ray origRay, Intersection intersect, Ray& refractedRay) {
	Vec3Df n = calculateSurfaceNormal(intersect.triangle);
	n.normalize();

	// calculate refracted direction vector
	Vec3Df v = origRay.direction;
	float v_dot_n = Vec3Df::dotProduct(v, n);

	// determine if we are at material-air or air-material interface from the dotproduct of the normal and the ray
	// NOTE: this code assumes all transitions to be with air on one side, so leave space between objects!
	float n1, n2;
	if (v_dot_n > 0) {
		n1 = intersect.material.Ni();
		n2 = 1.0f;
		n = -n;
		v_dot_n = Vec3Df::dotProduct(v, n);
	}
	else {
		n1 = 1.0f;
		n2 = intersect.material.Ni();
	}

	refractedRay.origin = intersect.point;

	// if the sqrt is negative, this will produce a NaN value, not cause an error
	refractedRay.direction = n1 / n2 * (v - v_dot_n * n) - n * sqrt(1 - (n1*n1*(1 - v_dot_n * v_dot_n)) / (n2*n2));
	refractedRay.direction.normalize();

	// check if sqrt returned something negative
	if (std::isnan(refractedRay.direction.p[0]) || std::isnan(refractedRay.direction.p[1]) || std::isnan(refractedRay.direction.p[2])) return false;

	return true;
}

/**
 * Moller and Trumbore algorithnm: point(u,v) = (1-u-v)*p0 + u*p1 + v*p2
 * returns true if the ray intersects with the triangle.
*/
bool rayIntersectionPointTriangle(Ray r, const Triangle & triangle, const Triangle & ignoreTriangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection)
{
	if (&ignoreTriangle == &triangle) return false;
	float epsilon = 0.001f;
	if (fabs(triangle.v[0] - ignoreTriangle.v[0]) < epsilon)
		if (fabs(triangle.v[1] - ignoreTriangle.v[1]) < epsilon)
			if (fabs(triangle.v[2] - ignoreTriangle.v[2]) < epsilon)
				return false;


	Vec3Df edges[2];
	float t;
	float u;
	float v;

	edges[0] = MyMesh.vertices[triangle.v[1]].p - MyMesh.vertices[triangle.v[0]].p;
	edges[1] = MyMesh.vertices[triangle.v[2]].p - MyMesh.vertices[triangle.v[0]].p;
	Vec3Df n = Vec3Df::crossProduct(edges[0], edges[1]);
	n.normalize();

	Vec3Df T = r.origin - MyMesh.vertices[triangle.v[0]].p;

	Vec3Df P = Vec3Df::crossProduct(r.direction, edges[1]);
	Vec3Df Q = Vec3Df::crossProduct(T, edges[0]);

	float invDet = 1.0 / Vec3Df::dotProduct(P, edges[0]); // 1 / P . E1

	t = invDet * Vec3Df::dotProduct(Q, edges[1]);
	u = invDet * Vec3Df::dotProduct(P, T);
	v = invDet * Vec3Df::dotProduct(Q, r.direction);

	// check whether there is an intersection
	if (u < 0.0 || v < 0.0 || u + v > 1.0) return false; // u >= 0, v >= 0, u + v <= 1

														 //	pointOfIntersection = MyMesh.vertices[triangle.v[0]].p + u * edges[0] + v * edges[1];
	pointOfIntersection = (1 - u - v) * MyMesh.vertices[triangle.v[0]].p + u * MyMesh.vertices[triangle.v[1]].p + v * MyMesh.vertices[triangle.v[2]].p;
	distanceLightToIntersection = t;

	// std::cout << "Found intersecting triangle (" << pointOfIntersection << ") at distance " << t << std::endl;

	return true;
}

/**
* pin - point of intersection going *in* the box
* pout - point of intersection going *out* the box
* returns true if the ray intersects with the box
* Reference: https://www.youtube.com/watch?v=USjbg5QXk3g (Math for Game Developers - Bullet Collision (Vector/AABB Intersection))
*/
bool rayIntersectionPointBox(Ray r, Box box, Vec3Df& pin, Vec3Df& pout)
{
	Vec3Df min = box.min;
	Vec3Df max = box.max;

	float tminX, tminY, tminZ, tmaxX, tmaxY, tmaxZ, tinX, tinY, tinZ, toutX, toutY, toutZ, tin, tout;

	tminX = (min.p[0] - r.origin.p[0]) / r.direction.p[0]; // (x - Ox) / Dx
	tminY = (min.p[1] - r.origin.p[1]) / r.direction.p[1]; // (y - Oy) / Dy
	tminZ = (min.p[2] - r.origin.p[2]) / r.direction.p[2]; // (z - Oz) / Dz

	tmaxX = (max.p[0] - r.origin.p[0]) / r.direction.p[0]; // (x - Ox) / Dx
	tmaxY = (max.p[1] - r.origin.p[1]) / r.direction.p[1]; // (y - Oy) / Dy
	tmaxZ = (max.p[2] - r.origin.p[2]) / r.direction.p[2]; // (z - Oz) / Dz

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

	pin = r.origin + r.direction * tin;
	pout = r.origin + r.direction * tout;
	return true;
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

void drawRay(Vec3Df rayOrig, Vec3Df rayDest, Vec3Df color) {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(0, 1, 1);
	glVertex3f(rayOrig[0], rayOrig[1], rayOrig[2]);
	glColor3f(color[0], color[1], color[2]);
	glVertex3f(rayDest[0], rayDest[1], rayDest[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();
}

// calculates the centroid of the triangle
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
	{
		glVertex3fv(MyLightPositions[i].pointer());
	}
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state.
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous
	//state after the pop.


	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	// glBegin(GL_LINES);
	// glColor3f(0, 1, 1);
	// glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	// glColor3f(1, 0, 0);
	// glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	// glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();

	// draw the bounding boxes
	//for (Box box : root->getBoxes())
	//{
	//	box.highlightEdges();
	//}

	for (int i = 0; i < 20; i++) {
		drawRay(recurseTestRayOrigins[i], recurseTestRayDestinations[i], Vec3Df(0, 1, 0));
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

	// We use a seed so that every scene will be the same for all lights,
	// even though we are using random points.
	std::mt19937 seed(light_speed_sphere);
	std::uniform_real_distribution<double> rndFloat(0.0, 1.0);

	// Retrieve the values of the current light.
	Vec3Df lightPosition = MyLightPositions[MyLightPositions.size() - 1];
	float lightSphereWidth = MyLightPositionRadius[MyLightPositionRadius.size() - 1];
	int lightSphereAmount = MyLightPositionAmount[MyLightPositionAmount.size() - 1];

	// Create the list of points.
	std::vector<Vec3Df> currentLightSphere;

	// We only calculate the lightSphere if it is actually needed, else we just use the MyLightPositions.
	if (MyLightPositionAmount[MyLightPositionAmount.size() - 1] > 1 && MyLightPositionRadius[MyLightPositionRadius.size() - 1] > 0) {

		// Calculate position for every surface light.
		for (int i = 0; i < lightSphereAmount; i++) {
			double theta = 2 * M_PI * rndFloat(seed);
			double phi = acos(1 - 2 * rndFloat(seed));
			double x = lightPosition[0] + sin(phi) * cos(theta) * lightSphereWidth;
			double y = lightPosition[1] + sin(phi) * sin(theta) * lightSphereWidth;
			double z = lightPosition[2] + cos(phi) * lightSphereWidth;
			Vec3Df offset = Vec3Df(x, y, z);
			MyLightPositions.push_back(offset);
		}
	}
	else {
		// We just add the normal light position.
		MyLightPositions.push_back(lightPosition);
	}
	// We add list of points around the sphere into the list.
	MySphereLightPositions.push_back(currentLightSphere);
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
	}
	else {
		return minimum;
	}
}

/**
 * Method to calculate the diffuse light.
 *
 * @param vertexPos the selected position.
 * @param normal The normal.
 * @param material The material.
 * @param lightPos Our lightposition.
 * @return The diffuse light.
 */
Vec3Df diffuseFunction(const Vec3Df &vertexPos, Vec3Df &normal, Material *material, Vec3Df lightPos) {

	Vec3Df vectorLight = (lightPos - vertexPos);
	vectorLight.normalize();
	return material->Kd() * std::max(0.0f, Vec3Df::dotProduct(normal, vectorLight));

}

/**
 * Method to calculate the specular light.
 *
 * @param vertexPosition the selected position.
 * @param normal the normal.
 * @param material the material.
 * @param lightPosition the lightposition.
 * @return The specular light.
 */
Vec3Df specularFunction(const Vec3Df &vertexPosition, Vec3Df &normal, Material *material, Vec3Df lightPosition) {

	Vec3Df vectorView = vertexPosition - MyCameraPosition;
	vectorView.normalize();

	Vec3Df vectorLight = lightPosition - vertexPosition;
	vectorLight.normalize();

	normal.normalize();
	Vec3Df reflection = vectorLight - 2 * (Vec3Df::dotProduct(normal, vectorLight)) * normal;
	reflection.normalize();

	float dotProduct = std::max(Vec3Df::dotProduct(vectorView, reflection), 0.0f);

	// take the power
	return material->Ks() * pow(dotProduct, material->Ns());

}

/**
 * Method to calculate the shading.
 *
 * @param ray the ray.
 * @param vertexPos  the vertex position.
 * @param normal the normal.
 * @param material the material.
 * @return The shading on impact.
 */
Vec3Df calculateShading(const Vec3Df &vertexPosition, Vec3Df &normal, Material *material) {

	// initially black
	Vec3Df calculatedColor(0, 0, 0);

	// Add the ambient shading
	if (ambient && material->has_Ka()) {
		calculatedColor = calculatedColor + material->Ka();
	}

	// Calculate the shading for every light source.
	for (int i = 0; i < MyLightPositions.size(); i++) {

		// Add the diffuse shading.
		if (diffuse && material->has_Kd()) {
			calculatedColor = calculatedColor + material->Tr() * diffuseFunction(vertexPosition, normal, material, MyLightPositions[i]);
		}

		// Add the specular shading.
		if (specular && material->has_Ks() && material->has_Ns()) {
			calculatedColor = calculatedColor + material->Tr() * specularFunction(vertexPosition, normal, material, MyLightPositions[i]);
		}
	}

	return calculatedColor;
}

std::vector<Box>* getIntersectionBoxes(Ray r, Tree<Box>* curr)
{
	Vec3Df pin, pout;
	if (curr == NULL) return &std::vector<Box>();
	if (!rayIntersectionPointBox(r, curr->data, pin, pout)) return &std::vector<Box>();
	std::vector<Box>* boxes = getIntersectionBoxes(r, curr->left);
	std::vector<Box>* rightBoxes = getIntersectionBoxes(r, curr->right);
	boxes->insert(boxes->end(), rightBoxes->begin(), rightBoxes->end());
	boxes->push_back(curr->data);
	return boxes;
}

std::vector<Box>* getIntersectionLeaves(Ray r, Tree<Box>* curr)
{
	Vec3Df pin, pout;
	if (curr == NULL) return &std::vector<Box>();
	if (!rayIntersectionPointBox(r, curr->data, pin, pout)) return &std::vector<Box>();
	std::vector<Box>* boxes = getIntersectionBoxes(r, curr->left);
	std::vector<Box>* rightBoxes = getIntersectionBoxes(r, curr->right);
	boxes->insert(boxes->end(), rightBoxes->begin(), rightBoxes->end());
	if (curr->left == NULL && curr->right == NULL) boxes->push_back(curr->data);
	return boxes;
}

/**
 * Gets the first box int the tree "curr", that intersects with the ray.
 * This Function should only be called within these pre-conditions:
 * 1 - BoxTree *curr is not NULL
 * 2 - The ray, atleast, intersects the overlapping boundingbox, or in other words the initial "curr->data".
**/
Box getFirstIntersectedBox(Ray r, Tree<Box>* curr, Vec3Df& pin, Vec3Df& pout)
{
	if (curr->left == NULL && curr->right == NULL) // if we hit a leaf
	{
		return curr->data;
	}
	else if (curr->left == NULL)
	{
		if (rayIntersectionPointBox(r, curr->right->data, pin, pout))
		{
			return getFirstIntersectedBox(r, curr->right, pin, pout);
		}
		return curr->data;
	}
	else if (curr->right == NULL)
	{
		if (rayIntersectionPointBox(r, curr->left->data, pin, pout))
		{
			return getFirstIntersectedBox(r, curr->left, pin, pout);
		}
		return curr->data;
	}
	else if (rayIntersectionPointBox(r, curr->left->data, pin, pout))
	{
		Vec3Df testpin, pout;
		if (rayIntersectionPointBox(r, curr->right->data, testpin, pout) && Vec3Df::distance(testpin, r.origin) < Vec3Df::distance(pin, r.origin))
		{
			return getFirstIntersectedBox(r, curr->right, pin, pout);
		}
		return getFirstIntersectedBox(r, curr->left, pin, pout);
	}
	else if (rayIntersectionPointBox(r, curr->right->data, pin, pout))
	{
		return getFirstIntersectedBox(r, curr->right, pin, pout);
	}
	return curr->data;
}

const Tree<Box>* getFirstIntersectedBoxFast(Ray r, const Tree<Box>* curr, Vec3Df& pin, Vec3Df& pout)
{
	Vec3Df testpinL, testpinR;
	if (curr->left == NULL || curr->right == NULL)
	{
		return curr;
	}

	bool intersectL = rayIntersectionPointBox(r, curr->left->data, testpinL, pout);
	bool intersectR = rayIntersectionPointBox(r, curr->right->data, testpinR, pout);

	if (intersectL && intersectR)
	{
		if (Vec3Df::distance(testpinL, r.origin) < Vec3Df::distance(testpinR, r.origin))
			return getFirstIntersectedBoxFast(r, curr->left, pin, pout);
		return getFirstIntersectedBoxFast(r, curr->right, pin, pout);
	}
	else if (intersectL)
	{
		return getFirstIntersectedBoxFast(r, curr->left, pin, pout);
	}
	else if (intersectR)
	{
		return getFirstIntersectedBoxFast(r, curr->right, pin, pout);
	}

	return curr;
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
			const Triangle* triangle = &MyMesh.triangles[i];

			Vec3Df pointOfIntersection;
			float distanceRay;
			Vec3Df color;
			Vec3Df hit;
			//				if (rayIntersectionPointPlane(rayOrigin, normRayDirection, triangle, pointOfIntersection))
			//				{
			//					std::cout << pointOfIntersection << std::endl;
			//				}
			Ray r = { rayOrigin, normRayDirection };

			if (rayIntersectionPointTriangle(r, *triangle, Triangle(), pointOfIntersection, distanceRay))
			{
				std::cout << "Ray InterSects Triangle: " << pointOfIntersection << std::endl;
				//					Material mat = MyMesh.materials[MyMesh.triangleMaterials[i]];
				//					mat.set_Kd(0, 0, 0);
				//					std::cout << (mat.Kd()) << std::endl;
			}
		}
	}
	break;
	case 's':
	{
		Triangle t;
		Vec3Df intersection = DebugRay(testRayOrigin, testRayDestination, t, root);
		testRayDestination = intersection;
		/*std::cout << "Intersection Point: " << testRayDestination <<std::endl;
		std::cout << "Intersection Point: " << intersection <<std::endl;*/

		if (isInShadow(intersection, t, root)) {
			std::cout << "Intersection lies in the : SHADOW" << std::endl;

		}
		else {
			std::cout << "Intersection lies in the : LIGHT" << std::endl;
		}

	}
	break;
	case 'd': // constructs axis aligned bounding boxes
	{
		root->print(0);
		std::cout << root->size() << std::endl;
		root->highlightEdges();
	}
	break;
	case 'f':
	{

		Ray r;
		r.origin = rayOrigin;
		r.direction = rayDestination - rayOrigin;

		// Drawing the boxes
		//showBoxes(&tree);
		//showLeavesOnly(&tree);
		//showIntersectionBoxOnly(rayOrigin, normRayDirection, &tree);
		//showIntersectionLeafOnly(r, &tree);
		// getAllIntersectedLeafs(r, &tree, Vec3Df(), Vec3Df(), boxes);
	}
	break;
	case 'x':
		createLightPointer();
		break;
	case 'w':
		MyLightPositions[MyLightPositionsPointer] = MyCameraPosition;
		setupMySphereLightPositions();
		break;
	case 'l':
		if (MyLightPositionsPointer < MyLightPositions.size() - 1) {
			MyLightPositionsPointer++;
		}
		else {
			MyLightPositionsPointer = 0;
		}
		break;
	case 'p':
		MyLightPositionPower[MyLightPositionsPointer] += 50;
		break;
	case 'P':
		MyLightPositionPower[MyLightPositionsPointer] -= 50;
		if (MyLightPositionPower[MyLightPositionsPointer] < 0) {
			MyLightPositionPower[MyLightPositionsPointer] = 0;
		}
		break;
	case 'o':
		MyLightPositionRadius[MyLightPositionsPointer] += .01f;
		setupMySphereLightPositions();
		break;
	case 'O':
		MyLightPositionRadius[MyLightPositionsPointer] -= .01f;
		if (MyLightPositionRadius[MyLightPositionsPointer] < 0) {
			MyLightPositionRadius[MyLightPositionsPointer] = 0;
		}
		else {
			setupMySphereLightPositions();
		}
		break;
	case 'a':
		MyLightPositionAmount[MyLightPositionsPointer] += 1;
		setupMySphereLightPositions();
		break;
	case 'A':
		MyLightPositionAmount[MyLightPositionsPointer] -= 1;
		if (MyLightPositionAmount[MyLightPositionsPointer] < 1) {
			MyLightPositionAmount[MyLightPositionsPointer] = 1;
		}
		else {
			setupMySphereLightPositions();
		}
		break;
	case 'q':
	{
		drawRecurseRays = true;

		resetRecurseTestRays();

		Ray recurseTestRay;

		Vec3Df resultColor;

		recurseTestRay.origin = rayOrigin;
		recurseTestRay.direction = normRayDirection;

		std::cout << "Starting recursive raytrace..." << std::endl;

		Trace(0, recurseTestRay, resultColor, Triangle());

		std::cout << "Recursive ray trace done." << std::endl;

		drawRecurseRays = false;
	}
	break;
	default:
		break;
	}

	std::cout << t << " pressed! The mouse was in location " << x << "," << y << "!" << std::endl;
}