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

Vec3Df recurseTestRayOrigins[20];
Vec3Df recurseTestRayDestinations[20];

unsigned int recurseTestRayCount;
// define function before implementation
BoxTree initBoxTree();
void initAccelerationStructure();

std::vector<AABB> boxes;
BoxTree tree = BoxTree(AABB());

// std::map<std::string, unsigned int> materialIndex;

unsigned int maxRecursionLevel;

// Ray structure
struct Ray {
	Vec3Df origin;
	Vec3Df direction;
	bool insideMaterial;
};

struct Intersection {
	Vec3Df point;
	Triangle triangle;
	Material material;
	float distance;
};

// Some forward declarations needed because of recursive functions
void ComputeDirectLight(Vec3Df pointOfIntersection, Vec3Df& directColor);
bool ComputeReflectedRay(Ray origRay, Vec3Df pointOfIntersection, Triangle triangleOfIntersection, Ray& reflectedRay);
bool ComputeRefractedRay(Ray origRay, Intersection intersect, Ray& refractedRay);
void Trace(unsigned int level, Ray ray, Vec3Df& color, Triangle ignoreTriangle);

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
	// MyMesh.loadMtl("test_scene_1.mtl", materialIndex);
	MyMesh.loadMesh("test_scene_1.obj", true);
	MyMesh.computeVertexNormals();

	tree = initBoxTree();
	initAccelerationStructure();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
    //MyLightPositions.push_back(MyCameraPosition);
    createLightPointer();

	maxRecursionLevel = 2;
	recurseTestRayCount = 0;

	for (int i=0; i < 20; i++) {
		recurseTestRayOrigins[i] = Vec3Df(0,0,0);
		recurseTestRayDestinations[i] = Vec3Df(0,0,0);
	}
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
			for (int i = 0; i < box.triangles.size(); i++)
			{
				Vec3Df pointOfIntersection;
				float distanceRay;
				Triangle triangle = box.triangles[i];
				if (rayIntersectionPointTriangle(origin, direction, triangle, Triangle(), pointOfIntersection, distanceRay))
				{
					return Vec3Df(dest[0], dest[1], dest[2]);
				}
			}
		}
	}
	return Vec3Df(0, 0, 0);
	//return Vec3Df(dest[0], dest[1], dest[2]);
}

// returns whether the ray hit something or not
bool Intersect(unsigned int level, const Ray ray, Intersection& intersect, Triangle ignoreTriangle) {
	if (level > maxRecursionLevel) return false;

	intersect.distance = INFINITY;

	// trace the ray through all the triangles, optimization structure needs to go here still
	for (int i = 0; i < MyMesh.triangles.size(); i++) {
		Vec3Df intersectionPoint;
		float distance;
		Triangle triangle = MyMesh.triangles[i];
		// if an intersection gets found, put the resulting point and triangle in the result vars
		if (rayIntersectionPointTriangle(ray.origin, ray.direction, triangle, ignoreTriangle, intersectionPoint, distance)) {
			// check if distance is smaller than previous result and larger than zero
			if (distance < intersect.distance && distance > 0) {
				intersect.point = intersectionPoint;
				intersect.triangle = triangle;
				intersect.distance = distance;
				intersect.material = MyMesh.materials[MyMesh.triangleMaterials[i]];
			}
			// if (distanceRay < 0) pointOfIntersection = pointOfIntersection + 5*ray.direction;
		}
	}

	if (intersect.distance < 1000000) return true;
	else return false;
}

void Shade(unsigned int level, Ray origRay, Intersection intersect, Vec3Df& color) {
	Vec3Df directColor, reflectedColor, refractedColor;
	directColor = Vec3Df(0,0,0);
	reflectedColor = Vec3Df(0,0,0);
	refractedColor = Vec3Df(0,0,0);

	Ray reflectedRay, refractedRay;
	float reflection, transmission;

	// temporary reflection and transmission values of 0.5
	// TODO: calculate correct values according to material and angle
	reflection = 0.5;
	transmission = 0.5;

	bool computeDirect, computeReflect, computeRefract;
	computeDirect = computeReflect = computeRefract = false;
	Vec3Df mirrorReflectance = Vec3Df(0,0,0);

	// determine which contributions need to be computed for which materials
	switch (intersect.material.illum()) {
		// color on, ambient on (basically only diffuse)
		case 1:
			computeDirect = true;
			break;
		// highlights on (diffuse + specular highlights)
		case 2:
			computeDirect = true;
			computeReflect = true;
			break;
		// pure mirror
		case 3:
			computeDirect = true;
			computeReflect = true;
			mirrorReflectance = intersect.material.Ka();
			break;
		// realistic glass (both reflection and refraction)
		case 4:
			computeReflect = true;
			computeRefract = true;
			break;
		// general glossy material (reflection/refraction, ray trace on, fresnel off)
		case 6:
			computeReflect = true;
			computeRefract = true;
			break;
		// pure refractive material (no reflections)
		case 9:
			computeRefract = true;
			break;
		// by default enable all three, but give a warning in the terminal
		default:
			computeDirect = true;
			computeReflect = true;
			computeRefract = true;
			printf("Warning: unknown material type %d, please check...", intersect.material.illum());
			break;
	}

	if (computeDirect) ComputeDirectLight(intersect.point, directColor);

	if (computeReflect) {
		if (ComputeReflectedRay(origRay, intersect.point, intersect.triangle, reflectedRay)) Trace(level + 1, reflectedRay, reflectedColor, intersect.triangle);
	}

	if (computeRefract) {
		if (ComputeRefractedRay(origRay, intersect, refractedRay)) Trace(level + 1, refractedRay, refractedColor, intersect.triangle);
	}

	// TODO: figure out proper mirrorReflectance/specular usage, both should be handled differently
	color = (Vec3Df(1,1,1) - mirrorReflectance)*intersect.material.Kd()*directColor + mirrorReflectance*reflection*reflectedColor + intersect.material.Ks()*transmission*refractedColor;

	return;
}

void Trace(unsigned int level, Ray ray, Vec3Df& color, Triangle ignoreTriangle) {
	Intersection intersect;

	if (Intersect(level, ray, intersect, ignoreTriangle)) {
		if (intersect.distance < 0) {
			Vec3Df newRayDir = ray.origin - intersect.point;
			newRayDir.normalize();
			intersect.point = ray.origin + 5*newRayDir;
		}

		recurseTestRayOrigins[recurseTestRayCount] = ray.origin;
		recurseTestRayDestinations[recurseTestRayCount] = intersect.point;
		recurseTestRayCount++;

		std::cout << "  Traced a ray on level " << level << " from " << recurseTestRayOrigins[recurseTestRayCount - 1] << " to " << recurseTestRayDestinations[recurseTestRayCount - 1] << ". Travelled " << intersect.distance << std::endl;

		if (intersect.distance >= 0) Shade(level, ray, intersect, color);
	} else {
		color = Vec3Df(0,0,0);
	}

	return;
}

void ComputeDirectLight(Vec3Df pointOfIntersection, Vec3Df& directColor) {
	directColor = Vec3Df(0,0,0);

	return;
}

bool ComputeReflectedRay(Ray origRay, Vec3Df pointOfIntersection, Triangle triangleOfIntersection, Ray& reflectedRay) {
	Vec3Df n = calculateSurfaceNormal(triangleOfIntersection);
	n.normalize();

	reflectedRay.origin = pointOfIntersection;
	// calculate reflected direction vector
	reflectedRay.direction = origRay.direction - 2*n*(Vec3Df::dotProduct(origRay.direction, n));
	reflectedRay.direction.normalize();
	reflectedRay.insideMaterial = origRay.insideMaterial;

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
		v_dot_n = Vec3Df::dotProduct(v,n);
	} else {
		n1 = 1.0f;
		n2 = intersect.material.Ni();
	}

	refractedRay.origin = intersect.point;

	// if the sqrt is negative, this will produce a NaN value, not cause an error
	refractedRay.direction = n1/n2*(v - v_dot_n*n) - n*sqrt(1 - (n1*n1*(1 - v_dot_n*v_dot_n))/(n2*n2));
	refractedRay.direction.normalize();

	// check if sqrt returned something negative
	if (std::isnan(refractedRay.direction.p[0]) || std::isnan(refractedRay.direction.p[1]) || std::isnan(refractedRay.direction.p[2])) return false;

	refractedRay.insideMaterial = !(origRay.insideMaterial);

	return true;
}

/**
 * Moller and Trumbore algorithnm: point(u,v) = (1-u-v)*p0 + u*p1 + v*p2
 * returns true if the ray intersects with the triangle.
*/
bool rayIntersectionPointTriangle(Vec3Df rayOrigin, Vec3Df rayDirection, Triangle triangle, Triangle ignoreTriangle, Vec3Df& pointOfIntersection, float& distanceLightToIntersection)
{
	if (ignoreTriangle == triangle) return false;

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

	// std::cout << "Found intersecting triangle (" << pointOfIntersection << ") at distance " << t << std::endl;

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
	// glPushAttrib(GL_ALL_ATTRIB_BITS);
	// glDisable(GL_LIGHTING);
	// glBegin(GL_LINES);
	// glColor3f(0, 1, 1);
	// glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	// glColor3f(0, 0, 1);
	// glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	// glEnd();
	// glPointSize(10);
	// glBegin(GL_POINTS);
	// glVertex3fv(MyLightPositions[0].pointer());
	// glEnd();
	// glPopAttrib();

	// draw the bounding boxes
	for (AABB box : boxes)
	{
		box.highlightBoxEdges();
	}

	for (int i=0; i < 20; i++) {
		drawRay(recurseTestRayOrigins[i], recurseTestRayDestinations[i], Vec3Df(0,1,0));
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
	printf("%ld\n", curr->data.triangles.size());
	rec[depth] = 1;
	printTree(curr->right, depth + 1);
	rec[depth] = 0;
	printTree(curr->left, depth + 1);
}

// Recursively adds the nodes of the BoxTree to "boxes" whom elements will be drawn on the screen.
void showBoxes(struct BoxTree* curr)
{
	if (curr == NULL)return;
	boxes.push_back(curr->data);
	showBoxes(curr->left);
	showBoxes(curr->right);
}

void showLeavesOnly(struct BoxTree* curr)
{
	if (curr == NULL)return;
	if (curr->left == NULL && curr->right == NULL) boxes.push_back(curr->data);
	showLeavesOnly(curr->left);
	showLeavesOnly(curr->right);
}

void showIntersectionBoxOnly(Vec3Df rayOrigin, Vec3Df rayDirection, struct BoxTree* curr)
{
	Vec3Df pin, pout;

	if (curr == NULL)return;
	if (rayIntersectionPointBox(rayOrigin, rayDirection, curr->data, pin, pout)) boxes.push_back(curr->data);
	showIntersectionBoxOnly(rayOrigin, rayDirection, curr->left);
	showIntersectionBoxOnly(rayOrigin, rayDirection, curr->right);
}

void showIntersectionLeafOnly(Vec3Df rayOrigin, Vec3Df rayDirection, struct BoxTree* curr)
{
	Vec3Df pin, pout;

	if (curr == NULL)return;
	if (curr->left == NULL && curr->right == NULL && rayIntersectionPointBox(rayOrigin, rayDirection, curr->data, pin, pout)) boxes.push_back(curr->data);
	showIntersectionLeafOnly(rayOrigin, rayDirection, curr->left);
	showIntersectionLeafOnly(rayOrigin, rayDirection, curr->right);
}

BoxTree initBoxTree()
{
	std::pair<Vec3Df, Vec3Df> minMax = getMinAndMaxVertex();
	AABB aabb = AABB(minMax.first, minMax.second);
	return BoxTree(aabb);
}

void initAccelerationStructure()
{
	tree.splitAvg(4000);
	//showBoxes(&tree);
	showBoxes(&tree);
	printTree(&tree, 0);
}

/**
 * Gets the first box int the tree "curr", that intersects with the ray.
 * This Function should only be called within these pre-conditions:
 * 1 - BoxTree *curr is not NULL
 * 2 - The ray, atleast, intersects the overlapping boundingbox, or in other words the initial "curr->data".
**/
AABB getFirstIntersectedBox(Vec3Df rayOrigin, Vec3Df rayDirection, BoxTree* curr, Vec3Df& pin, Vec3Df& pout)
{
	if (curr->left == NULL && curr->right == NULL) // if we hit a leaf
	{
		return curr->data;
	}
	else if (curr->left == NULL)
	{
		return getFirstIntersectedBox(rayOrigin, rayDirection, curr->right, pin, pout);
	}
	else if (curr->right == NULL)
	{
		return getFirstIntersectedBox(rayOrigin, rayDirection, curr->left, pin, pout);
	}
	else if (rayIntersectionPointBox(rayOrigin, rayDirection, curr->left->data, pin, pout))
	{
		Vec3Df testpin, pout;
		if (rayIntersectionPointBox(rayOrigin, rayDirection, curr->right->data, testpin, pout) && Vec3Df::distance(testpin, rayOrigin) < Vec3Df::distance(pin, rayOrigin))
		{
			return getFirstIntersectedBox(rayOrigin, rayDirection, curr->right, pin, pout);
		}
		return getFirstIntersectedBox(rayOrigin, rayDirection, curr->left, pin, pout);
	}
	return getFirstIntersectedBox(rayOrigin, rayDirection, curr->right, pin, pout);
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

				if (rayIntersectionPointTriangle(rayOrigin, normRayDirection, triangle, Triangle(), pointOfIntersection, distanceRay))
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
		boxes.clear();

		// Splitting the tree
		//tree.splitMiddle(4000);
		tree.splitAvg(4000);

		// Drawing the boxes
		showBoxes(&tree);
		//showLeavesOnly(&tree);
		//showIntersectionBoxOnly(rayOrigin, normRayDirection, &tree);
		//showIntersectionLeafOnly(rayOrigin, normRayDirection, &tree);

		// Draw smallest Intersected Box
		//Vec3Df pin, pout;
		//boxes.push_back(getFirstIntersectedBox(rayOrigin, rayDestination - rayOrigin, &tree, pin, pout));

		printTree(&tree, 0);

		/*if (rayIntersectionPointBox(rayOrigin, normRayDirection, boxes[0], pin, pout))
		{
			std::cout << "Ray InterSects Box: \n" << pin << "\n" << pout << std::endl;
		}*/
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
            } else {
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
    case 'r':
            MyLightPositionRadius[MyLightPositionsPointer] += .01f;
            setupMySphereLightPositions();
            break;
    case 'R':
            MyLightPositionRadius[MyLightPositionsPointer] -= .01f;
            if (MyLightPositionRadius[MyLightPositionsPointer] < 0) {
                MyLightPositionRadius[MyLightPositionsPointer] = 0;
            } else {
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
            } else {
                setupMySphereLightPositions();
            }
            break;
	case 'q':
		{
			recurseTestRayCount = 0;

			Ray recurseTestRay;

			Vec3Df resultColor;

			recurseTestRay.origin = rayOrigin;
			recurseTestRay.direction = normRayDirection;
			recurseTestRay.insideMaterial = false;

			std::cout << "Starting recursive raytrace..." << std::endl;

			Trace(0, recurseTestRay, resultColor, Triangle());

			std::cout << "Recursive ray trace done." << std::endl;
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
