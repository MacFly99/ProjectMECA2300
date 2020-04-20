#ifndef POINT_H
#define POINT_H

#include "BOV.h"
#include <time.h>
#include <math.h>

extern int NPTS;
extern double SOURCE_TEMP;
extern double INIT_TEMP;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// Structure to represent a point
typedef struct point {
	double x;
	double y;
	double vx;
	double vy;
	double val;
	double div;
	double grad_x;
	double grad_y;
	double lapl;
}point;

void fillPointsRand(point* points);

void fillPointsGrid(point* points);

void updateData(point* points, GLfloat(*data)[8]);

static void colormap(float v, float color[3]);

#endif