#ifndef ACTU_H
#define ACTU_H

#include "BOV.h"
#include "./neighborhood_search.h"
#include "./point.h"
#include <time.h>
#include <math.h>

extern int NPTS;
extern double MASS;
extern double SOURCE_TEMP;
extern double INIT_TEMP;
extern double dt;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// Structure to represent a point

void updateDensity(point* points, neighborhood* nh, double kh);

void setDensity(point* points, neighborhood* nh, double kh);

void updatePressure(point* points);

void updateVelocity(point* points, neighborhood* nh, double kh);

void updatePosition(point* points);

#endif
