#include "actu.h"
#include <math.h>
void setDensity(point* points, neighborhood* nh, double kh)
{
	for (int i = 0; i < NPTS; i++)
	{
		double density = MASS * w_lucy(0, kh);
		neighbours* point = nh[i].list;
		for (int j = 0; j < (int)nh[i].nNeighbours; j++)
		{
			density += MASS * w_lucy(point->distance, kh);
			point = point->next;
		}
		points[i].density = density;
	}
}
void updateDensity(point* points, neighborhood* nh, double kh)
{
	float New_rho[2400];
	double dt = 0.1;
	for (int i = 0; i < NPTS; i++)
	{
		double der_density = 0;
		neighbours* point = nh[i].list;
		double speed_i = sqrt(points[i].vx * points[i].vx + points[i].vy * points[i].vy);
		for (int j = 0; j < (int)nh[i].nNeighbours; j++)
		{
			double speed_j = sqrt(points[point->index].vx * points[point->index].vx + points[point->index].vy * points[point->index].vy);
			der_density += MASS * (speed_i - speed_j) * grad_w_lucy(point->distance, kh);
			point = point->next;
		}
		New_rho[i] = dt * der_density;//Euler pour le moment
		// Maintenant qu'on a les Drho/dt pour chaque points, on peut faire évoluer la valeur vie leapfrog ou autre.
		// Cette valeur peut ensuite être corrigée avec la méthode CSPM
	}
	// METHODE CSPM
	for (int i = 0; i < NPTS; i++)
	{
		double num = 1;
		double denom = 1 / New_rho[i]; 
		neighbours* point = nh[i].list;
		for (int j = 0; j < (int)nh[i].nNeighbours; j++)
		{
			num += w_lucy(point->distance, kh);
			denom += w_lucy(point->distance, kh) / New_rho[point->index];
			point = point->next;
		}
		points[i].density = num / denom;
	}
}

void updatePressure(point* points)
{
	double rho_0 = 1000; // kg/m^3
	double B = 101300; // Pa
	double gamma = 7;
	for (int i = 0; i < NPTS; i++)
	{
		points[i].pressure = B * (pow(points[i].density / rho_0, gamma) - 1);
	}
}

void updateVelocity(point* points, neighborhood* nh, double kh)
{
	float vx[2400];
	float vy[2400];
	double kin_v = 0.000001;// kinematic viscosity m^2/s
	double g = 9.81; // gravity m/s^2
	for (int i = 0; i < NPTS; i++)
	{
		double dvx = 0;
		double dvy = g;
		neighbours* point = nh[i].list;
		for (int j = 0; j < (int)nh[i].nNeighbours; j++)
		{
			int l = point->index;
			dvx += MASS * (points[i].pressure / (points[i].density * points[i].density) + points[l].pressure / (points[l].density * points[l].density)) * grad_w_lucy1D(point->distX, kh) + 2 * kin_v * MASS / (points[l].density) * (points[i].vx - points[l].vx) * grad_w_lucy1D(point->distX, kh) / point->distX;
			dvy += MASS * (points[i].pressure / (points[i].density * points[i].density) + points[l].pressure / (points[l].density * points[l].density)) * grad_w_lucy1D(point->distY, kh) + 2 * kin_v * MASS / (points[l].density) * (points[i].vy - points[l].vy) * grad_w_lucy1D(point->distY, kh) / point->distY;
			point = point->next;
		}
		vx[i] = points[i].vx + dvx * dt; // euler method
		vy[i] = points[i].vy + dvy * dt;
	}
	for (int i = 0; i < NPTS; i++)
	{
		points[i].vx = vx[i];
		points[i].vy = vy[i];
	}
}

void updatePosition(point* points)
{
	for (int i = 0; i < NPTS; i++)
	{
		points[i].x += points[i].vx * dt; //euler
		points[i].y += points[i].vy * dt;
	}
}
