#ifndef SWARM_H
#define SWARM_H

#include <stdarg.h>
#include "pdecomp.h"

#define MAX_INPUTS 3
#define MAX_PARTICLES 25
#define V_MAX 10

#define MAX_ITER    100

//Main PSO algorithm
int part_swarm(int p_index);

//Calculates particle velocity
void getVelocity(int gBestIndex);

//Updates particle position according to getVelocity

void  updateParticles(int gBestIndex);

//Initialises particle position - random
float  initialise();

//
float testProblem(float temp_pos[], float trace[]);
float testProblem2(int p_index, int index);
//
float gRand();

//
float   getRandNum(int low, int high);

//
int   minimum();

int setRange(int test_seg);

double r_norm(double mu, double sigma);
#endif

