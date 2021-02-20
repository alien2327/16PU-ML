#ifndef GENETIC
#define GENETIC

#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <ctime>
#include <complex>
#include <cmath>

#define N_Particle 200

using namespace std;
double mrate = 0.005;

double Uniform(int);
void oper(double *, double, int);
double rand_normal(double , double);
void initPos(vector<vector<double>> &, double, double);
void func(vector<double> &, double [16]);
void mutation(vector<vector<double>> &);

// Evalution function
double lossfunction(double *, double *);
double mse(double *, double *);
double ep(double *, double *);

#endif