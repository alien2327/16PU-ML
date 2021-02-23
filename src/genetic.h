#ifndef GENETIC
#define GENETIC

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <thread>
#include <vector>
#include <complex>
#include <time.h>
#include <cmath>

/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

using namespace Eigen;
typedef Matrix<double, 2, 1000> Beam_t;
typedef Matrix<double, 2, 8> tBeam_t;
typedef Matrix<double, 16, 1> Vol_t;

static unsigned long mt[MT_N]; /* the array for the state vector  */
static int mti=MT_N+1; /* mti==MT_N+1 means mt[MT_N] is not initialized */

void oper(Beam_t &, int, int, double);
void applyOper(Beam_t &, int, double);
void mutation(Beam_t &, int, double);
void applyMut(Beam_t &, double);
void initPos(Beam_t &);
void testPart(tBeam_t &);
MatrixXd func(MatrixXd, MatrixXd);

double mse(MatrixXd, MatrixXd);
double ep(MatrixXd, MatrixXd);
double lossmin(MatrixXd, MatrixXd::Index &, MatrixXd::Index &);

double mean(MatrixXd, int);
double meansqure(MatrixXd, int);
double variance(MatrixXd, int);
double standard_deviation(MatrixXd, int);

double Uniform(int);
double rand_normal(double, double);

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

#endif