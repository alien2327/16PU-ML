#ifndef GENETIC
#define GENETIC

#include <iostream>
#include <Eigen/Dense>
#include <thread>
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
typedef Matrix<double, 2, 100> Beam_t;

static unsigned long mt[MT_N]; /* the array for the state vector  */
static int mti=MT_N+1; /* mti==MT_N+1 means mt[MT_N] is not initialized */

void oper(Beam_t &, double, int);
void initPos(Beam_t &);
void mutation();

double lossfunction(MatrixXd, MatrixXd);
double mse(MatrixXd, MatrixXd);
double ep(MatrixXd, MatrixXd);

double mean(MatrixXd, int);
double variance(MatrixXd, int);
double standard_deviation(MatrixXd, int);

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

#endif