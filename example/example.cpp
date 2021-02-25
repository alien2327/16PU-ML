#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <time.h>
#include <iostream>
#include <cmath>
#include <thread>

using namespace std;
using namespace Eigen;

/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

typedef Matrix<double, 2, 1000> Beam_t;

Matrix<double, 16, 16> mat;
static unsigned long mt[MT_N]; /* the array for the state vector  */
static int mti=MT_N+1; /* mti==MT_N+1 means mt[MT_N] is not initialized */

void genSet(Beam_t &, double, double, double, double);
double Uniform(int);
double rand_normal(double, double);
double mean(MatrixXd, int);
double meansqure(MatrixXd, int);
double variance(MatrixXd, int);
double standard_deviation(MatrixXd, int);
void mutation(Beam_t &, int, double, double, double, double);
void applyMut(Beam_t &, double, double, double, double);
MatrixXd func(MatrixXd);
double mse(MatrixXd, MatrixXd);
void applyOper(Beam_t &, int, double);
double ep(MatrixXd, MatrixXd);
void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
void printOut(Beam_t, char *);

int main() {
    Uniform((unsigned)time(NULL)); init_genrand((unsigned)time(NULL));

    char trueName[30] = "../result/true.dat";
    char testName[30] = "../result/test.dat";

    Beam_t trueSet;
    Beam_t testSet[50];
    Beam_t testSort[50];

    MatrixXd trueVol;
    MatrixXd testVol[50];
    std::vector<std::pair<double, int>> loss;

    Matrix<double, 2, 1> trueMeanSquare;
    Matrix<double, 2, 1> testMeanSquare[50];

    double posMax = 30, posMin = -30;
    double sigMax = 15, sigMin = 5;

    double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
    double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();

    genSet(trueSet, mx, my, sx, sy);

    posMax = 60, posMin = -60;
    sigMax = 20, sigMin = 1;

    printOut(trueSet, trueName);

    trueVol = func(trueSet);
    trueMeanSquare(0,0) = standard_deviation(trueSet, 0);
    trueMeanSquare(1,0) = standard_deviation(trueSet, 1);

    printf("\n*\tBeam Position <x>      : %1.5f\n", mean(trueSet, 0));
    printf("*\tBeam Position <y>      : %1.5f\n", mean(trueSet, 1));
    printf("*\tBeam Standard deviation: %1.5f\n", standard_deviation(trueSet, 0));
    printf("*\tBeam Standard deviation: %1.5f\n*\n", standard_deviation(trueSet, 1));

    for (int i = 0; i < 50; i++) {
        double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
        double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();
        genSet(testSet[i], mx, my, sx, sy);
        testVol[i] = func(testSet[i]);
        testMeanSquare[i](0,0) = standard_deviation(testSet[i], 0);
        testMeanSquare[i](1,0) = standard_deviation(testSet[i], 1);
        loss.push_back(std::make_pair(mse(trueVol, testVol[i]) + ep(trueMeanSquare, testMeanSquare[i]), i));
        testSort[i] = testSet[i];
    }
    sort(loss.begin(), loss.end(), std::less<>());
    for (int i = 0; i < 50; i++) {
        testSet[i] = testSort[loss[(int)i/10].second];
    }

    int t = 0;
    //以前の損失より確実に改善された保証がしたい
    //それと一緒に、local minimaから脱出もしたい

    do{
        for (int i = 0; i < 50; i++) {
            applyOper(testSet[i], 9*genrand_real1(), 0.001 + (0.1 - 0.001) * genrand_real1());
            testVol[i] = func(testSet[i]);
            testMeanSquare[i](0,0) = standard_deviation(testSet[i], 0);
            testMeanSquare[i](1,0) = standard_deviation(testSet[i], 1);
            loss.push_back(std::make_pair(mse(trueVol, testVol[i]) + ep(trueMeanSquare, testMeanSquare[i]), i));
            testSort[i] = testSet[i];
        }
        sort(loss.begin(), loss.end(), std::less<>());
        printf("\r*\t%d-st\t%d-st\t%d-st\t%d-st\t%d-st\tbunch selected at\t%d\ttry", loss[0].second, loss[1].second, loss[2].second, loss[3].second, loss[4].second, t);
        for (int i = 0; i < 50; i++) {
            double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
            double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();
            testSet[i] = testSort[loss[(int)i/10].second];
            if (t < 998) applyMut(testSet[i], mx, my, sx, sy);
        }
        t++;
    }while(t < 200);

    printf("\n*\n*\tBeam Position <x>      : %1.5f\n", mean(testSet[loss[0].second], 0));
    printf("*\tBeam Position <y>      : %1.5f\n", mean(testSet[loss[0].second], 1));
    printf("*\tBeam Standard deviation: %1.5f\n", standard_deviation(testSet[loss[0].second], 0));
    printf("*\tBeam Standard deviation: %1.5f\n*\n", standard_deviation(testSet[loss[0].second], 1));
    printf("*\tLoss value of data   : %1.5f\n", loss[loss[0].second].first);

    printOut(testSet[loss[0].second], testName);

    return 0;

}

void mutation(Beam_t &_input, int i, double mx, double my, double sx, double sy) {
    double mrate = 0.01;
    if (mrate > genrand_real1()) {
        _input(0, i) = rand_normal(mx, sx);
        _input(1, i) = rand_normal(my, sy);
    } 
}

void applyMut(Beam_t &_input, double mx, double my, double sx, double sy) {
    int i, j;
    std::vector<std::thread> threads;
    for (i = 0; i < _input.rows(); i++) {
        threads.emplace_back(
            [i, &_input, mx, my, sx, sy]{mutation(_input, i, mx, my, sx, sy);}
        );
    }
    for(std::thread &th:threads) {
            th.join();
    }
    return;
}

void oper(Beam_t &_input, int i, int c, double d) {
    double x = _input(0, i);
    double y = _input(1, i);
    if (c == 0) x = x + d;
    else if (c == 1) x = x - d;
    else if (c == 2) x = x * (1 + d);
    else if (c == 3) x = x * (1 - d);
    else if (c == 2) y = y + d;
    else if (c == 3) y = y - d;
    else if (c == 6) y = y * (1 + d);
    else if (c == 7) y = y * (1 - d);
    _input(0, i) = x;
    _input(1, i) = y;
    return;
}

void applyOper(Beam_t &_input, int c, double d) {
    int i, j;
    std::vector<std::thread> threads;
    for (i = 0; i < _input.rows(); i++) {
        //c = 16 * genrand_real1();
        threads.emplace_back(
            [i, &_input, c, d]{oper(_input, i, c, d);}
        );
    }
    for(std::thread &th:threads) {
            th.join();
    }
    return;
}

double mean(MatrixXd _input, int xy) {
    int i;
    double res = 0.0;
    for (i = 0; i < _input.cols(); i++) res += _input(xy, i);
    return res/_input.cols();
}

double meansqure(MatrixXd _input, int xy) {
    int i;
    double res = 0.0;
    for (i = 0; i < _input.cols(); i++) res += std::pow(_input(xy, i), 2);
    return res/_input.cols();
}

double variance(MatrixXd _input, int xy) {
    return meansqure(_input, xy) + std::pow(mean(_input, xy), 2);
}

double standard_deviation(MatrixXd _input, int xy) {
    double res = std::sqrt(variance(_input, xy));
    return res;
}

double mse(MatrixXd real, MatrixXd test) {
    int i;
    double res = 0;
    for (i = 0; i < 16; i++) {
        res += std::pow((real(i, 0) - test(i, 0)), 2);
    }
    res /= (double)16;
    return std::sqrt(res);
}

double ep(MatrixXd real, MatrixXd test) {
    double w=1.0e-3, res=0.0, norm, skew;
    norm = exp(w*std::pow(test(0, 0) - real(0, 0),2));
    skew = exp(w*std::pow(test(1, 0) - real(1, 0),2));
    res += norm + skew;
    return res-2;
}

void genSet(Beam_t &part, double mx, double my, double sx, double sy) {
    int i;
    for (i = 0; i < part.cols(); i++) {
        part(0, i) = rand_normal(mx, sx);
        part(1, i) = rand_normal(my, sy);
    }
}

double Uniform(int x0) {
	static int x=x0;
	int a=1103515245,b=12345,c=2147483647;
	x = (a*x + b)&c;
	return ((double)x+1.0) / ((double)c+2.0);
}

double rand_normal(double mu, double sigma){
	double z=sqrt( -2.0*log((double) Uniform(1)) ) * sin( 2.0*M_PI*(double) Uniform(1) );
	return mu + sigma*z;
}

MatrixXd func(MatrixXd particle) {
    int i, j, k;
    double theta = 0.125 * M_PI;
    double ans;
    MatrixXd buf = MatrixXd::Constant(16, 1, 0);
    MatrixXd vol = MatrixXd::Constant(16, 1, 0);
    for (i = 0; i < 1000; i++) {
        std::complex<double> pos(particle(0, i), particle(1, i));
        for (j = 0; j < 16; j++) {
            if (j == 0) {
                buf(j, 0) = 1;
            } else if (j%2 == 0) {
                buf(j, 0) += std::imag(pow(pos, (int)(j+1)/2));
            } else {
                buf(j, 0) += std::real(pow(pos, (int)(j+1)/2));
            }
        }
    }
    for (i = 0; i < 16; i++) buf(i, 0) /= 1000;
    for (i = 0; i < 16; i++) {
        if (i%2 == 0) {
            vol(i, 0) += std::sin(theta * (i+1)/2) * buf(i, 0);
        } else {
            vol(i, 0) += std::cos(theta * (i+1)/2) * buf(i, 0);
        }
    }
    //vol = vol.normalized();
    return vol;
}


/* initializes mt[MT_N] with a seed */
void init_genrand(unsigned long s) {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MT_N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length) {
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (MT_N>key_length ? MT_N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MT_N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MT_N) { /* generate N words at one time */
        int kk;

        if (mti == MT_N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<MT_N-MT_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MT_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void) {
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void) {
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

void printOut(Beam_t particle, char fName[]) {
    int i;
    std::string pos;
    std::ofstream posFile;
    posFile.open(fName);
    for (i = 0; i < particle.cols(); i++) {
        pos = std::to_string(particle(0, i)) + "," + std::to_string(particle(1, i)) + "\n";
        posFile << pos.c_str();
    }
    posFile << "-82.5, -82.5\n";
    posFile << "-82.5, 82.5\n";
    posFile << "82.5, -82.5\n";
    posFile << "82.5, 82.5\n";
}