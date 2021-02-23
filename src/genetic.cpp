#include "genetic.h"

void initPos(Beam_t &part) {
    int i;
    Uniform(100*genrand_real1());
    double posMax = 10, posMin = -10;
    double sigMax = 10, sigMin = 1;
    for (i = 0; i < part.cols(); i++) {
        double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
        double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();
        part(0, i) = rand_normal(mx, sx);
        part(1, i) = rand_normal(my, sy);
    }
}

void mutation(Beam_t &_input, int i, double mrate) {
    double posMax = 10, posMin = -10;
    if (mrate > genrand_real1()) _input(0, i) = posMin + (posMax - posMin) * genrand_real1();
    if (mrate > genrand_real1()) _input(1, i) = posMin + (posMax - posMin) * genrand_real1();
}

void applyMut(Beam_t &_input, double mrate) {
    int i, j;
    std::vector<std::thread> threads;
    for (i = 0; i < _input.rows(); i++) {
        threads.emplace_back(
            [i, &_input, mrate]{mutation(_input, i, mrate);}
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
    else if (c == 1) x = x + (-1) * d;
    else if (c == 2) x = x * (1 + d);
    else if (c == 3) x = x * (1 + (-1) * d);
    else if (c == 4) x = x * (1 + y * d);
    else if (c == 5) x = x * (1 + y * (-1) * d);
    else if (c == 6) x = x * (1 + std::abs(y) * d);
    else if (c == 7) x = x * (1 + std::abs(y) * (-1) * d);
    else if (c == 8) y = y + d;
    else if (c == 9) y = y + (-1) * d;
    else if (c == 10) y = y * (1 + d);
    else if (c == 11) y = y * (1 + (-1) * d);
    else if (c == 12) y = y * (1 + x * d);
    else if (c == 13) y = y * (1 + x * (-1) * d);
    else if (c == 14) y = y * (1 + std::abs(x) * d);
    else if (c == 15) y = y * (1 + std::abs(x) * (-1) * d);
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
    int i;
    double res = 0.0, avg = mean(_input, xy);
    for (i = 0; i < _input.cols(); i++) {
        res += std::pow((_input(xy, i) + avg), 2);
    }
    return res/_input.cols();
}

double standard_deviation(MatrixXd _input, int xy) {
    double res = std::sqrt(variance(_input, xy));
    return res;
}

MatrixXd func(MatrixXd particle, MatrixXd monitor) {
    int i, j, k;
    double ans;
    MatrixXd buf = MatrixXd::Constant(16, 1, 0);
    MatrixXd vol = MatrixXd::Constant(16, 1, 0);
    for (i = 0; i < particle.rows(); i++) {
        std::complex<double> pos(particle(0, i), particle(1, i));
        for (j = 0; j < 16; j++) {
            if (j%2 == 0) {
                buf(j, 0) += std::imag(pow(pos, (int)(j+1)/2));
            } else {
                buf(j, 0) += std::real(pow(pos, (int)(j+1)/2));
            }
        }
    }
    for (i = 0; i < 16; i++) buf(i, 0) /= particle.rows();
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
            vol(i, 0) += monitor(i, j) * buf(j, 0);
        }
    }
    return vol;
}

double lossmin(MatrixXd loss, MatrixXd::Index &row, MatrixXd::Index &col) {
    return loss.minCoeff(&row, &col);
}

double mse(MatrixXd real, MatrixXd test) {
    int i;
    double res = 0;
    real = real.normalized();
    test = test.normalized();
    //std::cout << "real " << real << std::endl;
    //std::cout << "test " << test << std::endl;
    for (i = 0; i < 16; i++) {
        res += std::pow((real(i, 0) - test(i, 0)), 2);
    }
    res /= (double)16;
    return res;
}

double ep(MatrixXd real, MatrixXd test) {
    double w=1.0e-3, res=0.0, norm, skew;
    norm = exp(w*std::pow(test(0, 0) - real(0, 0),2));
    skew = exp(w*std::pow(test(1, 0) - real(1, 0),2));
    res += norm + skew;
    return res;
}

void testPart(tBeam_t &part) {
    part(0,0) = 0.0; part(1,0) = 1.0;
    part(0,1) = 1.0; part(1,1) = 3.0;
    part(0,2) = 1.0; part(1,2) = -1.0;
    part(0,3) = 2.0; part(1,3) = 1.0;
    part(0,4) = 0.5; part(1,4) = -0.75;
    part(0,5) = 0.5; part(1,5) = 2.75;
    part(0,6) = 1.5; part(1,6) = 2.75;
    part(0,7) = 1.5; part(1,7) = -0.75;
}

double Uniform(int x0) {
	static int x=x0;
	int a=1103515245,b=12345,c=2147483647;
	x = (a*x + b)&c;
	return ((double)x+1.0) / ((double)c+2.0);
}

double rand_normal(double mu, double sigma){
	double z=sqrt( -2.0*log((double) Uniform(1)) ) * sin( 2.0*3.14159265358979323846*(double) Uniform(1) );
	return mu + sigma*z;
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

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void) {
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void) {
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)  { 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */