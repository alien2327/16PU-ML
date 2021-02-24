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
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);
void printOut(Beam_t, char *);

int main() {
    Uniform((unsigned)time(NULL)); init_genrand((unsigned)time(NULL));
    mat << 6.306858559892658533e-02,1.482997160229145061e-03,1.932693099737171711e-06,1.742187044074983119e-05,1.919478059537842648e-08,2.039445135809482120e-07,4.062724154098598651e-10,2.404000958475349078e-09,2.339842577559197180e-12,2.895374829352571522e-11,-4.359881092591580478e-14,3.461598848692095130e-13,-1.104536346658481312e-15,4.324116984954027844e-15,-1.877348033525889903e-17,4.667593986546279126e-17,
            6.270411513991192010e-02,1.366821678011800835e-03,5.679607562713779364e-04,1.231494223063116206e-05,1.229232589311296353e-05,7.875069673991775359e-08,1.887556336444376929e-07,-5.032083051544871113e-12,2.398880314071163594e-09,-1.123155326351647815e-11,2.592339652536403823e-11,-2.433950167470381421e-13,2.345004335895843432e-13,-3.862269341327394615e-15,1.497430046122542905e-15,-4.305385850426997581e-17,
            6.340916415576981768e-02,1.053164140103440881e-03,1.056142333341545270e-03,-5.102247399827231990e-09,1.750632343247255447e-05,-1.459828295828703233e-07,1.459220485792591677e-07,-2.417322178417649681e-09,-1.366591659752695378e-12,-1.978666818392504647e-11,-1.974489837029362573e-11,-1.326793840982442943e-15,-3.315028214666960608e-13,2.885524864304359093e-15,-2.902042619418628090e-15,4.360999108874165451e-17,
            6.190812886472765036e-02,5.578166919608529154e-04,1.354376314478419668e-03,-1.224170836782736781e-05,1.215541120510567287e-05,-1.873546610147871706e-07,-7.884246011902757104e-08,6.398421858865230233e-12,-2.394429571951320636e-09,2.597657514537112233e-11,-1.122921828858256082e-11,2.446867132584526887e-13,2.359062566675133670e-13,-1.571496106921129635e-15,3.931121286931093269e-15,-4.458594063392197194e-17,
            6.293527489299885480e-02,-3.055581947056812476e-06,1.484449106605406242e-03,-1.743895172823152083e-05,-5.619545925245574001e-08,3.538750576085582306e-10,-2.041954666132018319e-07,2.407568763179304716e-09,-1.612862869535479174e-12,1.412720481079651672e-13,2.881871834707998344e-11,-3.455177467692698062e-13,8.890278700616894760e-16,-1.313293984357184782e-17,-4.287695218658002360e-15,4.548027085292727818e-17,
            6.271618995107886008e-02,-5.722460811459761107e-04,1.373188897648665698e-03,-1.236634508401917867e-05,-1.246894666326963317e-05,1.909314740052408077e-07,-7.862721539700389796e-08,-1.334334986736940588e-11,2.424007592884255551e-09,-2.610100096836635330e-11,-1.155211769960800620e-11,2.453346346780046041e-13,-2.339943739799943139e-13,1.556870525339142641e-15,3.888011786747420304e-15,-4.445888361033005342e-17,
            6.173996654341597917e-02,-1.034981204277181464e-03,1.030900601147954941e-03,7.291690404951458104e-08,-1.724938733941066136e-05,1.427682967953107730e-07,1.446795326209211735e-07,-2.379281831949312079e-09,-1.933720641948481538e-11,1.970280978079300877e-11,-1.934703579812932572e-11,-3.411480469709581532e-15,3.290707837992892395e-13,-2.843890644252795321e-15,-2.867643904133659029e-15,4.181266591685914254e-17,
            6.162953397117861692e-02,-1.352061123963206558e-03,5.583669106594901180e-04,1.224710299581592480e-05,-1.219790806374033798e-05,-7.926872864672618489e-08,1.878168743082293204e-07,1.589687150932056608e-11,-2.401176234510510965e-09,1.101025819896856207e-11,2.609803965111637711e-11,-2.391383318966157028e-13,-2.381828916779071182e-13,3.867156757347120041e-15,1.617634770518585866e-15,-4.485664498166422574e-17,
            6.246310825667692607e-02,-1.475578539453794745e-03,-2.307143103382068324e-07,1.735399046043943551e-05,-1.270528980912599585e-08,-2.034144074767865791e-07,-4.631811706182973941e-12,2.397194637581149417e-09,-3.493338473571088572e-13,-2.861285653495743111e-11,-2.845720060263237057e-14,3.424555004263747645e-13,-2.444756449223852830e-16,-4.281178544239468023e-15,1.069035889908997740e-17,4.550293282999504462e-17,
            6.231559700086032605e-02,-1.365822532877950793e-03,-5.642011702986558883e-04,1.234471494294785183e-05,1.226632872726865122e-05,-7.945554865664676296e-08,-1.889732491921172740e-07,3.602180727006952938e-12,2.413426636539993177e-09,1.138168309736826288e-11,-2.629751158337622919e-11,-2.471595992847066223e-13,2.403410382397360571e-13,3.975519728852784444e-15,-1.613123353294505262e-15,-4.549733108474238518e-17,
            6.296191645533864845e-02,-1.050964719441855851e-03,-1.050505782292582898e-03,-2.481096037266311780e-09,1.744045567777191910e-05,1.454283671110233661e-07,-1.448056937267736226e-07,-2.405189279466958343e-09,-1.011532489020485733e-11,1.963266603018630725e-11,1.982874340169821727e-11,7.445106596646840444e-16,-3.325532212214533258e-13,-2.924742555486989283e-15,2.917085346337418889e-15,4.336500925813843653e-17,
            6.293960219485614238e-02,-5.696826547674368416e-04,-1.377570048881338917e-03,-1.245289075793739363e-05,1.236430021826216526e-05,1.901687279641771566e-07,8.064993142411434527e-08,1.653320590304075394e-11,-2.432187935751331518e-09,-2.644950050443716061e-11,1.124207275932461716e-11,2.463512756062822767e-13,2.392615453932187881e-13,1.552798788090359964e-15,-3.927842126764910569e-15,-4.453137516366950648e-17,
            6.225867012787331817e-02,6.251679917651621279e-07,-1.467849442583106035e-03,-1.725150169810810760e-05,-6.533734571797750406e-08,-1.217431450872897076e-09,2.023085069189140391e-07,2.388978480583722754e-09,2.169418608926768804e-11,4.713546377224208667e-13,-2.874413031842088451e-11,-3.450911987065095274e-13,-6.913389655318749675e-15,-9.351462155043441792e-17,4.291839424060956589e-15,4.560400855148345223e-17,
            6.219011137123432037e-02,5.636843797853755174e-04,-1.359654186325118166e-03,-1.225057543616735032e-05,-1.232491104153247124e-05,-1.890021632313294588e-07,7.798513220737458545e-08,-1.289157039813612499e-11,2.401339779947685949e-09,2.597162929569837210e-11,1.139908290871086090e-11,2.447484865030399886e-13,-2.328740860278413931e-13,-1.509560282914026318e-15,-3.919859962653870638e-15,-4.398772878494509285e-17,
            6.290489859870060374e-02,1.048322371121852619e-03,-1.047911067516170717e-03,4.305839581556852308e-08,-1.748491458085870933e-05,-1.451359765918498747e-07,-1.458580114931782318e-07,-2.410665400548597858e-09,-1.188378820603817085e-11,-1.990085746296588750e-11,1.961657800453363281e-11,-3.782473362601699055e-15,3.334389646291709219e-13,2.855967826049600475e-15,2.943241877454974319e-15,4.320179107702170763e-17,
            6.185486385506931001e-02,1.350970689472659686e-03,-5.593939816164058242e-04,1.221033213282047835e-05,-1.218349151114654657e-05,7.845288891546974036e-08,-1.870969078429822285e-07,3.683613697340263756e-12,-2.389395616940263846e-09,-1.115538855323494723e-11,-2.594178603098093507e-11,-2.416012011496970778e-13,-2.360859426102660706e-13,-3.917547754173792385e-15,-1.553285656380155576e-15,-4.426755574033102561e-17;

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

    double posMax = 20, posMin = -20;
    double sigMax = 8, sigMin = 4;

    double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
    double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();

    genSet(trueSet, mx, my, sx, sy);

    posMax = 50, posMin = -50;
    sigMax = 20, sigMin = 1;

    printOut(trueSet, trueName);

    trueVol = func(trueSet);
    trueMeanSquare(0,0) = meansqure(trueSet, 0);
    trueMeanSquare(1,0) = meansqure(trueSet, 1);

    printf("\n*\tBeam Position <x>    : %1.5f\n", mean(trueSet, 0));
    printf("*\tBeam Position <y>    : %1.5f\n", mean(trueSet, 1));
    printf("*\tBeam Mean squre <x^2>: %1.5f\n", meansqure(trueSet, 0));
    printf("*\tBeam Mean squre <y^2>: %1.5f\n*\n", meansqure(trueSet, 1));

    for (int i = 0; i < 50; i++) {
        double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
        double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();
        genSet(testSet[i], mx, my, sx, sy);
        //testVol[i] = func(testSet[i]);
        //testMeanSquare[i](0,0) = meansqure(testSet[i], 0);
        //testMeanSquare[i](1,0) = meansqure(testSet[i], 1);
        //loss.push_back(std::make_pair(mse(trueVol, testVol[i]) + ep(trueMeanSquare, testMeanSquare[i]), i));
        //testSort[i] = testSet[i];
    }
    //sort(loss.begin(), loss.end(), std::less<>());
    //for (int i = 0; i < 50; i++) {
    //    testSet[i] = testSort[loss[(int)i/10].second];
    //    applyMut(particle[i], mrate);
    //}

    int t = 0;

    do{
        for (int i = 0; i < 50; i++) {
            applyOper(testSet[i], 5*genrand_real1(), 0.01);
            testVol[i] = func(testSet[i]);
            testMeanSquare[i](0,0) = meansqure(testSet[i], 0);
            testMeanSquare[i](1,0) = meansqure(testSet[i], 1);
            loss.push_back(std::make_pair(mse(trueVol, testVol[i]) + ep(trueMeanSquare, testMeanSquare[i]), i));
            testSort[i] = testSet[i];
        }
        sort(loss.begin(), loss.end(), std::less<>());
        printf("\r*\t%d-st\t%d-st\t%d-st\t%d-st\t%d-st\tbunch selected at\t%d\ttry", loss[0].second, loss[1].second, loss[2].second, loss[3].second, loss[4].second, t);
        for (int i = 0; i < 50; i++) {
            double mx = posMin + (posMax - posMin) * genrand_real1(), my = posMin + (posMax - posMin) * genrand_real1();
            double sx = sigMin + (sigMax - sigMin) * genrand_real1(), sy = sigMin + (sigMax - sigMin) * genrand_real1();
            testSet[i] = testSort[loss[(int)i/10].second];
            applyMut(testSet[i], mx, my, sx, sy);
        }
        t++;
    }while(t < 100);

    printf("\n*\n*\tBeam Position <x>    : %1.5f\n", mean(testSet[loss[0].second], 0));
    printf("*\tBeam Position <y>    : %1.5f\n", mean(testSet[loss[0].second], 1));
    printf("*\tBeam Mean squre <x^2>: %1.5f\n", meansqure(testSet[loss[0].second], 0));
    printf("*\tBeam Mean squre <y^2>: %1.5f\n", meansqure(testSet[loss[0].second], 1));
    printf("*\tLoss value of data   : %1.5e\n", loss[loss[0].second].first);

    printOut(testSet[0], testName);

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
    if (c == 0) x = x + 100*d;
    else if (c == 1) x = x - 100*d;
    else if (c == 2) x = x * (1 + d);
    else if (c == 3) x = x * (1 - d);
    else if (c == 2) y = y + 100*d;
    else if (c == 3) y = y - 100*d;
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

double mse(MatrixXd real, MatrixXd test) {
    int i;
    double res = 0;
    for (i = 0; i < 16; i++) {
        res += std::pow((real(i, 0) - test(i, 0)), 2);
    }
    res /= (double)16;
    return res;
}

double ep(MatrixXd real, MatrixXd test) {
    double w=1.0e-4, res=0.0, norm, skew;
    norm = exp(w*std::pow(test(0, 0) - real(0, 0),2));
    skew = exp(w*std::pow(test(1, 0) - real(1, 0),2));
    res += norm + skew;
    return res;
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
    vol = vol.normalized();
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

void printOut(Beam_t particle, char fName[]) {
    int i;
    std::string pos;
    std::ofstream posFile;
    posFile.open(fName);
    for (i = 0; i < particle.cols(); i++) {
        pos = std::to_string(particle(0, i)) + "," + std::to_string(particle(1, i)) + "\n";
        posFile << pos.c_str();
    }
    posFile << "-60, -60\n";
    posFile << "-60, 60\n";
    posFile << "60, -60\n";
    posFile << "60, 60\n";
}