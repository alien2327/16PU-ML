#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <ctime>
#include <cmath>

using namespace std;

// Defining for Numbers
#define N_Particle 100
#define N_Gen 1000

double matrix[16][16] = {{6.306873e-02,1.485052e-03,1.744601e-06,1.746443e-05,8.816988e-09,2.046298e-07,4.066876e-10,2.414934e-09,4.219697e-12,2.907648e-11,3.852436e-14,3.482063e-13,6.995238e-16,4.347983e-15,5.753933e-18,4.694390e-17},
                        {6.272149e-02,1.369322e-03,5.687868e-04,1.235001e-05,1.233067e-05,7.890840e-08,1.896040e-07,-7.845883e-12,2.411497e-09,-1.135950e-11,2.606194e-11,-2.450212e-13,2.359561e-13,-3.906964e-15,1.535925e-15,-4.380312e-17},
                        {6.347593e-02,1.055858e-03,1.059183e-03,-1.409663e-08,1.757272e-05,-1.467625e-07,1.464097e-07,-2.429393e-09,-6.188355e-12,-1.984063e-11,-1.996585e-11,6.654149e-16,-3.347436e-13,2.930826e-15,-2.929838e-15,4.424168e-17},
                        {6.186982e-02,5.577365e-04,1.355161e-03,-1.225982e-05,1.215804e-05,-1.876458e-07,-7.911695e-08,7.107163e-12,-2.403030e-09,2.613247e-11,-1.123816e-11,2.444583e-13,2.380597e-13,-1.601293e-15,3.931849e-15,-4.475268e-17},
                        {6.295156e-02,-3.355237e-06,1.487105e-03,-1.748347e-05,-6.824366e-08,5.708346e-10,-2.049700e-07,2.420107e-09,9.088590e-13,9.485498e-14,2.902655e-11,-3.481846e-13,9.568604e-16,-1.725551e-17,-4.337726e-15,4.617208e-17},
                        {6.268608e-02,-5.728343e-04,1.374223e-03,-1.238438e-05,-1.248725e-05,1.913761e-07,-7.880946e-08,-1.460045e-11,2.432732e-09,-2.625178e-11,-1.163502e-11,2.481536e-13,-2.363464e-13,1.557587e-15,3.957337e-15,-4.516764e-17},
                        {6.171772e-02,-1.035720e-03,1.032172e-03,6.390149e-08,-1.727472e-05,1.432045e-07,1.448009e-07,-2.385763e-09,-1.682846e-11,1.976437e-11,-1.946149e-11,-2.713607e-15,3.309071e-13,-2.870683e-15,-2.894161e-15,4.243431e-17},
                        {6.160264e-02,-1.353672e-03,5.591256e-04,1.226994e-05,-1.223342e-05,-7.931822e-08,1.884919e-07,9.563228e-12,-2.412664e-09,1.122314e-11,2.621865e-11,-2.438039e-13,-2.386122e-13,3.934523e-15,1.598021e-15,-4.556063e-17},
                        {6.243357e-02,-1.477396e-03,8.603782e-08,1.739635e-05,-2.521083e-08,-2.041444e-07,2.023270e-10,2.410341e-09,-1.801238e-12,-2.884369e-11,-2.441107e-14,3.455714e-13,3.096349e-16,-4.314602e-15,4.570211e-19,4.577892e-17},
                        {6.235479e-02,-1.368343e-03,-5.654166e-04,1.238001e-05,1.230077e-05,-7.979085e-08,-1.896123e-07,3.448732e-12,2.423293e-09,1.143588e-11,-2.638135e-11,-2.481000e-13,2.401398e-13,3.991085e-15,-1.590821e-15,-4.561604e-17},
                        {6.297562e-02,-1.052527e-03,-1.051931e-03,-3.247740e-09,1.748653e-05,1.458952e-07,-1.453789e-07,-2.414164e-09,-8.666592e-12,1.976141e-11,1.990301e-11,1.078097e-15,-3.346480e-13,-2.952790e-15,2.936374e-15,4.387060e-17},
                        {6.295844e-02,-5.706842e-04,-1.379649e-03,-1.248473e-05,1.238961e-05,1.908217e-07,8.093580e-08,1.552000e-11,-2.443413e-09,-2.655498e-11,1.127987e-11,2.466167e-13,2.411239e-13,1.576993e-15,-3.955167e-15,-4.495202e-17},
                        {6.228738e-02,7.492344e-07,-1.470255e-03,-1.729346e-05,-7.361513e-08,-1.257991e-09,2.029695e-07,2.398431e-09,2.003401e-11,3.807977e-13,-2.885088e-11,-3.464824e-13,-5.142345e-15,-7.151004e-17,4.321349e-15,4.608714e-17},
                        {6.215899e-02,5.640412e-04,-1.361035e-03,-1.227517e-05,-1.234533e-05,-1.896113e-07,7.830328e-08,-1.261822e-11,2.412533e-09,2.614051e-11,1.145343e-11,2.464207e-13,-2.353931e-13,-1.541729e-15,-3.947440e-15,-4.469955e-17},
                        {6.291075e-02,1.049754e-03,-1.049125e-03,4.153612e-08,-1.751984e-05,-1.455399e-07,-1.463443e-07,-2.417844e-09,-1.146205e-11,-1.994118e-11,1.968727e-11,-2.963624e-15,3.343593e-13,2.885355e-15,2.952590e-15,4.370714e-17},
                        {6.182623e-02,1.352029e-03,-5.601695e-04,1.223160e-05,-1.221897e-05,7.863085e-08,-1.877706e-07,2.907797e-12,-2.401350e-09,-1.121401e-11,-2.611385e-11,-2.438476e-13,-2.376769e-13,-3.947036e-15,-1.577811e-15,-4.467018e-17}};

double crate = 0.20;
double mrate = 0.01;

// Data reading part
int samplelength;
void read_data(char [], char [], int **, int **);
void extract_momentum(double **, double **, double *, double *, int, int, char);
void dot_matrix(double **, double **, char, char);
void goertzel_method(int **, double **, double, double);

// Evalution function
double lossfunction(double *, double *);
double mse(double *, double *);
double ep(double *, double *);

double Uniform(int);
double rand_normal(double , double);
void oper(double *, double, short);
void initPos(double [][2]);
void func0(double [2], double [16]);
void func1(double [2], double [16]);
void mutation(double [][2]);
void roulette(double *, double *);
void select(double [][2], double *);

double mean(double *, int);
double variance(double *, int);
double covariance(double *, double *, int);
double standard_deviation(double *, int);
double skewness(double *, int);
double kurtosis(double *, int);

void printStat(double p[][2], char fName[20]) {
    int i;
    double x[N_Particle], y[N_Particle];
    string pos;
    ofstream posFile;
    posFile.open(fName);
    for (i = 0; i < N_Particle; i++) {
        x[i] = p[i][0]; y[i] = p[i][1];
        pos = to_string(p[i][0]) + "," + to_string(p[i][1]) + "\n";
        posFile << pos.c_str();
    }
    printf("Mean x: %5.5e,  Mean y: %5.5e\n", mean(x, (int)N_Particle), mean(y, (int)N_Particle));
    printf(" Var x: %5.5e,   Var y: %5.5e\n", variance(x, (int)N_Particle), variance(y, (int)N_Particle));
    printf(" Std x: %5.5e,   Std y: %5.5e\n", standard_deviation(x, (int)N_Particle), standard_deviation(y, (int)N_Particle));
    printf("Skew x: %5.5e,  Skew y: %5.5e\n", skewness(x, (int)N_Particle), skewness(y, (int)N_Particle));
    printf("Kurt x: %5.5e,  Kurt y: %5.5e\n", kurtosis(x, (int)N_Particle), kurtosis(y, (int)N_Particle));
}

void find_min(double *z, int *zmin, int size)
{
    int i;
    double buf[2]={0,0};
    for (i=0;i<size;i++)
    {
        buf[0] = z[i];
        if (buf[0] < buf[1]) *zmin = i;
        buf[1] = buf[0];
    }
    return;
}

int main()
{
    Uniform(time(NULL));
    int i, j, k, l, m;
    char *file13 = "./wave13.dat", *file15 = "./wave15.dat";
    char finit[20] = "init.csv", ffinal[20] = "final.csv";
    int **buf13 = new int*[16];
    int **buf15 = new int*[16];
    double zmax = -10, global_zmax = -10;
    double **voltage_13 = new double*[16];
    double **voltage_15 = new double*[16];
    double moment_real_13[16];
    double moment_real_15[16];
    read_data(file13, file15, buf13, buf15);
    goertzel_method(buf13, voltage_13, 2, 0);
    goertzel_method(buf15, voltage_15, 2, 0);
    extract_momentum(voltage_13, voltage_15, moment_real_13, moment_real_15, 1, 0, 'g');
    double proton[N_Particle][2], z_part[N_Particle], cprob[N_Particle], v0, v13_0;
    v13_0 = voltage_13[0][0];
    for (l=0;l<16;l++)
    {
        voltage_13[0][l] /= v13_0;
    }
    initPos(proton);
    printStat(proton, finit);
    cout << endl;
    for (m=0;m<N_Gen;m++)
    {
        double mom_try[16];
        for (i=0;i<N_Particle;i++)
        {
            printf("\r %d-st Generation %d-st Particle processing\t", m, i);
            double v_buf[16], m_buf[16], p_buf[2], z_buf[2];
            p_buf[0] = proton[i][0]; p_buf[1] = proton[i][1];
            func0(p_buf, v_buf); func1(p_buf, m_buf);
            v0 = v_buf[0];
            for (l=0;l<16;l++)
            {
                v_buf[l] /= v0;
            }
            z_buf[0] = mse(moment_real_13, m_buf) + ep(moment_real_13, m_buf);
            oper(&p_buf[0], p_buf[1], (int)8*Uniform(1)); oper(&p_buf[1], p_buf[0], (int)8*Uniform(1));
            func0(p_buf, v_buf); func1(p_buf, m_buf);
            v0 = v_buf[0];
            for (l=0;l<16;l++)
            {
                v_buf[l] /= v0;
            }
            z_buf[1] = mse(moment_real_13, m_buf) + ep(moment_real_13, m_buf);
            if (z_buf[0] > z_buf[1])
            {
                proton[i][0] = p_buf[0]; proton[i][1] = p_buf[1];
                z_part[i] = 1/z_buf[1];
                break;
            }
            //cout << mse(voltage_13[0], v_buf) << "\t" << ep(moment_real_13, m_buf) << endl;
            //cout << z_part[i] << endl;
        }
        roulette(z_part, cprob);
        select(proton, cprob);
        mutation(proton);
    }
    cout << endl; cout << endl;
    printStat(proton, ffinal);
    cout << endl;
    double m_init, mom[16];
    for (i=0;i<N_Particle;i++)
    {
        double m_buf[16];
        func1(proton[i], m_buf);
        for (j=0;j<16;j++)
        {
            mom[j] += m_buf[j];
        }
    }
    m_init = mom[0];
    for (j=0;j<16;j++)
    {
        mom[j] /= m_init;
    }
    for (i=0;i<16;i++)
    {
        printf("%d\t%5.3e\t%5.3e\n", i, moment_real_13[i], mom[i]);
    }

    return 0;
}

void select0(double p[][2], double *val) {
    double newch[N_Particle][2];
    int i, j=0, ns[N_Particle/10];
    double buf[2]={0,0};
    for (i=0;i<N_Particle;i++)
    {
        buf[0] = val[i];
        if (buf[0] < buf[1]) 
        {
            ns[j] = i;
            j++;
            if (j==N_Particle/10) break;
        }
        buf[1] = buf[0];
    }

    for (i=0;i<N_Particle;i++)
    {
        for (j=0;j<2;j++) 
        {
            newch[i][j] = p[ns[(int)(i/10)]][j];
        }
    }
    for (i=0;i<N_Particle;i++)
    {
        for (j=0;j<2;j++) 
        {
            p[i][j] = newch[i][j];
        }
    }
    return;
}

void select(double p[][2], double *val) {
    double newch[N_Particle][2];
    int i, j, ns[N_Particle];
    float rnd0;
    for (i = 0; i < N_Particle; i++) {
        rnd0 = Uniform(1);
        for (j = 0; j < N_Particle; j++) {
            if (*(val+j) > (double)rnd0) {
                ns[i] = j;
                break;
            }
        }
    }
    for (i=0;i<N_Particle;i++)
    {
        for (j=0;j<2;j++) 
        {
            newch[i][j] = p[ns[i]][j];
        }
    }
    for (i=0;i<N_Particle;i++)
    {
        for (j=0;j<2;j++) 
        {
            p[i][j] = newch[i][j];
        }
    }
    return;
}

void roulette(double *val, double *q) {
    int i;
    double sum;
    sum = 0.0;
    for (i = 0; i < N_Particle; i++) {
        sum += *(val + i);
    }
    *val /= sum;
    *q = *val;
    for (i = 1; i < N_Particle; i++) {
        *(val + i) /= sum;
        *(q + i) = *(q + i - 1) + *(val + i);
    }
    return;
}

void mutation(double p[][2])
{   
    double posMax = 80, posMin = -80;
    for (int i=0;i<N_Particle;i++)
    {
        if (Uniform(1) <= mrate)
        {
            
            p[i][0] = posMin + (posMax - posMin) * Uniform(1); // Position of X
        }
        if (Uniform(1) <= mrate)
        {
            p[i][1] = posMin + (posMax - posMin) * Uniform(1); // Position of Y
        }
    }
    return;
}

void initPos(double p[][2])
{
    cout << "Initializing particle position" << endl;
    int i, j;
    double posMax = 80, posMin = -80;
    for (i=0;i<N_Particle;i++)
    {
        p[i][0] = posMin + (posMax - posMin) * Uniform(1); // Position of X
        p[i][1] = posMin + (posMax - posMin) * Uniform(1); // Position of Y
    }
    return;
}

void func0(double p[2], double v[16]) {
    int i, j;
    double ans;
    for (i = 0; i < 16; i++) {
        ans = 0;
        for (j = 0; j < 8; j++) {
            complex<double> pos(p[0], p[1]);
            if (j == 0) {
                ans += matrix[i][j];
            } else if (j == 7) {
                ans += matrix[i][2*j-1]*real(pow(82.5, -1*j) * pow(pos, j));
                ans += matrix[i][2*j]*imag(pow(82.5, -1*j) * pow(pos, j));
                ans += matrix[i][15]*real(pow(82.5, -8) * pow(pos, 8));
            } else {
                ans += matrix[i][2*j-1]*real(pow(82.5, -1*j) * pow(pos, j));
                ans += matrix[i][2*j]*imag(pow(82.5, -1*j) * pow(pos, j));
            }
        }
        v[i] = ans;
    }
}

void func1(double p[2], double mom[16]) {
    int i, j;
    double ans;
    complex<double> pos(p[0], p[1]);
    mom[0] = 1.0;
    mom[1] = real(pow(pos, 1));
    mom[2] = imag(pow(pos, 1));
    mom[3] = real(pow(pos, 2));
    mom[4] = imag(pow(pos, 2));
    mom[5] = real(pow(pos, 3));
    mom[6] = imag(pow(pos, 3));
    mom[7] = real(pow(pos, 4));
    mom[8] = imag(pow(pos, 4));
    mom[9] = real(pow(pos, 5));
    mom[10] = imag(pow(pos, 5));
    mom[11] = real(pow(pos, 6));
    mom[12] = imag(pow(pos, 6));
    mom[13] = real(pow(pos, 7));
    mom[14] = imag(pow(pos, 7));
    mom[15] = real(pow(pos, 8));
    return;
}

double Uniform(int x0) {
	static int x=x0;
	int a=1103515245,b=12345,c=2147483647;
	x = (a*x + b)&c;
	return ((double)x+1.0) / ((double)c+2.0);
}

double rand_normal(double mu, double sigma){
	double z=sqrt( -2.0*log(Uniform(1)) ) * sin( 2.0*M_PI*Uniform(1) );
	return mu + sigma*z;
 }

 void oper(double *x, double y, short c)
{
    int i, j;
    double d=1.0e-2;
    if (c==1) *x = *x + d;
    else if (c==2) *x = *x + (-1)*d;
    else if (c==3) *x = *x*(1+d);
    else if (c==4) *x = *x*(1+(-1)*d);
    else if (c==5) *x = *x*(1+y*d);
    else if (c==6) *x = *x*(1+y*(-1)*d);
    else if (c==7) *x = *x*(1+abs(y)*d);
    else if (c==8) *x = *x*(1+abs(y)*(-1)*d);
    return;
}

double mse(double *signal, double *simul)
{
    int i;
    double res=0;
    for (i = 0; i < 16; i++)
    {
        res += (signal[i] - simul[i])*(signal[i] - simul[i]);
    }
    res /= (double)16;
    return res;
}

double ep(double *signal, double *simul)
{
    int i=2, ii, iii;
    double w=1.0e-5, res=0.0, norm, skew;
    //ii = i-1;
    //iii = i;
    //norm = exp(w*pow((simul[ii] - signal[ii]),2));
    //skew = exp(w*pow((simul[iii] - signal[iii]),2));
    //res += norm + skew;
    ii = 2*i-1;
    iii = 2*i;
    norm = exp(w*pow((simul[ii] - signal[ii]),2));
    skew = exp(w*pow((simul[iii] - signal[iii]),2));
    res += norm + skew;
    return res;
}

double lossfunction(double *signal, double *simul)
{
    return 1.0e+3*mse(signal, simul);//+ep(signal, simul));
}

double mean(double *set, int size)
{
    double res = 0.0;
    for (int i=0;i<size;i++)
    {
        res += set[i];
    }
    res /= size;
    return res;
}

double variance(double *set, int size)
{
    double res = 0.0, avg = mean(set, size);
    for (int i=0;i<size;i++)
    {
        res += (set[i] + avg)*(set[i] + avg);
    }
    res /= size;
    return res;
}

double covariance(double *setA, double *setB, int size)
{
    double res = 0.0;
    double cov[52];
    for (int i=0;i<size;i++)
    {
        cov[i] = (setA[i] - mean(setA, size)) * (setB[i] - mean(setB, size));
    }
    res = mean(cov, size);
    return res;
}

double standard_deviation(double *set, int size) {
    double res = sqrt(variance(set, size));
    return res;
}

double skewness(double *set, int size)
{
    int i;
    double res, avg = mean(set, size), r = (double)size/(((double)size-1)*((double)size-2)), t;
    for (i=0;i<size;i++)
    {
        t = pow((set[i]-avg)/standard_deviation(set, size), 3);
        res += t;
    }
    res *= r;
    return res;
} 

double kurtosis(double *set, int size)
{
    int i;
    double r1, r2, t, res, avg = mean(set, size);
    r1 = ((double)size*((double)size+1))/(((double)size-1)*((double)size-2)*((double)size-3));
    r2 = (3*((double)size-1)*((double)size-1))/(((double)size-2)*((double)size-3));
    for (i=0;i<size;i++)
    {
        t = pow(set[i]-avg, 4)/pow(standard_deviation(set, size), 4);
        res += t;
    }
    res = r1*res - r2;
    return res;
}

void read_data(char data13[], char data15[], int ** _output13, int ** _output15)
{
/*
* COMMANTARY FOR WAVE DATA READING MODE:
*
* 16PickUp Beam Position Monitor(16PU) have a unique readout system, which contain bunch information.
* This function is for reading raw wave binary data from 16pu, and due to the limitation of ADC, the data only have just about 140 turn of bunch.
*
* The structer of wave binary data :
*   Header 19bytes -> wave_0000(year)_00(month)_00(day)_00(hour)_00(min)_00(sec)_address00(13 or 15).dat
*   Wave data [header 4bytes -> wave
*              channel number 2bytes -> 0x00 0xff ~ 0x0f 0xff; no meaning in 0xff
*              wave data 65528*2bytes -> 2bytes is one sample
*              footer 4bytes -> data] 131066btyes -> for 1 channel, and there are 16 channel datas.
*
* ISSUES
* 1) Channel index of 16PU at #13 of Main Ring(MR): Channel 10 and Channel 11 are swapped.
* 2) ADC of 16PU outputs 14 bits for 1 sample. To make easy to calculate, it has meaningless 0b00 bit so that ADC outputs 16 bits(= 2bytes).
* 3) Can this 140turns limitation be fixed with DMA(Direct Memory Access)?
*/
    ifstream fin13, fin15;
    fin13.open(data13, ios::in);
    fin15.open(data15, ios::in);
    const double gain13[16] = {1.000000000000000000e0, 9.864749932865497506e-1, 1.008554910983899822e0, 9.909215189952645941e-1, 1.016551297946559407e0, 1.004441025921927633e0, 1.023606422028664920e0, 1.007893430992459427e0,
    1.019020740443200568e0, 1.002733038130677601e0, 1.021101905384611097e0, 1.002615873152301385e0, 1.045142626872535008e0, 9.986874348446839189e-1, 1.022369446580688512e0, 9.918986411869963327e-1};
    const double gain15[16] = {1.000000000000000000e0, 9.927984352486749486e-1, 9.886965960947093901e-1, 9.844666881364172450e-1, 9.811551707487992102e-1, 9.783521372729635512e-1, 9.896903610596464729e-1, 9.869302062198921366e-1,
    9.906593504726489696e-1, 9.859330582601781856e-1, 9.980913383872332956e-1, 9.882277496998294053e-1, 1.000800169181874200e0, 9.924633771768432977e-1, 9.944628454760420233e-1, 1.001591146042800640e0};
 
    if (fin13.is_open() && fin15.is_open()) {
        cout << "File opened successfully" << endl;
        for (int j = 0; j < 16; j++) {
            int vol13[65528];
            int vol15[65528];
            fin13.seekg(25 + j*131066);
            fin15.seekg(25 + j*131066);
            char buffer13[131056], buffer15[131056];
            fin13.read(buffer13, 131056);
            fin15.read(buffer15, 131056);
            for (int i = 0; i < 65528; i++) {
                vol13[i] = gain13[j] * (((unsigned char) buffer13[0 + i*2])*64 + ((unsigned char) buffer13[1 + i*2])/4);
                vol15[i] = gain15[j] * (((unsigned char) buffer15[0 + i*2])*64 + ((unsigned char) buffer15[1 + i*2])/4);
            }
            int size = 468;
            int num = 0, mean = 0, lag = 0;
            for (int i = 0; i < 65528; i++) {
                num += vol13[i];
            }
            mean = num / 65528;
            int threshold = mean + 50;
            for (int i = 0; i < size; i++) {
                if ((vol13[i] + vol13[i+1] + vol13[i+2])/3 > threshold) {
                    lag += i-4;
                    break;
                }
            }
            samplelength = 65528 - lag;
            _output13[j] = new int[samplelength];
            for (int i = lag; i < 65528; i++) {
                _output13[j][i-lag] = vol13[i];
            }
            num = 0, mean = 0, lag = 0;
            for (int i = 0; i < 65528; i++) {
                num += vol15[i];
            }
            mean = num / 65528;
            threshold = mean + 50;
            for (int i = 0; i < size; i++) {
                if ((vol15[i] + vol15[i+1] + vol15[i+2])/3 > threshold) {
                    lag += i-4;
                    break;
                }
            }
            samplelength = 65528 - lag;
            _output15[j] = new int[samplelength];
            for (int i = lag; i < 65528; i++) {
                _output15[j][i-lag] = vol15[i];
            }
        }
        fin13.close();
        fin15.close();
        int *voltage_buf = new int[65528];
        voltage_buf = _output13[10];
        _output13[10] = _output13[11];
        _output13[11] = voltage_buf;
        return;
    } else {
        cout << "File Not Found" << endl;
        cout << "File Opening error occured" << endl;
        return;
    }
}

void dot_matrix(double ** _input_double, double ** _output, char address, char method)
{
    if (address == 13) {
        double monitor_matrix[16][16] = 
        {{9.999413998070250109e-01,9.999460435510404421e-01,1.000036939761390409e+00,9.998124202563855034e-01,1.000271488239349749e+00,9.998693920504219124e-01,1.000111480706963318e+00,9.998456629174904409e-01,9.999773275920876836e-01,1.000162081686137849e+00,9.998997463053272972e-01,1.000186489373741328e+00,9.997085318373660767e-01,1.000132905206968958e+00,9.999092115071314124e-01,1.000191401767124733e+00},
        {8.402690869148989350e+01,7.855569013495049546e+01,5.868037260691073698e+01,3.282522159378924442e+01,-7.166513112806041086e-01,-3.193109828147644436e+01,-6.147302662371219384e+01,-7.772471435037130050e+01,-8.545037779220223229e+01,-7.759246361210819032e+01,-6.062931303009958839e+01,-3.158140091037386199e+01,-1.446554977717688650e-01,3.261407689753472994e+01,5.948528506596413479e+01,7.875739354011824389e+01},
        {-3.184580464191904547e-01,3.228656153914549520e+01,5.995639198359812383e+01,7.820350744773398333e+01,8.476220176545791674e+01,7.757651271415326732e+01,6.031722793075426381e+01,3.220030029154117557e+01,4.337302640585219415e-01,-3.240834556455450866e+01,-5.952089441992973917e+01,-7.802962192320494239e+01,-8.487089293163909076e+01,-7.812708433428798571e+01,-5.990049734393068093e+01,-3.223788653676192695e+01},
        {7.137783296339693152e+03,5.112365477098272095e+03,-7.851964983994369618e+01,-5.081306374642089395e+03,-7.195729819149987634e+03,-5.035351301641623650e+03,5.802079164518411147e+01,5.085147097495299931e+03,7.223963968785255929e+03,5.052380434869687633e+03,3.093289319604604870e+01,-5.086032890206775846e+03,-7.219966991685299035e+03,-5.069979210138622875e+03,-4.370857235080535474e+00,5.135855289893223926e+03},
        {-4.128864486652737753e+01,5.087607693875425866e+03,7.136140344759595791e+03,5.099855744849646726e+03,-4.102723224737526664e+01,-5.040457598271681491e+03,-7.316189874531491114e+03,-5.057047774689844118e+03,-2.612702648095267577e+01,5.048683556148226671e+03,7.205862535473165735e+03,5.006432193559158804e+03,-2.916023524624746344e+01,-5.076237364001140122e+03,-7.202681364551679508e+03,-5.072362732333781423e+03},
        {6.049011323557938449e+05,2.396483519005203561e+05,-4.382757016334204818e+05,-5.618301118111484684e+05,-1.571629856130904727e+03,5.632854633836707799e+05,4.344627039328728570e+05,-2.369006016064683208e+05,-6.115075652575897984e+05,-2.343046354074297997e+05,4.341596774890003726e+05,5.562759611219625222e+05,-2.106216762353617469e+03,-5.657101597251122585e+05,-4.308280296919397661e+05,2.393924747899266367e+05},
        {-2.662443326703010825e+03,5.646673200289527886e+05,4.279343994867376168e+05,-2.366162533193831914e+05,-6.090826484717677813e+05,-2.330150638983670797e+05,4.439187612618173589e+05,5.646857144498244161e+05,-2.418693814792285139e+03,-5.565914658578806557e+05,-4.390198956733581144e+05,2.405489163061128347e+05,6.118203734236392193e+05,2.341006563270064944e+05,-4.343683409772803425e+05,-5.654667673793170834e+05},
        {5.087124571121729910e+07,1.119067933127190452e+06,-5.272955327068960667e+07,8.680483496539668413e+05,5.111064435668140650e+07,5.821473845092523843e+05,-5.345720835225655884e+07,1.006413324495706591e+06,5.141201192359492183e+07,7.237471585985765560e+05,-5.296820200364012271e+07,1.323752426483716350e+06,5.158727877359198779e+07,2.446428717375808046e+05,-5.230182401809504628e+07,8.855952866770752007e+05},
        {-3.838756986744244932e+05,5.241258514119739085e+07,-6.772074381322815316e+05,-5.168932041957750171e+07,-5.462655290799565846e+05,5.222136688028974086e+07,-8.780164960629242705e+05,-5.219173148805026710e+07,5.996072744596879929e+05,5.105952124851896614e+07,4.918003338072042679e+05,-5.171546453306925297e+07,3.770671803330594557e+05,5.209912062451072782e+07,-4.198797444921328570e+05,-5.208505744134169817e+07},
        {4.304301365307842255e+09,-1.599427431200844288e+09,-3.150392094739248276e+09,4.117372586963991642e+09,1.440313067957474105e+07,-4.085007020188084126e+09,3.254280087345081806e+09,1.586597720981162310e+09,-4.327028667759906769e+09,1.603631306629714966e+09,3.173473113201637268e+09,-4.082679482509153366e+09,3.369625028280577809e+07,4.095570413339972973e+09,-3.184360063912318707e+09,-1.595812367223905325e+09},
        {-2.102333120259914547e+07,4.125340881569543839e+09,-3.180174857363269329e+09,-1.615202780518994808e+09,4.320756676590899467e+09,-1.623865802434098244e+09,-3.183257980723783493e+09,4.136961603416942120e+09,-3.777085914116869122e+07,-4.054970198830155849e+09,3.158694890463657379e+09,1.568760697542835236e+09,-4.363833028166501999e+09,1.656159573613260984e+09,3.124283466705462456e+09,-4.116477574751970768e+09},
        {3.604535387839906006e+11,-2.554800833334613037e+11,-1.207473219268465191e+08,2.585142133768315430e+11,-3.620215949421724854e+11,2.553867624233882751e+11,-2.375439299573695183e+09,-2.564574463999008484e+11,3.632380267932874756e+11,-2.560974503151954346e+11,1.396151675565687180e+09,2.512340434692117004e+11,-3.650361730757618408e+11,2.594215007235646667e+11,-3.893049071341799259e+09,-2.551705562390465393e+11},
        {2.143926967192851603e+08,2.665171992605645752e+11,-3.731634033223012085e+11,2.675063575278382263e+11,-6.247197764797392488e+07,-2.626556799043166199e+11,3.794801826879353638e+11,-2.692885313205568848e+11,2.790921686277306080e+09,2.623359363631508484e+11,-3.718211762495717773e+11,2.636847362527774353e+11,-1.689168636435148001e+09,-2.652555795022042542e+11,3.749221298006658325e+11,-2.695171217885629578e+11},
        {2.937521841317168750e+13,-2.745575140157400000e+13,2.125136478914786328e+13,-1.207462281257746680e+13,5.157480427142016602e+10,1.174723629863906055e+13,-2.147299785263950781e+13,2.753898628755692188e+13,-2.945194560154652734e+13,2.721046700469135938e+13,-2.110487736917602344e+13,1.182805539971532422e+13,-7.890405831542990112e+10,-1.188735452162009375e+13,2.120539213286899609e+13,-2.754674298073397266e+13},
        {5.591433517171932220e+10,1.195678964676207422e+13,-2.114999901277765234e+13,2.748473708065037109e+13,-2.919485683621238672e+13,2.703571649562666016e+13,-2.157422929283521875e+13,1.216316281476184961e+13,-1.794486419635328064e+11,-1.169456655227083203e+13,2.090942587083592188e+13,-2.675360293263904297e+13,2.942612648068332422e+13,-2.744143813504088672e+13,2.142687804228785156e+13,-1.221366315172071094e+13},
        {1.411664770215943750e+15,-1.425581589371573500e+15,1.401387927769147000e+15,-1.424603075057991750e+15,1.396734144148960750e+15,-1.395123838367636750e+15,1.417306274271324000e+15,-1.428928753756292000e+15,1.410059836059073250e+15,-1.402477989716191000e+15,1.376925171344775250e+15,-1.376245229453973500e+15,1.400404994852157000e+15,-1.414332897004222250e+15,1.406113595484454750e+15,-1.432234529429722500e+15}};
        if (method == 'f') {
            double ** temp = new double*[16];
            for (int i = 0; i < 16; i++) {
                _output[i] = new double[samplelength/52];
                temp[i] = new double[samplelength-samplelength%52];
                for (int j = 0; j < samplelength; j++) {
                    for (int k = 0; k < 16; k++) {
                        temp[i][j] += monitor_matrix[i][k] * _input_double[k][j];
                    }
                }
            }
            double *sig, *del_normal, *del_skew, *normal, *skew;
            sig = new double[52];
            del_normal = new double[52];
            del_skew = new double[52];
            normal = new double[samplelength/52];
            skew = new double[samplelength/52];   
            for (int i = 0; i < samplelength; i++) {
                sig[i%52] = temp[0][i];
                del_normal[i%52] = temp[1][i];
                del_skew[i%52] = temp[2][i];
                if (i % 52 == 51) {
                    _output[1][(i+1)/52 - 1] = covariance(sig, del_normal, 52)/variance(sig, 52);
                    _output[2][(i+1)/52 - 1] = covariance(sig, del_skew, 52)/variance(sig, 52);
                }
            }
            delete [] _input_double;
        } else if (method == 'g') {
            for (int i = 0; i < 16; i++) {
                _output[i] = new double[samplelength/52];
                for (int j = 0; j < samplelength/52; j++) {
                    for (int k = 0; k < 16; k++) {
                        _output[i][j] += monitor_matrix[i][k] * _input_double[k][j];
                    }
                }
            }
            //delete [] _input_double;
        }
    } else if (address == 15) {
        double monitor_matrix[16][16] =
        {{1.000003842111640173e+00,9.999339155758250053e-01,9.999613485487530706e-01,9.999245728760640661e-01,1.000169133027317470e+00,9.999887853817465144e-01,9.999875367206880750e-01,9.999073921609406002e-01,9.999856232906959574e-01,1.000073356874412234e+00,1.000017455261239752e+00,1.000071559698644696e+00,9.998422062121793408e-01,1.000044355136718854e+00,9.999446488560410629e-01,1.000140151052085535e+00},
        {8.380796005233882795e+01,7.771145478628254466e+01,5.962960636531391856e+01,3.198826049418720885e+01,4.229671786544587087e-01,-3.218673900362998808e+01,-5.958301014935136664e+01,-7.845578013674658280e+01,-8.456861448758583322e+01,-7.793976855751184019e+01,-6.089376232309970760e+01,-3.217744786098418075e+01,-3.785840078649872070e-02,3.238496513906805774e+01,6.003392960722874960e+01,7.779868772824023893e+01},
        {9.622955467608648172e-01,3.199108449008893373e+01,6.197620014722254211e+01,7.764073970663888247e+01,8.666517794977795575e+01,7.749201832231049991e+01,6.112833626807520915e+01,3.181451876045507277e+01,7.090234546601680288e-01,-3.312628849943363463e+01,-5.958684108496061782e+01,-7.806307359498298126e+01,-8.412913589212907084e+01,-7.821944695427288252e+01,-5.922829736390616517e+01,-3.306211160653684544e+01},
        {7.018994189263068620e+03,5.080357160639074209e+03,-1.837478689975672523e+02,-5.080218438777588744e+03,-7.366222460595396115e+03,-5.062213150046626652e+03,-6.330767817191470925e+01,5.136716320112727772e+03,7.194531204534600874e+03,5.060009952214029909e+03,4.143149338333218168e+01,-5.025487138487653283e+03,-7.214379876305093603e+03,-5.014752753263587692e+03,-2.908876472634371169e+01,5.075261854226389005e+03},
        {4.719073310244203157e+01,5.068121272655576831e+03,7.272526950383818985e+03,5.026391884011711227e+03,3.344636575995555461e+01,-5.056848282609405942e+03,-7.228762161704264145e+03,-5.107356320504989526e+03,1.881793136202152184e+01,5.064331720763235353e+03,7.364170437630719789e+03,5.003562853610916136e+03,2.785927761432990835e+01,-5.070438958997464397e+03,-7.198569966895035577e+03,-5.081946293932432127e+03},
        {5.931658195265893592e+05,2.338815616208376305e+05,-4.476826783846946782e+05,-5.602819832985377871e+05,-1.639226467128859667e+03,5.638122697122534737e+05,4.375989044537804439e+05,-2.385500737585728930e+05,-6.142313662606854923e+05,-2.341003643493511190e+05,4.422583155467043398e+05,5.564866324803766329e+05,1.553429189010804066e+03,-5.593735328294278588e+05,-4.361454774954200839e+05,2.352126008516644652e+05},
        {1.733590102775183823e+03,5.636934467244467232e+05,4.260196134581808583e+05,-2.405033210676578165e+05,-6.225448630826633889e+05,-2.370911855171922653e+05,4.360538632523380802e+05,5.695775586681684945e+05,-3.431410642116015424e+03,-5.617377246291603660e+05,-4.498123185435839114e+05,2.349160907993594592e+05,6.111627218046365306e+05,2.317092636066973209e+05,-4.345437074614539160e+05,-5.605178167293900624e+05},
        {4.987092356801528484e+07,4.447518450176328188e+05,-5.295405933303615451e+07,1.458359062735761749e+06,5.220777192180593312e+07,9.405731605342011899e+05,-5.344339141771052033e+07,9.279005807395989541e+05,5.188593210407877713e+07,6.072565646811689949e+05,-5.424851930740482360e+07,8.082962582123586908e+05,5.162304938693358004e+07,3.637025648001093650e+05,-5.258282670158503950e+07,6.282779393843616126e+05},
        {1.977806225841316627e+05,5.169505295087503642e+07,-1.130016897219566396e+06,-5.212425412875749916e+07,1.357198185928351886e+05,5.213947339811725169e+07,-8.250188782278122380e+04,-5.270661495766038448e+07,2.963350411871391698e+05,5.208375103712174296e+07,2.079337748516576830e+05,-5.132615877065725625e+07,-2.330371215148026240e+05,5.161127541797751188e+07,-8.508464904555617250e+04,-5.154316436626729369e+07},
        {4.202112711203301430e+09,-1.634195666702368498e+09,-3.106909989698691845e+09,4.149427056621729374e+09,-2.202770379153316468e+07,-4.123984889188008308e+09,3.228289590372787476e+09,1.621709320575987101e+09,-4.400534968101716042e+09,1.639161894067622900e+09,3.265042470884566307e+09,-4.060240162135756016e+09,-8.020696722210187465e+06,4.041783743209918499e+09,-3.155110322165650368e+09,-1.596973679787877798e+09},
        {3.388537066449384391e+07,4.030980185253132343e+09,-3.205676374572342873e+09,-1.580978307620392799e+09,4.430568934738506317e+09,-1.631733733018624544e+09,-3.193446903317361832e+09,4.156792061131393433e+09,-4.499298766080981120e+06,-4.133607259863242626e+09,3.253068793993239403e+09,1.589668265389823675e+09,-4.350170441346466064e+09,1.605656955020836115e+09,3.176162853574741364e+09,-4.074707600644390583e+09},
        {3.491852905899704590e+11,-2.530942426647469788e+11,3.492940201083008766e+09,2.569342206703104248e+11,-3.726711443960348511e+11,2.604355857435857544e+11,-3.275830177691141605e+09,-2.585089495585712585e+11,3.693232092650638428e+11,-2.597955839410419617e+11,-6.870804176606903076e+08,2.535960551532515564e+11,-3.641042681879528198e+11,2.526051257965307312e+11,6.335898619488921165e+08,-2.521209547659824524e+11},
        {3.420532074989078522e+09,2.576447819911822815e+11,-3.709081933430524902e+11,2.707143839164924622e+11,-2.477324968804492950e+09,-2.664755235872942810e+11,3.794037234303944092e+11,-2.698843546482242737e+11,-6.309554232246243395e+06,2.685989224207263184e+11,-3.850580728768602905e+11,2.633616879409482117e+11,5.472732178109289408e+08,-2.612775108673136902e+11,3.735300257592744141e+11,-2.641227410725647278e+11},
        {2.793837643239760938e+13,-2.656053689846369531e+13,2.106818956430916797e+13,-1.220541616123223242e+13,1.410053179974940491e+11,1.192487904183614062e+13,-2.146026365579735938e+13,2.771347116307352734e+13,-3.001522284013966016e+13,2.777374808534591797e+13,-2.191481854220090625e+13,1.191792536429475586e+13,-8.449977767506344604e+10,-1.155660444583321289e+13,2.093750252803382422e+13,-2.663338881366574219e+13},
        {1.559331277284637146e+11,1.149561165512746289e+13,-2.090364987720534766e+13,2.756545234813671484e+13,-3.019842706768454297e+13,2.770540350745921875e+13,-2.176459105054109766e+13,1.225192226734455664e+13,-6.997186266261053467e+10,-1.204215544655462500e+13,2.185911061233155469e+13,-2.703374156063121875e+13,2.938449759477989062e+13,-2.673242217819759375e+13,2.100472118451342969e+13,-1.179704810393715625e+13},
        {1.323100685038486250e+15,-1.360341993819534250e+15,1.376796562643399500e+15,-1.423658727584073000e+15,1.440792882199037250e+15,-1.428285215714143250e+15,1.423892139812061250e+15,-1.434509674775447500e+15,1.435978372553446750e+15,-1.436521328263131500e+15,1.443796679699884000e+15,-1.397169795814290500e+15,1.397732743069400000e+15,-1.367627397169950500e+15,1.367615245144060750e+15,-1.367057611243086750e+15}};
        if (method == 'f') {
            double ** temp = new double*[16];
            for (int i = 0; i < 16; i++) {
                _output[i] = new double[samplelength/52];
                temp[i] = new double[samplelength-samplelength%52];
                for (int j = 0; j < samplelength; j++) {
                    for (int k = 0; k < 16; k++) {
                        temp[i][j] += monitor_matrix[i][k] * _input_double[k][j];
                    }
                }
            }
            double *sig, *del_normal, *del_skew, *normal, *skew;
            sig = new double[52];
            del_normal = new double[52];
            del_skew = new double[52];
            normal = new double[samplelength/52];
            skew = new double[samplelength/52];
            for (int i = 0; i < samplelength; i++) {
                sig[i%52] = temp[0][i];
                del_normal[i%52] = temp[1][i];
                del_skew[i%52] = temp[2][i];
                if (i % 52 == 51) {
                    _output[1][(i+1)/52 - 1] = covariance(sig, del_normal, 52)/variance(sig, 52);
                    _output[2][(i+1)/52 - 1] = covariance(sig, del_skew, 52)/variance(sig, 52);
                }
            }
            delete [] _input_double;
        } else if (method == 'g') {
            for (int i = 0; i < 16; i++) {
                _output[i] = new double[samplelength/52];
                for (int j = 0; j < samplelength/52; j++) {
                    for (int k = 0; k < 16; k++) {
                        _output[i][j] += monitor_matrix[i][k] * _input_double[k][j];
                    }
                }
            }
            //delete [] _input_double;
        }
    }
}

void extract_momentum(double ** _input_double13, double ** _input_double15, double *mom13, double *mom15, int order, int bunch, char method)
{
    double ** _output13 = new double*[16];
    double ** _output15 = new double*[16];
    const double align13[2] = {0.7799954861648153, -0.17114181002093012};
    const double align15[2] = {-0.25069580086359916, 1.0363901454673603};
    dot_matrix(_input_double13, _output13, 13, method);
    dot_matrix(_input_double15, _output15, 15, method);
    for (int i = 0; i < 16; i++) {
        mom13[i] = _output13[i][0];
        mom15[i] = _output15[i][0];
    }
    double mono13=mom13[0], mono15=mom15[0];
    for (int i = 0; i < 16; i++) {
        mom13[i] /= mono13;
        mom15[i] /= mono15;
    }
    return;
}

void goertzel_method(int ** _input, double ** _output, double k, double x_1)
{
/*
* It will get rid of all other non-signal frequency.
* But we should figure out how to handle the signal when it have other oscillation.
*/
    const int num = 52;
    double omega = 2 * M_PI * k / num;
    double ans;
    double s[num];
    for (int ch = 0; ch < 16; ch++) {
        _output[ch] = new double[samplelength/52];
        for (int i = 0; i < samplelength; i++) {
            if (i == 0) {
                s[0] = _input[ch][i] + 2 * cos(omega) * x_1;
            } else if (i == 1) {
                s[1] = _input[ch][i] + 2 * cos(omega) * s[0] - x_1;
            } else {
                s[i%52] = _input[ch][i] + 2 * cos(omega) * s[i%52-1] - s[i%52-2];
            }
            if (i%52 == 51) {
                ans = sqrt(pow(s[num-1] - cos(omega) * s[num-2], 2) + pow(sin(omega) * s[num-2], 2))/(4*num*num);
                _output[ch][(i+1)/52 - 1] = ans;
            }
        }
    }
    delete [] _input;
}