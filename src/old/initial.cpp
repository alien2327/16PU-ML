#define _USE_MATH_DEFINES
#define N 20
#define NBIT 20
#define MAX_GEN 100

#include <iostream>
#include <time.h>
#include <cmath>

using namespace std;

float rnd(short int);
void init(int [][NBIT]);
void func(int [][NBIT], double *, double *);
void roulette(double *, double *);
void cross(int [][NBIT]);
void select(int [][NBIT], double *);
void mutation(int [][NBIT]);
void find_max(double *, double *, double *, double *);
void showChrom(int [][NBIT]);
void showRes(int [][NBIT], double [N], double [N]);
void showPq(double *, double *);

double crate = 0.15;
double mrate = 0.01;

int main() {
    /*
    1. initialize chromosome
    2. calculate using f(x)
    3. setting roulette
    4. select chromosome
    5. cross chromosome
    6. mutate chromosome
    7. back to 2 (repeat MAX_GEN)
    */ 
    int gen, chrom[N][NBIT];
    double x[N], z[N], cprob[N];
    double xmax = 0.0, zmax = -10.0, global_xmax = 0.0, global_zmax = -10.0;
    clock_t start = clock();
    init(chrom);
    showChrom(chrom);
    for (gen = 0; gen < MAX_GEN; gen++) {
        func(chrom, x, z);
        //showRes(chrom, x, z);
        find_max(x, z, &xmax, &zmax);
        printf("   xmax=%10.6f, zmax=%10.6f  ", xmax, zmax);
        if (zmax >= global_zmax) {
            global_zmax = zmax;
            global_xmax = xmax;
        }
        printf("  gxmax=%10.6f, gzmax=%10.6f\n", global_xmax, global_zmax);
        if (gen == MAX_GEN - 1) break;
        roulette(z, cprob);
        //showPq(z, cprob);
        select(chrom, cprob);
        //showChrom(chrom);
        cross(chrom);
        //showChrom(chrom);
        mutation(chrom);
        //showChrom(chrom);
    }
    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    printf("\ntime %lf[s]\n", time);
    return 0;
}

float rnd(short int sd) {
    static short int ix = 1, init_on = 0;
    if (((sd % 2) != 0) && (init_on) == 0) {
        ix = sd;
        init_on = 1;
    }
    ix = 899 * ix;
    if (ix < 0) {
        ix = ix + 32767 + 1;
    }
    return ((float) ix/32768.0);
}

double f_x(double x) {
    /*
    Target function: y(x) = 2 + sin(2.2*pi*x) + cos(3.1*pi*x) (0.0 <= x <= 3.0)
    */
    double y;
    y = 2.0 + sin(2.2 * M_PI * x) + cos(3.1 * M_PI * x);
    return y;
}

void func(int chrm[][NBIT], double *x, double *z) {
    int i, j, dec, dec_max = 1, dec_p = 1;
    double x_min = 0.0, x_max = 3.0, x0;
    for (i = 1; i < NBIT; i++) {
        dec_p *= 2;
        dec_max += dec_p;
    }
    for (i = 0; i < N; i++) {
        dec = chrm[i][0];
        dec_p = 1;
        for (j = 1; j < NBIT; j++) {
            dec_p *= 2;
            dec += chrm[i][j] * dec_p;
        }
        x0 = x_min + (x_max - x_min) * (double)dec / (double)dec_max;
        *x++ = x0;
        *z++ = f_x(x0);
    }
    return;
}

void mutation(int chrm[][NBIT]) {
    int i, j, *adr_x;
    for (i = 0; i < N; i++) {
        for (j = 0; j < NBIT; j++) {
            if (rnd(1) <= mrate) {
                adr_x = *chrm + NBIT * i + j;
                printf("Mutation point=(%3d-th chrom, %3d-th locus )\n", i, j);
                if (*adr_x == 1) *adr_x = 0;
                else *adr_x = 1;
            }
        }
    }
    printf("\n");
    return;
}

void find_max(double *x, double *z, double *xmax, double *zmax) {
    int i;
    *zmax = -10.0;
    for (i = 0; i < N; i++) {
        if (*(z + i) > *zmax) {
            *zmax = *(z + i);
            *xmax = *(x + i);
        }
    }
    return;
}

void roulette(double *val, double *q) {
    int i;
    double sum;
    sum = 0.0;
    for (i = 0; i < N; i++) {
        sum += *(val + i);
    }
    *val /= sum;
    *q = *val;
    for (i = 1; i < N; i++) {
        *(val + i) /= sum;
        *(q + i) = *(q + i - 1) + *(val + i);
    }
    return;
}

void cross(int chrm[][NBIT]) {
    int newch[N][NBIT], ns[N];
    int i, j, nc = 0, cross_point, i1, i2, *chng1, *chng2;
    for (i = 0; i < N; i++) {
        if (rnd(1) <= crate) {
            ns[nc] = i;
            nc++;
        }
    }
    if (nc % 2 == 1) {
        ns[nc] = (int)(rnd(1) * (float)N);
        nc++;
    }
    for (i = 0; i < nc / 2; i++) {
        i1 = 2 * i;
        i2 = i1 + 1;
        printf("target chromosome %3d, %3d\n", ns[i1], ns[i2]);
        cross_point = (int)(rnd(1) * (float)(NBIT - 1)) + 1;
        printf("cross point=%3d\n", cross_point);
        for (j = 0; j < cross_point; j++) {
            chng1 = *chrm + NBIT * ns[i1] + j;
            chng2 = *chrm + NBIT * ns[i2] + j;
            newch[i1][j] = *chng1;
            newch[i2][j] = *chng2;
        }
        for (j = cross_point; j < NBIT; j++) {
            chng1 = *chrm + NBIT * ns[i1] + j;
            chng2 = *chrm + NBIT * ns[i2] + j;
            newch[i1][j] = *chng2;
            newch[i2][j] = *chng1;
        }
    }
    for (i = 0; i < nc; i++) {
        for (j = 0; j < NBIT; j++) {
            chng1 = *chrm + NBIT * ns[i] + j;
            *chng1 = newch[i][j];
        }
    }
    printf("\n");
    return;
}

void init(int chrm[][NBIT]) {
    int i, j;
    rnd((unsigned)time(NULL)); //乱数の初期化
    for (i = 0; i < N; i++) {
        for (j = 0; j < NBIT; j++) {
            if (rnd(1) <= 0.5) chrm[i][j] = 0; //乱数が0.5以下なら染色体に0
            else chrm[i][j] = 1; //乱数が0.5以上なら染色体に1
        }
    }
    return;
}

void select(int chrm[][NBIT], double *val) {
    int newch[N][NBIT], ns[N];
    int i, j, *chng;
    float rnd0;
    for (i = 0; i < N; i++) {
        rnd0 = rnd(1);
        for (j = 0; j < N; j++) {
            if (*(val+j) > (double)rnd0) {
                ns[i] = j;
                break;
            }
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < NBIT; j++) {
            chng = *chrm + NBIT * ns[i] + j;
            newch[i][j] = *chng;
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < NBIT; j++) {
            chng = *chrm + NBIT * i + j;
            *chng = newch[i][j];
        }
    }
    return;
}

void showChrom(int chrm[][NBIT]) {
    int i, j;
    for (i = 0; i < N; i++) {
        cout << i << ":\t";
        for (j = 0; j < NBIT; j++) cout << chrm[i][j];
        cout << endl;
    }
    cout << endl;
    return;
}

void showRes(int chrm[][NBIT], double x[N], double z[N]) {
    int i, j;
    for (i = 0; i < N; i++) {
        cout << i << ":\t";
        for (j = 0; j < NBIT; j++) cout << chrm[i][j];
        printf("  %10.6f  %10.6f", x[i], z[i]);
        cout << endl;
    }
    cout << endl;
    return;
}

void showPq(double *sp, double *cp) {
    int i;
    printf("          p_i         q_i\n");
    for (i = 0; i < N; i++) {
        cout << i << ":\t";
        printf("%10.6f  %10.6f", sp[i], cp[i]);
        cout << endl;
    }
    cout << endl;
    return;
}