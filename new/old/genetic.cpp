#define testParticle 100
#define candParticle 20
#define mRate 0.01
#define cRate 0.15
#define _USE_MATH_DEFINES

#include <iostream>
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

struct Chromo {
    bool x[8];
    bool y[8];
};

void genChrome(int, Chromo &);

int main() {
    int i, j;
    boost::numeric::ublas::matrix<Chromo> p(candParticle,testParticle);
    for (i = 0; i < candParticle; i++) {
        for (j = 0; j < testParticle; j++) {
            genChrome(i*j, p(i, j));
        }
    }
    return 0;
}

void genChrome(int seed, Chromo &p) {
    int i;
    boost::random::mt19937 gen(seed);
    boost::random::uniform_int_distribution<> dist(0, 1);
    for (i = 0; i < 8; i++) {
        p.x[i] = dist(gen); p.y[i] = dist(gen);
    }
    return;
}