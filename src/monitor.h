#ifndef MONITOR_16PU
#define MONITOR_16PU

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#define _Channel    16
#define _bufLen     65536
#define _Turn       140
#define _Bunch      9
#define _Sampling   52
#define M_PI        3.14159265358979323846

using namespace Eigen;

class monitor16pu {
    private:
    MatrixXd gain13 = MatrixXd::Constant(_Channel, 1, 0);
    MatrixXd align13 = MatrixXd::Constant(2, 1, 0);
    MatrixXd beta13 = MatrixXd::Constant(2, 1, 0);
    MatrixXd matInvNorm13 = MatrixXd::Constant(_Channel, _Channel, 0);

    MatrixXd gain15 = MatrixXd::Constant(_Channel, 1, 0);
    MatrixXd align15 = MatrixXd::Constant(2, 1, 0);
    MatrixXd beta15 = MatrixXd::Constant(2, 1, 0);
    MatrixXd matInvNorm15 = MatrixXd::Constant(_Channel, _Channel, 0);

    public:
    monitor16pu();

    //DAQ
    MatrixXd getVol(int, char *, double, double);

    //Analysis
    MatrixXd getMoment(int, MatrixXd);
    void calEmit(MatrixXd, MatrixXd);
};

#endif