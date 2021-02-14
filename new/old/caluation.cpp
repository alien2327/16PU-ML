#define _USE_MATH_DEFINES

#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/qvm/mat_operations.hpp>

#include "daq.cpp"

//次回やること。
//daq.cppで計算した結果をもとにビームサイズのパラメータを求める。
//以前作った遺伝的アルゴリズムを改良し、作り直す。
//サイズ情報と信号情報を組み合わせて、ビームのRMS分布を再構成する。

boost::numeric::ublas::matrix<double> getEmittance(boost::numeric::ublas::matrix<double>, boost::numeric::ublas::matrix<double>);

int main() {

    char waveData13[100] = "../../data/beam/wave/wave_2019_12_20_11_00_45_address13.dat";
    char waveData15[100] = "../../data/beam/wave/wave_2019_12_20_11_00_45_address15.dat";
    boost::numeric::ublas::matrix<double> volData13 = getVoltage(13, waveData13);
    boost::numeric::ublas::matrix<double> volData15 = getVoltage(15, waveData15);
    boost::numeric::ublas::matrix<double> res13 = getMoment(13, waveData13);
    boost::numeric::ublas::matrix<double> res15 = getMoment(15, waveData15);

    boost::numeric::ublas::matrix<double> emittance = getEmittance(res13, res15);

    return 0;
}

boost::numeric::ublas::matrix<double> getEmittance(boost::numeric::ublas::matrix<double> mom13, boost::numeric::ublas::matrix<double> mom15) {
    addr13 monitor13;
    addr15 monitor15;
    boost::numeric::ublas::matrix<double> e(2,1);
    boost::numeric::ublas::matrix<double> beta(2,2);
    beta(0,0) = monitor13.beta[0]; beta(0,1) = monitor13.beta[1];
    beta(1,0) = monitor15.beta[0]; beta(1,1) = monitor15.beta[1];

    boost::qvm::inverse(beta);
    std::cout << beta << std::endl;

    return e;
}