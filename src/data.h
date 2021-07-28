#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <complex>
#include <string>
#include <bitset>
#include <cmath>

using namespace std;

class Data {
private:
	int samplelength;
	const double gain13[16] = { 1.000000000000000000e0, 9.864749932865497506e-1, 1.008554910983899822e0, 9.909215189952645941e-1, 1.016551297946559407e0, 1.004441025921927633e0, 1.023606422028664920e0, 1.007893430992459427e0,
	1.019020740443200568e0, 1.002733038130677601e0, 1.021101905384611097e0, 1.002615873152301385e0, 1.045142626872535008e0, 9.986874348446839189e-1, 1.022369446580688512e0, 9.918986411869963327e-1 };
	const double gain15[16] = { 1.000000000000000000e0, 9.927984352486749486e-1, 9.886965960947093901e-1, 9.844666881364172450e-1, 9.811551707487992102e-1, 9.783521372729635512e-1, 9.896903610596464729e-1, 9.869302062198921366e-1,
	9.906593504726489696e-1, 9.859330582601781856e-1, 9.980913383872332956e-1, 9.882277496998294053e-1, 1.000800169181874200e0, 9.924633771768432977e-1, 9.944628454760420233e-1, 1.001591146042800640e0 };


public:
	void read_data(std::string, std::string, int**, int**);
	void extract_momentum(double**, double**, double*, double*, int, int, char);
	void dot_matrix(double**, double**, char, char);
	void goertzel_method(int**, double**, double, double);

	double mean(double*, int);
	double variance(double*, int);
	double covariance(double*, double*, int);
	double standard_deviation(double*, int);
	double skewness(double*, int);
	double kurtosis(double*, int);
};