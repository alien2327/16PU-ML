#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

double mean(boost::numeric::ublas::vector<int>, int);
double mean(boost::numeric::ublas::vector<double>, int);

double variance(boost::numeric::ublas::vector<int>, int);
double variance(boost::numeric::ublas::vector<double>, int);

double standard_deviation(boost::numeric::ublas::vector<int>, int);
double standard_deviation(boost::numeric::ublas::vector<double>, int);

double skewness(boost::numeric::ublas::vector<int>, int);
double skewness(boost::numeric::ublas::vector<double>, int);

double kurtosis(boost::numeric::ublas::vector<int>, int);
double kurtosis(boost::numeric::ublas::vector<double>, int);

void goertzel_method(std::vector<boost::numeric::ublas::vector<int>> &, std::vector<boost::numeric::ublas::vector<double>> &, double, double);

double mean(boost::numeric::ublas::vector<int> _input, int size) {
    int i;
    double res = 0.0;
    for (i = 0; i < size; i++) res += (double)_input(i);
    return res/(double)size;
}
double mean(boost::numeric::ublas::vector<double> _input, int size) {
    int i;
    double res = 0.0;
    for (i = 0; i < size; i++) res += _input(i);
    return res/(double)size;
}

double variance(boost::numeric::ublas::vector<int> _input, int size) {
    int i;
    double res = 0.0, avg = mean(_input, size);
    for (i = 0; i < size; i++) {
        res += pow((_input(i) + avg), 2);
    }
    return res/(double)size;
}
double variance(boost::numeric::ublas::vector<double> _input, int size) {
    int i;
    double res = 0.0, avg = mean(_input, size);
    for (i = 0; i < size; i++) {
        res += pow((_input(i) + avg), 2);
    }
    return res/(double)size;
}

double standard_deviation(boost::numeric::ublas::vector<int> _input, int size) {
    double res = sqrt(variance(_input, size));
    return res;
}
double standard_deviation(boost::numeric::ublas::vector<double> _input, int size) {
    double res = sqrt(variance(_input, size));
    return res;
}

double skewness(boost::numeric::ublas::vector<int> _input, int size)
{
    int i;
    double res, avg = mean(_input, size), r = (double)size/(((double)size-1)*((double)size-2)), t;
    for (i = 0; i < size; i++) {
        t = pow((_input(i)-avg)/standard_deviation(_input, size), 3);
        res += t;
    }
    res *= r;
    return res;
}
double skewness(boost::numeric::ublas::vector<double> _input, int size)
{
    int i;
    double res, avg = mean(_input, size), r = (double)size/(((double)size-1)*((double)size-2)), t;
    for (i = 0; i < size; i++) {
        t = pow((_input(i)-avg)/standard_deviation(_input, size), 3);
        res += t;
    }
    res *= r;
    return res;
} 

double kurtosis(boost::numeric::ublas::vector<int> _input, int size)
{
    int i;
    double r1, r2, t, res, avg = mean(_input, size);
    r1 = ((double)size*((double)size+1))/(((double)size-1)*((double)size-2)*((double)size-3));
    r2 = (3*((double)size-1)*((double)size-1))/(((double)size-2)*((double)size-3));
    for (i = 0; i < size; i++) {
        t = pow(_input(i)-avg, 4)/pow(standard_deviation(_input, size), 4);
        res += t;
    }
    res = r1*res - r2;
    return res;
}
double kurtosis(boost::numeric::ublas::vector<double> _input, int size)
{
    int i;
    double r1, r2, t, res, avg = mean(_input, size);
    r1 = ((double)size*((double)size+1))/(((double)size-1)*((double)size-2)*((double)size-3));
    r2 = (3*((double)size-1)*((double)size-1))/(((double)size-2)*((double)size-3));
    for (i = 0; i < size; i++) {
        t = pow(_input(i)-avg, 4)/pow(standard_deviation(_input, size), 4);
        res += t;
    }
    res = r1*res - r2;
    return res;
}

