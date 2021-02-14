#include <iostream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;

int main() {
    accumulator_set<double, stats<tag::min, tag::mean, tag::sum> > acc;

    acc(3.0);
    acc(1.0);
    acc(4.0);
    acc(2.0);
    acc(5.0);

    std::cout << extract::min(acc) << std::endl;  // 最小値
    std::cout << extract::mean(acc) << std::endl; // 平均値
    std::cout << extract::sum(acc) << std::endl;  // 合計値
}