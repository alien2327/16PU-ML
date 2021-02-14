#include <iostream>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>

using namespace boost::qvm;

int main() {
    using mat = mat<double, 2, 2>;
    auto out = [](mat const& m) {
        std::cout << "(("
            << A00(m) << ", "
            << A01(m) << "), ("
            << A10(m) << ", "
            << A11(m) << "))" << std::endl;
    };

    mat a, b;
    A00(a) = 1;
    A01(a) = 2;
    A10(a) = 3;
    A11(a) = 4;

    out (a);

    b = 2 * a + a * 3 - a / 4;
    out(b);

    b += a;
    out(b);
    b -= a;
    out(b);
    b *= a;
    out(b);
    b /= 3;
    out(b);

    return 0;
}