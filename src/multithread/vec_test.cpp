#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <ctime>
#include <complex>
#include <cmath>

using namespace std;

int main()
{
    vector<int> a(2, 0);
    vector<int> b(2, 0);
    a[0] = 23;
    a[1] = 11;
    b = a;
    printf("%d, %d, %d, %d", a[0], a[1], b[0], b[1]);
    return 0;
}