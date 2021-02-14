#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

boost::numeric::ublas::matrix<int> CreateMatrix(const int n_cols, const int n_rows, int *, int);
boost::numeric::ublas::matrix<double> CreateMatrix(const int n_cols, const int n_rows, double *, int);

boost::numeric::ublas::matrix<int> CreateMatrix(const int n_cols, const int n_rows, int *arr, int size) {
    boost::numeric::ublas::matrix<int> m(n_cols, n_rows);
    if (size != n_cols * n_rows) {
        printf("Size error");
        return m;
    }
    for(size_t x = 0; x < n_cols; x++) {
        for(size_t y = 0; y < n_rows; y++) {
            m(x, y) = arr[x * n_cols + y];
        }
    }    
    return m;
}

boost::numeric::ublas::matrix<double> CreateMatrix(const int n_cols, const int n_rows, double *arr, int size) {
    boost::numeric::ublas::matrix<double> m(n_cols, n_rows);
    if (size != n_cols * n_rows) {
        printf("Size error");
        return m;
    }
    for(size_t x = 0; x < n_cols; x++) {
        for(size_t y = 0; y < n_rows; y++) {
            m(x, y) = arr[x * n_cols + y];
        }
    }    
    return m;
}