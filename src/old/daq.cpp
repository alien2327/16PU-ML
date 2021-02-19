#define _USE_MATH_DEFINES

#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "parameter.hpp"
#include "static.cpp"
#include "matOper.cpp"
#include "dsp.cpp"

void readMatrix();
void readData(int, char *, double *, boost::numeric::ublas::matrix<double> &);
void extractData(char *, boost::numeric::ublas::matrix<double>);
boost::numeric::ublas::matrix<double> getMoment(int, char *);
boost::numeric::ublas::matrix<double> getVoltage(int, char *);
/*
int main() {
    int i;
    addr13 monitor13;
    addr15 monitor15;
    char waveData13[100] = "../../data/beam/wave/wave_2019_12_20_11_00_45_address13.dat";
    char waveData15[100] = "../../data/beam/wave/wave_2019_12_20_11_00_45_address15.dat";
    char fName13[100] = "./moment13.dat";
    char fName15[100] = "./moment15.dat";
    char fTune13[100] = "./tune13.dat";
    char fTune15[100] = "./tune15.dat";

    boost::numeric::ublas::matrix<double> volData13(16, 140*9);
    boost::numeric::ublas::matrix<double> volData15(16, 140*9);

    readData(13, waveData13, monitor13.gain, volData13);
    readData(15, waveData15, monitor15.gain, volData15);

    boost::numeric::ublas::matrix<double> matInverse13 = CreateMatrix(16, 16, monitor13.mat, 16*16);
    boost::numeric::ublas::matrix<double> matInverse15 = CreateMatrix(16, 16, monitor15.mat, 16*16);

    boost::numeric::ublas::matrix<double> res13 = boost::numeric::ublas::prod(matInverse13, volData13);
    boost::numeric::ublas::matrix<double> res15 = boost::numeric::ublas::prod(matInverse15, volData15);
    boost::numeric::ublas::matrix<double> tuneData13(2, 140);
    boost::numeric::ublas::matrix<double> tuneData15(2, 140);

    std::vector<double> x13Real(140, 0); std::vector<double> x13Imag(140, 0);
    std::vector<double> y13Real(140, 0); std::vector<double> y13Imag(140, 0);

    std::vector<double> x15Real(140, 0); std::vector<double> x15Imag(140, 0);
    std::vector<double> y15Real(140, 0); std::vector<double> y15Imag(140, 0);

    for (i = 0; i < 140; i++) {
        res13(1, i*9) = res13(1, i*9) - monitor13.align[0]; x13Real[i] = res13(1, i*9);
        res13(2, i*9) = res13(2, i*9) - monitor13.align[1]; y13Real[i] = res13(2, i*9);
        res15(1, i*9) = res15(1, i*9) - monitor15.align[0]; x15Real[i] = res15(1, i*9);
        res15(2, i*9) = res15(2, i*9) - monitor15.align[1]; y15Real[i] = res15(2, i*9);
    }

    transform(x13Real, x13Imag); transform(y13Real, y13Imag);
    transform(x15Real, x15Imag); transform(y15Real, y15Imag);

    for (i = 0; i < 140; i++) {
        tuneData13(0, i) = sqrt(pow(x13Real[i],2)+pow(x13Imag[i],2));
        tuneData13(1, i) = sqrt(pow(y13Real[i],2)+pow(y13Imag[i],2));
        tuneData15(0, i) = sqrt(pow(x15Real[i],2)+pow(x15Imag[i],2));
        tuneData15(1, i) = sqrt(pow(y15Real[i],2)+pow(y15Imag[i],2));
    }

    std::cout << "Exctracting moment data" << std::endl;
    extractData(fName13, res13);
    extractData(fName15, res15);
    std::cout << "Moment data successfully exctrated" << std::endl;
    std::cout << "Exctracting tune data" << std::endl;
    extractData(fTune13, tuneData13);
    extractData(fTune15, tuneData15);
    std::cout << "Tune data successfully exctrated" << std::endl;

    return 0;
}
*/

boost::numeric::ublas::matrix<double> getMoment(int addr, char fName[]) {
    int i, j;
    if (addr == 13) {
        addr13 monitor13;
        boost::numeric::ublas::matrix<double> volData13(16, 140*9);
        readData(13, fName, monitor13.gain, volData13);
        boost::numeric::ublas::matrix<double> matInverse13 = CreateMatrix(16, 16, monitor13.mat, 16*16);
        boost::numeric::ublas::matrix<double> res13 = boost::numeric::ublas::prod(matInverse13, volData13);
        for (i = 0; i < res13.size2(); i++) {
            double norm = res13(0, i);
            for (j = 0; j < res13.size1(); j++) {
                res13(j, i) /= norm;
            }
        }
        return res13;
    } else if (addr == 15) {
        addr15 monitor15;
        boost::numeric::ublas::matrix<double> volData15(16, 140*9);
        readData(15, fName, monitor15.gain, volData15);
        boost::numeric::ublas::matrix<double> matInverse15 = CreateMatrix(16, 16, monitor15.mat, 16*16);
        boost::numeric::ublas::matrix<double> res15 = boost::numeric::ublas::prod(matInverse15, volData15);
        for (i = 0; i < res15.size2(); i++) {
            double norm = res15(0, i);
            for (j = 0; j < res15.size1(); j++) {
                res15(j, i) /= norm;
            }
        }
        return res15;
    } else {
        boost::numeric::ublas::matrix<double> error(1, 1);
        return error;
    }
}

boost::numeric::ublas::matrix<double> getVoltage(int addr, char fName[]) {
    if (addr == 13) {
        addr13 monitor13;
        boost::numeric::ublas::matrix<double> volData13(16, 140*9);
        readData(13, fName, monitor13.gain, volData13);
        return volData13;
    } else if (addr == 15) {
        addr15 monitor15;
        boost::numeric::ublas::matrix<double> volData15(16, 140*9);
        readData(15, fName, monitor15.gain, volData15);
        return volData15;
    } else {
        boost::numeric::ublas::matrix<double> error(1, 1);
        return error;
    }
}

void extractData(char _fName[], boost::numeric::ublas::matrix<double> input) {
    size_t i, j;
    std::ofstream outfile;
    outfile.open(_fName);
    if (outfile.is_open()) {
        if (input.size2() > 140) {
            for (i = 0; i < input.size2()/9; i++) {
                outfile << std::to_string(i) + " ";
                for (j = 0; j < input.size1(); j++) {
                    outfile <<  std::to_string(input(j, i*9)) + " ";
                }
                outfile << "\n";
            }
        } else {
            for (i = 0; i < input.size2(); i++) {
                outfile << std::to_string((double)i/input.size2()) + " ";
                for (j = 0; j < input.size1(); j++) {
                    outfile <<  std::to_string(input(j, i)) + " ";
                }
                outfile << "\n";
            }
        }
        outfile.close();
    }
}

void readData(int addr, char _input[], double *gain, boost::numeric::ublas::matrix<double> &_output) {
    int i, j, k, size = 468, samplelength = 0;
    boost::numeric::ublas::matrix<int> bufData(16, 65536);
    std::ifstream fin;
    fin.open(_input, std::ios::in);
    if (fin.is_open()) {
        for (i = 0; i < 16; i++) {
            boost::numeric::ublas::vector<int> vol(65536);
            char buf[131056];
            fin.seekg(25 + i*131066);
            fin.read(buf, 131056);
            for (j = 0; j < 65536; j++) {
                if (addr == 13) {
                    vol(j) = gain[i] * ((unsigned char) buf[0 + j*2])*64 + ((unsigned char) buf[1 + j*2])/4;
                } else if (addr == 15) {
                    vol(j) = gain[i] * ((unsigned char) buf[0 + j*2])*64 + ((unsigned char) buf[1 + j*2])/4;
                }
            }
            int lag = 0, threshold = 0;
            double avg = mean(vol, 65536);
            threshold = avg + 50;
            for (j = 0; j < size; j++) {
                if ((vol(j) + vol(j+1) + vol(j+2))/3 > threshold) {
                    lag += j - 4;
                    break;
                }
            }
            for (j = lag; j < 65536; j++) bufData(i, j-lag) = vol(j);
        }
        if (addr == 13) {
            int vol[65536];
            for (i = 0; i < 65536; i++) {
                vol[i] = bufData(10, i);
                bufData(10, i) = bufData(11, i);
                bufData(11, i) = vol[i];
            }
        }
    } else {
        puts("File not opened");
    }
    fin.close();
    goertzel_method(bufData, _output, 1, 0);
    return;
}
