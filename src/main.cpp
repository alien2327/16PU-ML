#include "monitor.h"
#include "genetic.h"

#include <string>
#include <vector>
#include <algorithm>

//次回やること。
//daq.cppで計算した結果をもとにビームサイズのパラメータを求める。
//以前作った遺伝的アルゴリズムを改良し、作り直す。
//サイズ情報と信号情報を組み合わせて、ビームのRMS分布を再構成する。

//次回やること。
//バンチを変える演算子の導入 v
//従来の遺伝的アルゴリズムにするか、論文に従うか
//評価関数・損失関数の導入 v
//報告書・発表資料の作成

void dataRead(MatrixXd &, MatrixXd &, MatrixXd &);
void displayStat(Beam_t *, MatrixXd);
void doGenetic(MatrixXd, MatrixXd);
void testSet(MatrixXd &, MatrixXd &);
void printOut(Beam_t);

double mrate = 0.005;

monitor16pu beamMonitor;

int main() {
    init_genrand((unsigned)time(NULL));
    std::cout << "*****************************************************************************************" << std::endl;
    std::cout << "*" << std::endl;
    //MatrixXd vol13 = MatrixXd::Constant(_Channel, _Turn*_Bunch, 0);
    //MatrixXd vol15 = MatrixXd::Constant(_Channel, _Turn*_Bunch, 0);
    //MatrixXd beamsize = MatrixXd::Constant(4, _Turn*_Bunch, 0);
    //dataRead(vol13, vol15, beamsize);

    MatrixXd vol13 = MatrixXd::Constant(_Channel, 1, 0);
    MatrixXd beamsize = MatrixXd::Constant(4, 1, 0);

    testSet(vol13, beamsize);

    doGenetic(vol13, beamsize);
    std::cout << "*\tDone" << std::endl;
    std::cout << "*" << std::endl;
    std::cout << "*****************************************************************************************" << std::endl;
    return 0;
}

void testSet(MatrixXd &vol, MatrixXd &beamsize) {
    tBeam_t testparticle;
    testPart(testparticle);
    vol = func(testparticle, beamMonitor.mat13);
    beamsize(0,0) = meansqure(testparticle, 0);
    beamsize(1,0) = meansqure(testparticle, 1);
    printf("*\n*\tTest Beam Position <x>: %1.5f\n", mean(testparticle, 0));
    printf("*\tTest Beam Position <y>: %1.5f\n", mean(testparticle, 1));
    printf("*\tTest Beam Size sigma_x: %1.5f\n", beamsize(0, 0));
    printf("*\tTest Beam Size sigma_y: %1.5f\n*\n", beamsize(1, 0));
}

void doGenetic(MatrixXd v, MatrixXd bs) {
    int i, t = 0;
    Beam_t particle[50];
    Beam_t particlesort[50];
    MatrixXd::Index row, col;
    Vol_t vol[50];
    //MatrixXd loss = MatrixXd::Constant(50, 1, 0);
    MatrixXd beamsize = MatrixXd::Constant(2, 1, 0);
    std::cout << "*\tInitializing particle position" << std::endl;
    for (i = 0; i < 50; i++) {
        initPos(particle[i]);
        vol[i] = func(particle[i], beamMonitor.mat13);
        beamsize(0,0) = meansqure(particle[i], 0);
        beamsize(1,0) = meansqure(particle[i], 1);
    }
    //displayStat(particle, loss);

    do {
        std::vector<std::pair<int, double>> loss;
        for (i = 0; i < 50; i++) {
            applyOper(particle[i], 16 * genrand_real1(), 0.01);
            vol[i] = func(particle[i], beamMonitor.mat13);
            beamsize(0,0) = meansqure(particle[i], 0);
            beamsize(1,0) = meansqure(particle[i], 1);
            loss.push_back(std::make_pair(i, mse(v.block(0,0,16,1), vol[i]) + ep(bs.block(0,0,3,1), beamsize)));
            particlesort[i] = particle[i];
        }
        //displayStat(particle, loss);
        //
        //particlesort = particle;
        sort(loss.begin(), loss.end());
        printf("\r*\tMinimum Loss: %1.10e", loss[0].second);
        for (i = 0; i < 50; i++) {
            particle[i] = particlesort[loss[(int)i/10].first];
            //loss[i] = mse(v.block(0,0,16,1), vol[i]) + ep(bs.block(0,0,3,1), beamsize);
            applyMut(particle[i], mrate);
        }
        t++;
    }while(t < 1000);
    std::cout << std::endl;

    printf("*\n*\tBeam Position <x>: %1.5f\n", mean(particle[0], 0));
    printf("*\tBeam Position <y>: %1.5f\n", mean(particle[0], 1));
    printf("*\tBeam Size sigma_x: %1.5f\n", meansqure(particle[0], 0));
    printf("*\tBeam Size sigma_y: %1.5f\n*\n", meansqure(particle[0], 1));
    printOut(particle[0]);

    //displayStat(particle, loss);

    return;
}

void dataRead(MatrixXd &vol13, MatrixXd &vol15, MatrixXd &beamsize) {
    std::string ymd, hms;
    std::string fName13, fName15;
    std::cout << "*\tInput Year/Month/Day: ";
    std::cin >> ymd;
    std::cout << "*\tInput Hour/Minute/Second: ";
    std::cin >> hms;
    fName13 = "../16PU/data/beam/wave/wave_" + ymd.substr(0, 4) + "_" + ymd.substr(4, 2) + "_" + ymd.substr(6, 2) + "_" + hms.substr(0, 2) + "_" + hms.substr(2, 2) + "_" + hms.substr(4, 2) + "_address13.dat";
    fName15 = "../16PU/data/beam/wave/wave_" + ymd.substr(0, 4) + "_" + ymd.substr(4, 2) + "_" + ymd.substr(6, 2) + "_" + hms.substr(0, 2) + "_" + hms.substr(2, 2) + "_" + hms.substr(4, 2) + "_address15.dat";
    const char * waveData13 = fName13.c_str();
    const char * waveData15 = fName15.c_str();
    vol13 = beamMonitor.getVol(13, waveData13, 2, 0);
    vol15 = beamMonitor.getVol(15, waveData15, 2, 0);
    MatrixXd mom13 = beamMonitor.getMoment(13, vol13);
    MatrixXd mom15 = beamMonitor.getMoment(15, vol15);
    beamsize = beamMonitor.calEmit(mom13, mom15);
    printf("*\n*\tBeam Position <x>: %1.5f\n", mom13(1, 0));
    printf("*\tBeam Position <y>: %1.5f\n", mom13(2, 0));
    printf("*\tBeam Size sigma_x: %1.5f\n", beamsize(0, 0));
    printf("*\tBeam Size sigma_y: %1.5f\n*\n", beamsize(1, 0));
    return;
}

void displayStat(Beam_t *particle, MatrixXd loss) {
    int i, j;
    double mx, my, sx, sy;
    std::cout << "*\t" << std::endl;
    for (i = 0; i < 4; i++) {
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            std::cout << 4*i+j+1 << "-st bunch\t\t";
        }
        std::cout << "" << std::endl;
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            mx = mean(particle[4*i+j], 0);
            printf("mean x : %1.5e\t", mx);
        }
        std::cout << std::endl;
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            my = mean(particle[4*i+j], 1);
            printf("mean y : %1.5e\t", my);
        }
        std::cout << std::endl;
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            sx = standard_deviation(particle[4*i+j], 0);
            printf("sigma x: %1.5e\t", sx);
        }
        std::cout << std::endl;
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            sy = standard_deviation(particle[4*i+j], 1);
            printf("sigma y: %1.5e\t", sy);
        }
        std::cout << std::endl;
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            printf("loss   : %1.5e\t", loss(4*i+j, 0));
        }
        std::cout << std::endl;
        std::cout << "*" << std::endl;
    }

    return;
}

void printOut(Beam_t particle) {
    int i;
    std::string pos;
    std::ofstream posFile;
    posFile.open("./result/final.dat");
    for (i = 0; i < particle.cols(); i++) {
        pos = std::to_string(particle(0, i)) + "," + std::to_string(particle(1, i)) + "\n";
        posFile << pos.c_str();
    }
    posFile << "-82.5, -82.5\n";
    posFile << "-82.5, 82.5\n";
    posFile << "82.5, -82.5\n";
    posFile << "82.5, 82.5\n";
}