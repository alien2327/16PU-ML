#include "monitor.h"
#include "genetic.h"

//次回やること。
//daq.cppで計算した結果をもとにビームサイズのパラメータを求める。
//以前作った遺伝的アルゴリズムを改良し、作り直す。
//サイズ情報と信号情報を組み合わせて、ビームのRMS分布を再構成する。

//次回やること。
//バンチを変える演算子の導入
//従来の遺伝的アルゴリズムにするか、論文に従うか
//評価関数・損失関数の導入
//報告書・発表資料の作成

void dataRead(MatrixXd &, MatrixXd &, MatrixXd &);
void displayStat(Beam_t *, MatrixXd);
void doGenetic(MatrixXd, MatrixXd);

monitor16pu beamMonitor;

int main() {
    init_genrand((unsigned)time(NULL));
    std::cout << "*****************************************************************************************" << std::endl;
    std::cout << "*" << std::endl;
    MatrixXd vol13 = MatrixXd::Constant(_Channel, _Turn*_Bunch, 0);
    MatrixXd vol15 = MatrixXd::Constant(_Channel, _Turn*_Bunch, 0);
    MatrixXd beamsize = MatrixXd::Constant(4, _Turn*_Bunch, 0);
    dataRead(vol13, vol15, beamsize);
    doGenetic(vol13, beamsize);
    std::cout << "*\tDone" << std::endl;
    std::cout << "*" << std::endl;
    std::cout << "*****************************************************************************************" << std::endl;
    return 0;
}

void doGenetic(MatrixXd v, MatrixXd bs) {
    int i;
    Beam_t particle[20];
    Vol_t vol[20];
    MatrixXd loss = MatrixXd::Constant(20, 1, 0);
    MatrixXd beamsize = MatrixXd::Constant(2, 1, 0);
    std::cout << "*\tInitializing particle position" << std::endl;
    for (i = 0; i < 20; i++) {
        initPos(particle[i]);
        vol[i] = func(particle[i], beamMonitor.mat13);
        beamsize(0,0) = standard_deviation(particle[i], 0);
        beamsize(1,0) = standard_deviation(particle[i], 1);
        loss(i, 0) = mse(v.block(0,0,16,1), vol[i]) + ep(bs.block(0,0,3,1), beamsize);
    }
    displayStat(particle, loss);
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
    return;
}

void displayStat(Beam_t *particle, MatrixXd loss) {
    int i, j;
    double mx, my, sx, sy;
    std::cout << "*\t" << std::endl;
    for (i = 0; i < 5; i++) {
        std::cout << "*\t";
        for (j = 0; j < 4; j++) {
            std::cout << 4*i+j << "-st bunch\t\t";
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