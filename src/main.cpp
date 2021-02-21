#include "monitor.h"
#include "genetic.h"

//次回やること。
//daq.cppで計算した結果をもとにビームサイズのパラメータを求める。
//以前作った遺伝的アルゴリズムを改良し、作り直す。
//サイズ情報と信号情報を組み合わせて、ビームのRMS分布を再構成する。

void dataRead();
void doGenetic();

int main() {
    std::cout << "*******************************************************************************" << std::endl;
    doGenetic();
    std::cout << "*\tDone" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    return 0;
}

void doGenetic() {
    int i;
    Beam_t particle[20];
    std::cout << "*\tInitializing particle position" << std::endl;
    for (i = 0; i < 20; i++) {
        initPos(particle[i]);
    }
    return;
}

void dataRead() {
    monitor16pu beamMonitor;
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
    MatrixXd mom13 = beamMonitor.getMoment(13, beamMonitor.getVol(13, waveData13, 2, 0));
    MatrixXd mom15 = beamMonitor.getMoment(15, beamMonitor.getVol(15, waveData15, 2, 0));
    beamMonitor.calEmit(mom13, mom15);
    return;
}