#include "monitor.h"

//次回やること。
//daq.cppで計算した結果をもとにビームサイズのパラメータを求める。
//以前作った遺伝的アルゴリズムを改良し、作り直す。
//サイズ情報と信号情報を組み合わせて、ビームのRMS分布を再構成する。

int main() {

    std::cout << "*******************************************************************************" << std::endl;

    monitor16pu beamMonitor;

    char waveData13[100] = "../16PU/data/beam/wave/wave_2019_12_20_11_00_45_address13.dat";
    char waveData15[100] = "../16PU/data/beam/wave/wave_2019_12_20_11_00_45_address15.dat";

    MatrixXd mom13 = beamMonitor.getMoment(13, beamMonitor.getVol(13, waveData13, 2, 0));
    MatrixXd mom15 = beamMonitor.getMoment(15, beamMonitor.getVol(15, waveData15, 2, 0));

    beamMonitor.calEmit(mom13, mom15);

    std::cout << "*\tDone" << std::endl;

    std::cout << "*******************************************************************************" << std::endl;
    
    return 0;
}