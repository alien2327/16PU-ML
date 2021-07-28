#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <conio.h>
#include <stdlib.h>
#include <filesystem>
#include <thread>
#include <random>
#include <vector>
#include <array>
#include <cmath>

#include "param.h"
#include "data.h"

std::uniform_real_distribution<> d_ds(-0.01, 0.01);
std::uniform_real_distribution<> d_dtheta(-1.5, 1.5);
std::uniform_real_distribution<> d_dsteer(-0.0005, 0.0005);

std::uniform_real_distribution<> meanForTest(-10, 10);
std::uniform_real_distribution<> stdForTest(1, 10);
std::uniform_real_distribution<> meanForTrue(-1, 1);
std::uniform_real_distribution<> stdForTrue(1, 10);

std::uniform_int_distribution<> mode_selecting(0, 6);

std::uniform_real_distribution<> percentage(0.0, 1.0);

class Space {
private:
	double multipole[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	double voltage[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	double* grid;

	Mon13 mon;

public:
	double bunch[BUNCH_SIZE][2];
	std::string parentDirName = "./result/";

	Space(void);
	Space(double, double);
	Space(double, double, double, double);
	Space(std::string);

	void transform(int, double);
	void setMoment(int);
	void setMoment(int, double, double);
	void setVoltage(void);
	void calFDTD(void);
	void calFDTD(double, double);
	void saveBunch(int);

	double* at(int, int);
	double* intersection(double, double);
	double* getCoord(int);
	double* getMoment(void);
	double* getVoltage(void);
};


Space::Space(void)
{
	double mx = meanForTest(engine);
	double my = meanForTest(engine);
	double sx = stdForTest(engine);
	double sy = stdForTest(engine);
	std::normal_distribution<> dx(mx, sx);
	std::normal_distribution<> dy(my, sy);
	for (size_t i = 0; i < BUNCH_SIZE; ++i) {
		*(*(bunch + i)) = dx(engine);
		*(*(bunch + i) + 1) = dy(engine);
	}
}


Space::Space(double mx, double my)
{
	double sx = stdForTest(engine);
	double sy = stdForTest(engine);
	std::normal_distribution<> dx(mx, sx);
	std::normal_distribution<> dy(my, sy);
	for (size_t i = 0; i < BUNCH_SIZE; ++i) {
		*(*(bunch + i)) = dx(engine);
		*(*(bunch + i) + 1) = dy(engine);
	}
}


Space::Space(double mx, double sx, double my, double sy)
{
	std::normal_distribution<> dx(mx, sx);
	std::normal_distribution<> dy(my, sy);
	for (size_t i = 0; i < BUNCH_SIZE; ++i) {
		*(*(bunch + i)) = dx(engine);
		*(*(bunch + i) + 1) = dy(engine);
	}
}

Space::Space(std::string fname)
{
	std::ifstream is(fname, std::ifstream::binary);
	if (is) {
		is.seekg(0);
		for (size_t i = 0; i < BUNCH_SIZE; ++i) is.read((char*)&bunch[i][0], sizeof(double));
		for (size_t i = 0; i < BUNCH_SIZE; ++i) is.read((char*)&bunch[i][1], sizeof(double));
	}
	else {
		double mx = meanForTest(engine);
		double my = meanForTest(engine);
		double sx = stdForTest(engine);
		double sy = stdForTest(engine);
		std::normal_distribution<> dx(mx, sx);
		std::normal_distribution<> dy(my, sy);
		for (size_t i = 0; i < BUNCH_SIZE; ++i) {
			*(*(bunch + i)) = dx(engine);
			*(*(bunch + i) + 1) = dy(engine);
		}
	}
}

void Space::transform(int mode, double d) {
	int j;

	double dx = 100 * d_ds(engine);
	double dy = 100 * d_ds(engine);
	double dsx = d * d_ds(engine);
	double dsy = d * d_ds(engine);
	double ddx = d * d_ds(engine);
	double ddy = d * d_ds(engine);
	double dddx = d * d_ds(engine);
	double dddy = d * d_ds(engine);
	double theta = d_dtheta(engine);

	std::normal_distribution<> dx_r(meanForTest(engine), stdForTest(engine));
	std::normal_distribution<> dy_r(meanForTest(engine), stdForTest(engine));

	switch (mode)
	{
	case MODE_TRANS:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)); double y = *(*(bunch + j) + 1);
			*(*(bunch + j)) += dx + d_dsteer(engine);
			*(*(bunch + j) + 1) += dy + d_dsteer(engine);
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_SCALE:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)); double y = *(*(bunch + j) + 1);
			*(*(bunch + j)) += x * dsx + d_dsteer(engine);
			*(*(bunch + j) + 1) += y * dsy + d_dsteer(engine);
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_SHEER_1:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)); double y = *(*(bunch + j) + 1);
			*(*(bunch + j)) += dx * y + d_dsteer(engine);
			*(*(bunch + j) + 1) += dy * x + d_dsteer(engine);
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_SHEER_2:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)) + d_dsteer(engine); double y = *(*(bunch + j) + 1) + d_dsteer(engine);
			*(*(bunch + j)) = std::cos(theta) * x - std::sin(theta) * y;
			*(*(bunch + j) + 1) = std::sin(theta) * x + std::cos(theta) * y;
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_BEND_1:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)); double y = *(*(bunch + j) + 1);
			*(*(bunch + j)) += dx * std::abs(y) + d_dsteer(engine);
			*(*(bunch + j) + 1) += dy * std::abs(x) + d_dsteer(engine);
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_BEND_2:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			double x = *(*(bunch + j)); double y = *(*(bunch + j) + 1);
			*(*(bunch + j)) += ddx * std::pow(y, 2) + d_dsteer(engine);
			*(*(bunch + j) + 1) += ddy * std::pow(x, 2) + d_dsteer(engine);
			if (std::pow(*(*(bunch + j)), 2) + std::pow(*(*(bunch + j) + 1), 2) > std::pow(RADIUS - 5, 2)) {
				double* temp = intersection(*(*(bunch + j)), *(*(bunch + j) + 1));
				*(*(bunch + j)) = *temp; *(*(bunch + j) + 1) = *(temp + 1);
			}
		}
		break;
	case MODE_RAND:
		for (j = 0; j < BUNCH_SIZE; ++j) {
			*(*(bunch + j)) = dx_r(engine) + d_dsteer(engine);
			*(*(bunch + j) + 1) = dy_r(engine) + d_dsteer(engine);
		}
		break;
	default:
		break;
	}
}

void Space::setMoment(int mode)
{
	int i, j, k;
	int mono;
	switch (mode)
	{
	case FDTD_MODE:
		calFDTD();
		for (i = 0; i < GRID_NUMBER; ++i) {
			for (j = 0; j < GRID_NUMBER; ++j) {
				std::complex<double> pos(at(i, j)[0], at(i, j)[1]);
				int point = GRID_NUMBER * i + j;
				if (*(grid + point) == 0) continue;
				for (k = 0; k < 16; ++k) {
					if (k == 0) *(multipole + k) += *(grid + point);
					else if (k % 2 == 0) *(multipole + k) += *(grid + point) * std::imag(pow(pos, (int)(k + 1) / 2));
					else *(multipole + k) += *(grid + point) * std::real(pow(pos, (int)(k + 1) / 2));
				}
			}
		}
		mono = multipole[0];
		for (i = 0; i < 16; ++i) *(multipole + i) /= mono;
		delete[]grid;
		break;
	case DIRECT_MODE:
		for (i = 0; i < BUNCH_SIZE; ++i) {
			std::complex<double> pos(*(*(bunch + i)), *(*(bunch + i) + 1));
			for (j = 0; j < 16; j++) {
				if (j == 0) *(multipole + j) += 1;
				else if (j % 2 == 0) *(multipole + j) += std::imag(pow(pos, (int)(j + 1) / 2));
				else *(multipole + j) += std::real(pow(pos, (int)(j + 1) / 2));
			}
		}
		mono = *multipole;
		for (i = 0; i < 16; ++i) *(multipole + i) /= mono;
		break;
	default:
		break;
	}
}

void Space::setMoment(int mode, double offX, double offY)
{
	int i, j, k;
	int mono;
	switch (mode)
	{
	case FDTD_MODE:
		calFDTD();
		for (i = 0; i < GRID_NUMBER; ++i) {
			for (j = 0; j < GRID_NUMBER; ++j) {
				std::complex<double> pos(at(i, j)[0] + offX, at(i, j)[1] + offY);
				int point = GRID_NUMBER * i + j;
				if (*(grid + point) == 0) continue;
				for (k = 0; k < 16; ++k) {
					if (k == 0) *(multipole + k) += *(grid + point);
					else if (k % 2 == 0) *(multipole + k) += *(grid + point) * std::imag(pow(pos, (int)(k + 1) / 2));
					else *(multipole + k) += *(grid + point) * std::real(pow(pos, (int)(k + 1) / 2));
				}
			}
		}
		mono = multipole[0];
		for (i = 0; i < 16; ++i) *(multipole + i) /= mono;
		delete[]grid;
		break;
	case DIRECT_MODE:
		for (i = 0; i < BUNCH_SIZE; ++i) {
			std::complex<double> pos(*(*(bunch + i)) + offX, *(*(bunch + i) + 1) + offY);
			for (j = 0; j < 16; j++) {
				if (j == 0) *(multipole + j) += 1;
				else if (j % 2 == 0) *(multipole + j) += std::imag(pow(pos, (int)(j + 1) / 2));
				else *(multipole + j) += std::real(pow(pos, (int)(j + 1) / 2));
			}
		}
		mono = *multipole;
		for (i = 0; i < 16; ++i) *(multipole + i) /= mono;
		break;
	default:
		break;
	}

}


void Space::setVoltage(void)
{
	size_t i, k;
	for (i = 0; i < 16; i++) {
		for (k = 0; k < 16; k++) {
			*(voltage + i) += *(mon.param + (i * 16 + k)) * *(multipole + k);
		}
	}
	double norm = 0;
	for (i = 0; i < 16; ++i) norm += std::pow(*(voltage + i), 2);
	norm = std::sqrt(norm);
	for (i = 0; i < 16; ++i) *(voltage + i) /= norm;
}

void Space::calFDTD(void)
{
	size_t num = BUNCH_SIZE;
	int tot = GRID_NUMBER * GRID_NUMBER;
	double L = RADIUS*2;
	double dL = L / GRID_NUMBER;
	grid = new double[tot];
	for (size_t i = 0; i < tot; ++i) *(grid + i) = 0;
	for (size_t i = 0; i < num; ++i) {
		double x0 = *(*(bunch + i)), y0 = *(*(bunch + i) + 1);
		if (std::pow(x0, 2) + std::pow(y0, 2) > std::pow(RADIUS - 5, 2)) {
			double* inter = intersection(x0, y0);
			x0 = *inter; y0 = *(inter + 1);
		}
		int x, y, n1, n2, n3, n4;
		y = GRID_NUMBER - 1 - (y0 + (L - dL) / 2) / dL;
		x = (x0 + (L - dL) / 2) / dL;
		n1 = GRID_NUMBER * y + x;
		n2 = GRID_NUMBER * y + x + 1;
		n3 = GRID_NUMBER * (y - 1) + x;
		n4 = GRID_NUMBER * (y - 1) + x + 1;
		*(grid + n1) += (abs(at(y - 1, x + 1)[0] - x0) * abs(at(y - 1, x + 1)[1] - y0)) / (dL * dL);
		*(grid + n2) += (abs(at(y - 1, x)[0] - x0) * abs(at(y - 1, x)[1] - y0)) / (dL * dL);
		*(grid + n3) += (abs(at(y, x + 1)[0] - x0) * abs(at(y, x + 1)[1] - y0)) / (dL * dL);
		*(grid + n4) += (abs(at(y, x)[0] - x0) * abs(at(y, x)[1] - y0)) / (dL * dL);
	}
}


void Space::calFDTD(double offX, double offY)
{
	size_t num = BUNCH_SIZE;
	int tot = GRID_NUMBER * GRID_NUMBER;
	double L = RADIUS * 2;
	double dL = L / GRID_NUMBER;
	grid = new double[tot];
	for (size_t i = 0; i < tot; ++i) *(grid + i) = 0;
	for (size_t i = 0; i < num; ++i) {
		double x0 = *(*(bunch + i))-offX, y0 = *(*(bunch + i) + 1)-offY;
		if (std::pow(x0, 2) + std::pow(y0, 2) > std::pow(RADIUS - 5, 2)) {
			double* inter = intersection(x0, y0);
			x0 = *inter; y0 = *(inter + 1);
		}
		int x, y, n1, n2, n3, n4;
		y = GRID_NUMBER - 1 - (y0 + (L - dL) / 2) / dL;
		x = (x0 + (L - dL) / 2) / dL;
		n1 = GRID_NUMBER * y + x;
		n2 = GRID_NUMBER * y + x + 1;
		n3 = GRID_NUMBER * (y - 1) + x;
		n4 = GRID_NUMBER * (y - 1) + x + 1;
		*(grid + n1) += (abs(at(y - 1, x + 1)[0] - x0) * abs(at(y - 1, x + 1)[1] - y0)) / (dL * dL);
		*(grid + n2) += (abs(at(y - 1, x)[0] - x0) * abs(at(y - 1, x)[1] - y0)) / (dL * dL);
		*(grid + n3) += (abs(at(y, x + 1)[0] - x0) * abs(at(y, x + 1)[1] - y0)) / (dL * dL);
		*(grid + n4) += (abs(at(y, x)[0] - x0) * abs(at(y, x)[1] - y0)) / (dL * dL);
	}
}

void Space::saveBunch(int num)
{
	std::string chrName = "bunch_" + std::to_string(num) + ".dat";
	std::ofstream fout;
	fout.open(parentDirName + chrName, std::ios::out | std::ios::binary | std::ios::trunc);
	for (size_t i = 0; i < BUNCH_SIZE; ++i) fout.write((char*)&bunch[i][0], sizeof(double));
	for (size_t i = 0; i < BUNCH_SIZE; ++i) fout.write((char*)&bunch[i][1], sizeof(double));
	fout.close();
}

double* Space::at(int y, int x)
{
	double res[2];
	double L = RADIUS * 2;
	double dL = L / GRID_NUMBER;
	res[0] = (-L / 2 + dL / 2) + dL * y;
	res[1] = (L / 2 - dL / 2) - dL * x;
	return res;
}

double* Space::intersection(double x0, double y0)
{
	double res[2];
	double m, n;
	double A, B, C, D;
	double x_r, y_r, r = RADIUS - 5;
	x_r = 0;
	y_r = 0;
	if (x_r != x0) {
		m = (y_r - y0) / (x_r - x0);
		n = (y0 * x_r - x0 * y_r) / (x_r - x0);
		A = m * m + 1;
		B = (m * n - m * y_r - x_r);
		C = (x_r * x_r + y_r * y_r - r * r + n * n - 2 * n * y_r);
		D = B * B - A * C;
		if (x0 > x_r) res[0] = -(B - std::sqrt(D)) / A;
		else res[0] = -(B + std::sqrt(D)) / A;
		res[1] = m * res[0] + n;
	}
	else {
		res[0] = x0;
		if (y0 > y_r) res[1] = y_r + std::sqrt(r * r - (x0 - x_r) * (x0 - x_r));
		else res[1] = y_r - std::sqrt(r * r - (x0 - x_r) * (x0 - x_r));
	}
	return res;
}

double* Space::getCoord(int idx)
{
	double coord[2] = { *(*(bunch + idx)), *(*(bunch + idx) + 1) };
	return coord;
}

double* Space::getMoment(void)
{
	return multipole;
}

double* Space::getVoltage(void)
{
	return voltage;
}

void display(double*, std::pair<Space, double>);
void select(std::vector<std::pair<Space, double>>&);
void mutate(std::vector<std::pair<Space, double>>&);
double* evaluate(std::vector<std::pair<Space, double>>);
double calLoss(double *, double*);
void saveIndividual(std::vector<std::pair<Space, double>>, double, double);

void readMoment(std::string, double**);

int generation = 0;

int main(int argc, char* argv[]) {
	// Setting basic parameters
	size_t i, j, turn;

	Data readWave;
	std::string fname13, fname15;
	if (argv[1]) {
		fname13 = argv[1];
		if (argv[2]) { 
			fname15 = argv[2];
			double** mom = new double* [500];
			readMoment(fname13, mom);
			for (turn = std::atoi(argv[3]); turn < std::atoi(argv[4]); ++turn) {
				std::string dirName = "./result";
				while (!std::filesystem::create_directory(dirName)) {
					std::cout << "Directory already exist!" << std::endl;
					break;
				}
				std::ofstream data_csv = std::ofstream(dirName + "/data_result.csv");
				data_csv << "GEN,1-st,2-st,3-st,4-st,5-st,6-st,7-st,8-st,9-st,10-st,11-st,12-st,13-st,14-st,15-st,16-st,MeanLoss,BestLoss" << std::endl;

				// 0-step: Initialize Bunch
				//Space trueSpace(fname13);
				std::vector<std::pair<Space, double>> testSpace(POPULATION_SIZE);

				//trueSpace.setMoment(DIRECT_MODE);
				double* trueData = new double[16];
				for (i = 0; i < 16; ++i) *(trueData + i) = *(mom[turn] + i);
				for (i = 0; i < 16; ++i) {
					if (i == 0) data_csv << "TRUE," << *(trueData + i) << ",";
					else if (i < 15) data_csv << *(trueData + i) << ",";
					else if (i == 15) data_csv << *(trueData + i) << ",0,0" << std::endl;
				}
				double mx = trueData[1];
				double my = trueData[2];
				double sx = stdForTest(engine);
				double sy = std::sqrt(std::abs(sx * sx - (trueData[3] - (mx * mx - my * my))));

				// 1-step: Initialize Individual
				std::vector<std::thread> ths;
				for (i = 0; i < POPULATION_SIZE; ++i) {
					ths.emplace_back(
						[i, turn, &testSpace, trueData, mx, my, sx, sy]() {
							testSpace[i].first = Space(mx, sx, my, sy);
							testSpace[i].first.setMoment(DIRECT_MODE, mx, my);
							testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment());
						}
					);
				}
				for (std::thread& th : ths) th.join();
				ths.clear();
				std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
				for (i = 0; i < 16; ++i) {
					if (i == 0) data_csv << "INIT," << *(testSpace[0].first.getMoment() + i) << ",";
					else data_csv << *(testSpace[0].first.getMoment() + i) << ",";
				}
				data_csv << evaluate(testSpace)[0] << "," << evaluate(testSpace)[1] << std::endl;
				saveIndividual(testSpace, mx, my);
				select(testSpace);

				// 2-step: Main algorithm
				do {
					printf("\nGeneration %03u: Processing individual...\n", generation);
					for (i = 1; i < POPULATION_SIZE; ++i) {
						ths.emplace_back(
							[i, turn, mx, my, &testSpace, trueData]() {
								int count = 0;
								double prev_loss = testSpace[i].second;
								double now_loss;
								do {
									testSpace[i].first.transform(mode_selecting(engine), 1);
									testSpace[i].first.setMoment(DIRECT_MODE, mx, my);
									testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment());
									now_loss = testSpace[i].second;
									count++;
								} while (now_loss > prev_loss && count < 20);
							}
						);
					}
					for (std::thread& th : ths) th.join();
					ths.clear();
					select(testSpace);
					std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
					std::shuffle(testSpace.begin() + 1, testSpace.end(), engine);
					mutate(testSpace);
					std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
					for (i = 0; i < 16; ++i) {
						if (i == 0) data_csv << generation << "," << *(testSpace[0].first.getMoment() + i) << ",";
						else data_csv << *(testSpace[0].first.getMoment() + i) << ",";
					}
					data_csv << evaluate(testSpace)[0] << "," << testSpace[0].second << std::endl;
					display(trueData, testSpace[0]);
					generation++;
					std::shuffle(testSpace.begin() + 1, testSpace.end(), engine);
					saveIndividual(testSpace, mx, my);
				} while (generation < MAX_GENERATIONS);

				// 3-step: Display result
				display(trueData, testSpace[0]);
				std::string cmd = "py ./make_animation.py " + std::to_string(generation);
				std::system(cmd.c_str());
				std::cout << "Done !" << std::endl;
				generation = 0;
			}
		}
		else {
			std::string dirName = "./result";
			while (!std::filesystem::create_directory(dirName)) {
				std::cout << "Directory already exist!" << std::endl;
				break;
			}
			std::ofstream data_csv = std::ofstream(dirName + "/data_result.csv");
			data_csv << "GEN,1-st,2-st,3-st,4-st,5-st,6-st,7-st,8-st,9-st,10-st,11-st,12-st,13-st,14-st,15-st,16-st,MeanLoss,BestLoss" << std::endl;

			// 0-step: Initialize Bunch
			Space trueSpace(fname13);
			std::vector<std::pair<Space, double>> testSpace(POPULATION_SIZE);

			trueSpace.setMoment(DIRECT_MODE);
			double* trueData = trueSpace.getMoment();
			for (i = 0; i < 16; ++i) {
				if (i == 0) data_csv << "TRUE," << *(trueData + i) << ",";
				else if (i < 15) data_csv << *(trueData + i) << ",";
				else if (i == 15) data_csv << *(trueData + i) << ",0,0" << std::endl;
			}
			double mx = trueData[1];
			double my = trueData[2];
			double sx = stdForTest(engine);
			double sy = std::sqrt(std::abs(sx * sx - (trueData[3] - (mx * mx - my * my))));

			// 1-step: Initialize Individual
			std::vector<std::thread> ths;
			for (i = 0; i < POPULATION_SIZE; ++i) {
				ths.emplace_back(
					[i, turn, &testSpace, trueData, mx, my, sx, sy]() {
						testSpace[i].first = Space(mx, sx, my, sy);
						testSpace[i].first.setMoment(DIRECT_MODE, mx, my);
						testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment());
					}
				);
			}
			for (std::thread& th : ths) th.join();
			ths.clear();
			std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
			for (i = 0; i < 16; ++i) {
				if (i == 0) data_csv << "INIT," << *(testSpace[0].first.getMoment() + i) << ",";
				else data_csv << *(testSpace[0].first.getMoment() + i) << ",";
			}
			data_csv << evaluate(testSpace)[0] << "," << evaluate(testSpace)[1] << std::endl;
			saveIndividual(testSpace, mx, my);
			select(testSpace);

			// 2-step: Main algorithm
			do {
				printf("\nGeneration %03u: Processing individual...\n", generation);
				for (i = 1; i < POPULATION_SIZE; ++i) {
					ths.emplace_back(
						[i, turn, mx, my, &testSpace, trueData]() {
							int count = 0;
							double prev_loss = testSpace[i].second;
							double now_loss;
							do {
								testSpace[i].first.transform(mode_selecting(engine), 1);
								testSpace[i].first.setMoment(DIRECT_MODE, mx, my);
								testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment());
								now_loss = testSpace[i].second;
								count++;
							} while (now_loss > prev_loss && count < 20);
						}
					);
				}
				for (std::thread& th : ths) th.join();
				ths.clear();
				select(testSpace);
				std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
				std::shuffle(testSpace.begin() + 1, testSpace.end(), engine);
				mutate(testSpace);
				std::sort(testSpace.begin(), testSpace.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
				for (i = 0; i < 16; ++i) {
					if (i == 0) data_csv << generation << "," << *(testSpace[0].first.getMoment() + i) << ",";
					else data_csv << *(testSpace[0].first.getMoment() + i) << ",";
				}
				data_csv << evaluate(testSpace)[0] << "," << testSpace[0].second << std::endl;
				display(trueData, testSpace[0]);
				generation++;
				std::shuffle(testSpace.begin() + 1, testSpace.end(), engine);
				saveIndividual(testSpace, mx, my);
			} while (generation < MAX_GENERATIONS);

			// 3-step: Display result
			display(trueData, testSpace[0]);
			std::string cmd = "py ./make_animation.py " + std::to_string(generation);
			std::system(cmd.c_str());
			std::cout << "Done !" << std::endl;
		}
	}


	return 0;
}

/*
			for (i = 0; i < POPULATION_SIZE; ++i) {
				ths.emplace_back(
					[i, idx, d, turn, &testSpace, trueData]() {
						int count = 0;
						Space temp;
						for (size_t t = 0; t < BUNCH_SIZE; ++t) {
							*(*(temp.bunch + 1)) = *(*(testSpace[i].first.bunch + 1));
							*(*(temp.bunch + 1) + 1) = *(*(testSpace[i].first.bunch + 1) + 1);
						}
						double prev_loss = testSpace[i].second;
						double now_loss, temp_loss;
						do {
							int mode = mode_selecting(engine);
							testSpace[i].first.transform(mode, d);
							testSpace[i].first.setMoment(trueData[1], trueData[2]);
							testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment(), idx);
							now_loss = testSpace[i].second;
							count++;
						} while (now_loss > prev_loss && count < 20);
						if (now_loss > prev_loss) {
							for (size_t t = 0; t < BUNCH_SIZE; ++t) {
								*(*(testSpace[i].first.bunch + 1)) = *(*(temp.bunch + 1));
								*(*(testSpace[i].first.bunch + 1) + 1) = *(*(temp.bunch + 1) + 1);
							}
							testSpace[i].first.setMoment(trueData[1], trueData[2]);
							testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment(), idx);
						}
					}
				);
			}
* 
* 
			for (i = 0; i < POPULATION_SIZE; ++i) {
				ths.emplace_back(
					[i, idx, d, turn, &testSpace, trueData]() {
						int count = 0;
						Space temp;
						for (size_t t = 0; t < BUNCH_SIZE; ++t) {
							*(*(temp.bunch + 1)) = *(*(testSpace[i].first.bunch + 1));
							*(*(temp.bunch + 1) + 1) = *(*(testSpace[i].first.bunch + 1) + 1);
						}
						double prev_loss = testSpace[i].second;
						double now_loss, temp_loss;
						do {
							temp.transform(mode_selecting(engine), d);
							temp.setMoment(DIRECT_MODE);
							temp_loss = calLoss(trueData, temp.getMoment(), idx);
							now_loss = temp_loss;
							count++;
							if (now_loss > prev_loss) {
								*(*(temp.bunch + 1)) = *(*(testSpace[i].first.bunch + 1));
								*(*(temp.bunch + 1) + 1) = *(*(testSpace[i].first.bunch + 1) + 1);
							}
						} while (now_loss > prev_loss && count < 20);
						for (size_t t = 0; t < BUNCH_SIZE; ++t) {
							*(*(testSpace[i].first.bunch + 1)) = *(*(temp.bunch + 1));
							*(*(testSpace[i].first.bunch + 1) + 1) = *(*(temp.bunch + 1) + 1);
						}
						testSpace[i].first.setMoment(DIRECT_MODE);
						testSpace[i].second = calLoss(trueData, testSpace[i].first.getMoment(), idx);
					}
				);
			}
*/

double calLoss(double *true_b, double* test_b)
{
	double loss = 0;
	for (size_t j = 0; j < 16; ++j) {
		loss += std::pow(true_b[j] - test_b[j],2) / std::pow(8250, (int)((j + 1) / 2));
	} 
	return loss;
}

void saveIndividual(std::vector<std::pair<Space, double>> test_b, double mx, double my)
{
	size_t i, j;
	std::string chrName = "./result/individual_" + std::to_string(generation) + ".dat";
	std::ofstream fout;
	fout.open(chrName, std::ios::out | std::ios::binary | std::ios::trunc);
	for (j = 0; j < POPULATION_SIZE; ++j) {
		for (size_t i = 0; i < BUNCH_SIZE; ++i) {
			double x = test_b[j].first.bunch[i][0] + mx;
			fout.write((char*)&x, sizeof(double));
		}
		for (size_t i = 0; i < BUNCH_SIZE; ++i) {
			double y = test_b[j].first.bunch[i][1] + my;
			fout.write((char*)&y, sizeof(double));
		}
	}
	fout.close();
}

void readMoment(std::string fname, double** mom)
{
	std::ifstream is(fname, std::ifstream::binary);
	if (is) {
		is.seekg(0);
		for (size_t i = 0; i < 500; ++i) {
			mom[i] = new double[16];
			for (size_t j = 0; j < 16; ++j) {
				is.read((char*)&mom[i][j], sizeof(double));
			}
		}
	}
}

void mutate(std::vector<std::pair<Space, double>>& test_b) {
	std::vector<std::thread> ths;
	for (size_t i = 1; i < POPULATION_SIZE; ++i) {
		ths.emplace_back(
			[i, &test_b]() {
				printf("\rGeneration %03u: Mutating %u-st individual", generation, i);
				if (percentage(engine) < MUTATION_RATE) test_b[i].first.transform(mode_selecting(engine), 1);
			}
		);
	}
	for (std::thread& th : ths) th.join();
	ths.clear();
}

void select(std::vector<std::pair<Space, double>> &test_b) {
	size_t v_size = int(POPULATION_SIZE / TOURNAMENT_SIZE);
	size_t i;
	std::vector<std::thread> ths;
	std::shuffle(test_b.begin(), test_b.end(), engine);
	for (i = 0; i < v_size; ++i) {
		ths.emplace_back(
			[i, &test_b]() {
				printf("\rGeneration %03u: Selecting %u-st individual", generation, i);
				std::vector<std::pair<Space, double>> temp;
				for (size_t j = 0; j < TOURNAMENT_SIZE; ++j) temp.emplace_back(test_b[TOURNAMENT_SIZE * i + j]);
				std::sort(temp.begin(), temp.end(), [](std::pair<Space, double> const& iLhs, std::pair<Space, double> const& iRhs) {return iLhs.second < iRhs.second; });
				for (size_t j = 0; j < TOURNAMENT_SIZE; ++j) test_b[TOURNAMENT_SIZE * i + j] = temp[0];
				temp.clear();
			}
		);
	}
	for (std::thread& th : ths) th.join();
	ths.clear();
	std::cout << std::endl;
}

double* evaluate(std::vector<std::pair<Space, double>> test_b)
{
	double res[2];
	double mean = 0;
	for (int i = 0; i < POPULATION_SIZE; ++i) mean += test_b[i].second/POPULATION_SIZE;
	res[0] = mean; res[1] = test_b[0].second;
	return res;
}

void display(double* true_b, std::pair<Space, double> test_b) {
	printf("\n** Best Loss: %1.10e\n", test_b.second);
	std::cout << "00  True data\tBest data" << std::endl;
	for (size_t i = 0; i < 16; ++i) {
		printf("%02u %1.3e\t%1.3e\t%3.3f%%\n", i, true_b[i], test_b.first.getMoment()[i], 100 * std::abs((true_b[i] - test_b.first.getMoment()[i]) / true_b[i]));
	}
	return;
}