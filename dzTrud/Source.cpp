#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#define DELTA 0.01
#define TIME_START 0.
#define TIME_END 10
#define EPS 1e-6

#define N 13

class Matrix {
private:
	double** matrix;

	double countExpDif(double phiA, double phiK) {
		return 1e-12 / 0.026 * exp((phiA - phiK) / 0.026);
	}
public:
	Matrix() {
		matrix = new double*[N];
		for (int i = 0; i < N; i++) {
			matrix[i] = new double[N];
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				matrix[i][j] = 0.;
			}
		}
	}

	void setMatrix(double* vectorY) {
		double phi3 = vectorY[10];
		double phi4 = vectorY[11];

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				matrix[i][j] = 0.;
			}
		}

		// Инициализация производных переменных состояния
		for (int i = 0; i < 4; i++) {
			matrix[i][i] = 1;
			matrix[i][i + 4] = -1. / DELTA;
		}

		// Инициализация диагональных значений напряжений ёмкости
		for (int i = 4; i < 7; i++) {
			matrix[i][i] = 1.;
		}

		// Инициализация недиагональных значений напряжений ёмкости
		matrix[4][9] = -1.;
		matrix[5][11] = -1.;
		matrix[6][10] = -1.;
		matrix[6][11] = 1.;

		// Инициализация недиагональных значений токов индуктивности
		matrix[7][3] = 1e-3;
		matrix[7][8] = -1.;

		// Инициализация токов в узлах
		// phi1
		matrix[8][7] = -1.;
		matrix[8][12] = -1.;

		// phi2
		matrix[9][0] = 1e-3;
		matrix[9][9] = 1. / 20.;
		matrix[9][10] = -1. / 20.;
		matrix[9][12] = 1.;

		// phi3
		matrix[10][2] = 2e-12;
		matrix[10][9] = -1. / 20.;
		matrix[10][10] = 1. / 20. + countExpDif(phi3, phi4) + 1. / 1e6;
		matrix[10][11] = -countExpDif(phi3, phi4) - 1. / 1e6;

		// phi4
		matrix[11][1] = 1e-3;
		matrix[11][2] = -2e-12;
		matrix[11][10] = -countExpDif(phi3, phi4) - 1. / 1e6;
		matrix[11][11] = countExpDif(phi3, phi4) + 1. / 1e6;

		// Инициализация тока на источнике типа E
		matrix[12][8] = -1.;
		matrix[12][9] = 1.;
	}

	double** getMatrix() {
		return matrix;
	}
};

class VectorB {
private:
	double* vector;

	double countExp(double phiA, double phiK) {
		return 1e-12 * (exp((phiA - phiK) / 0.026) - 1.);
	}
public:
	VectorB() {
		vector = new double[N];
		
		for (int i = 0; i < N; i++) {
			vector[i] = 0.;
		}
	}

	void setVectorB(double* vectorY, double* vectorPrev, double t) {
		// Производные переменных состояния
		for (int i = 0; i < 4; i++) {
			vector[i] = vectorY[i] - (vectorY[i + 4] - vectorPrev[i]) / DELTA;
		}
		
		// Переменные состояния
		vector[4] = vectorY[4] - vectorY[9];
		vector[5] = vectorY[5] - vectorY[11];
		vector[6] = vectorY[6] - (vectorY[10] - vectorY[11]);
		vector[7] = 1e-3 * vectorY[3] - vectorY[8];

		// Узлы
		vector[8] = -vectorY[7] - vectorY[12];
		vector[9] = vectorY[12] + 1e-3 * vectorY[0] + (vectorY[9] - vectorY[10]) / 20.;
		vector[10] = -(vectorY[9] - vectorY[10]) / 20. + countExp(vectorY[10], vectorY[11]) + 2e-12 * vectorY[2] + 
			(vectorY[10] - vectorY[11]) / 1e6;
		vector[11] = -countExp(vectorY[10], vectorY[11]) - 2e-12 * vectorY[2] - (vectorY[10] - vectorY[11]) / 1e6 +
			1e-3 * vectorY[1];

		// Источник типа E
		vector[12] = vectorY[9] - vectorY[8] - sin(2 * M_PI * t);
	}

	double* getVectorB() {
		return vector;
	}

	void invertVectorB() {
		for (int i = 0; i < N; i++) {
			vector[i] = -vector[i];
		}
	}
};

class VectorY {
private:
	double* vector;
public:
	VectorY() {
		vector = new double[N];
		setZero();
	}

	void setZero() {
		for (int i = 0; i < N; i++) {
			vector[i] = 0.;
		}
	}

	void setVectorY(double* vectorY) {
		for (int i = 0; i < N; i++) {
			vector[i] = vectorY[i];
		}
	}

	void changeVectorY(double* vectorY) {
		for (int i = 0; i < N; i++) {
			vector[i] += vectorY[i];
		}
	}

	double* getVectorY() {
		return vector;
	}

	bool compareToZero() {
		for (int i = 0; i < N; i++) {
			if (fabs(vector[i]) > EPS) {
				return false;
			}
		}
		return true;
	}
};

class VectorPrev {
private:
	double* vector;
public:
	VectorPrev() {
		vector = new double[4];

		for (int i = 0; i < 4; i++) {
			vector[i] = 0.;
		}
	}

	void setVectorPrev(double Uc1, double Uc2, double Ucb, double Il) {
		vector[0] = Uc1;
		vector[1] = Uc2;
		vector[2] = Ucb;
		vector[3] = Il;
	}

	double* getVectorPrev() {
		return vector;
	}
};

double* gauss(double **matrA, double *vecB) {
	int i, j, row;
	double k;
	double* u = new double[N];

	for (row = 1; row < N; row++)
	{
		for (i = row; i < N; i++) {
			k = matrA[i][row - 1] / matrA[row - 1][row - 1];
			for (j = 0; j < N; j++) {
				matrA[i][j] -= k * matrA[row - 1][j];
			}
			vecB[i] -= k * vecB[row - 1];
		}
	}

	for (row = 0; row < N; row++) {
		k = matrA[row][row];
		for (j = row; j < N; j++) {
			matrA[row][j] /= k;
		}
		vecB[row] /= k;
	}

	for (row = N - 1; row > -1; row--)
	{
		for (i = N - 1; i > row; i--) {
			vecB[row] -= matrA[row][i] * u[i];
		}
		u[row] = vecB[row];
	}
	return u;
}

int main() {
	Matrix matrixA = Matrix();
	VectorB vectorB = VectorB();
	VectorY vectorY = VectorY();
	VectorY deltaVectorY = VectorY();
	VectorPrev prev = VectorPrev();

	// Инициализация НУ
	double t = TIME_START;
	prev.setVectorPrev(0, 0, 0, 0);
	ofstream fout("output.csv");
	fout << "t; U;" << endl;

	// Основной цикл
	while (t < TIME_END) {
		vectorY.setZero();
		do {
			deltaVectorY.setZero();
			vectorB.setVectorB(vectorY.getVectorY(), prev.getVectorPrev(), t);
			vectorB.invertVectorB();
			matrixA.setMatrix(vectorY.getVectorY());
			deltaVectorY.setVectorY(gauss(matrixA.getMatrix(), vectorB.getVectorB()));
			vectorY.changeVectorY(deltaVectorY.getVectorY());
		} while (!deltaVectorY.compareToZero());
		double* y = vectorY.getVectorY();
		prev.setVectorPrev(y[4], y[5], y[6], y[7]);
		fout << t << "; " << y[11] << ";" << endl;
		t += DELTA;
	}

	fout.close();
}