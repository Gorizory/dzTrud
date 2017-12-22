#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#define DELTA 1e-6
#define TIME_START 0.
#define TIME_END 5e-6
#define EPS 1e-6
#define EPSILON2
#define EPSILON1

#define N 13

#define L 1e-3
#define C1 1e-6
#define C2 1e-4
#define Rb 20
#define Ru 1e6
#define Cb 2e-12
#define It 1e-12
#define M_Phit 0.026



void printVector(double *A, int size)
{
	for (int j = 0; j < size; j++)
	{
		printf("% .3e\n", A[j]);
	}
	puts("");
}

void printMatrix(double **A, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			printf("% .3e\t", A[i][j]);
		}
		putchar('\n');
	}
	puts("");
}








class Matrix {
private:
	double** matrix;

	double countExpDif(double phiA, double phiK) {
		return It / M_Phit * exp((phiA - phiK) / M_Phit);
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

		// U', I'
		for (int i = 0; i < 4; i++) {
			matrix[i][i] = 1.;
			matrix[i][i + 4] = -1. / DELTA;
		}

		// U diagonal
		for (int i = 4; i < 7; i++) {
			matrix[i][i] = 1.;
		}

		// U non-diag
		matrix[4][9] = -1.;
		matrix[5][11] = -1.;
		matrix[6][10] = -1.;
		matrix[6][11] = 1.;

		// I non-diag
		matrix[7][3] = L;
		matrix[7][8] = -1.;

		// Nodes
		// phi1
		matrix[8][7] = -1.;
		matrix[8][12] = -1.;

		// phi2
		matrix[9][0] = C1;
		matrix[9][9] = 1. / Rb;
		matrix[9][10] = -1. / Rb;
		matrix[9][12] = 1.;

		// phi3
		matrix[10][2] = Cb;
		matrix[10][9] = -1. / Rb;
		matrix[10][10] = 1. / Rb + countExpDif(phi3, phi4) + 1. / Ru;
		matrix[10][11] = -countExpDif(phi3, phi4) - 1. / Ru;

		// phi4
		matrix[11][1] = C2;
		matrix[11][2] = -Cb;
		matrix[11][10] = -countExpDif(phi3, phi4) - 1. / Ru;
		matrix[11][11] = countExpDif(phi3, phi4) + 1. / Ru;

		// E
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
		return It * (exp((phiA - phiK) / M_Phit) - 1.);
	}
public:
	VectorB() {
		vector = new double[N];

		for (int i = 0; i < N; i++) {
			vector[i] = 0.;
		}
	}

	void setVectorB(double* vectorY, double* vectorPrev, double t) {
		// U', I'
		for (int i = 0; i < 4; i++) {
			vector[i] = vectorY[i] - (vectorY[i + 4] - vectorPrev[i]) / DELTA;
		}

		// U, I
		vector[4] = vectorY[4] - vectorY[9];
		vector[5] = vectorY[5] - vectorY[11];
		vector[6] = vectorY[6] - (vectorY[10] - vectorY[11]);
		vector[7] = L * vectorY[3] - vectorY[8];

		// Phi
		vector[8] = -vectorY[7] - vectorY[12];
		vector[9] = vectorY[12] + C1 * vectorY[0] + (vectorY[9] - vectorY[10]) / Rb;
		vector[10] = -(vectorY[9] - vectorY[10]) / Rb + countExp(vectorY[10], vectorY[11]) + Cb * vectorY[2] +
			(vectorY[10] - vectorY[11]) / Ru;
		vector[11] = -countExp(vectorY[10], vectorY[11]) - Cb * vectorY[2] - (vectorY[10] - vectorY[11]) / Ru +
			C2 * vectorY[1];

		// E
		vector[12] = vectorY[9] - vectorY[8] - 10. * sin(2 * M_PI / 1e4 * t);
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

class Gauss {
private:
	double **a, *y;
	void copy(double **A, double *Y) {
		for (int i = 0; i < N; i++) {
			y[i] = Y[i];
			for (int j = 0; j < N; j++)
				a[i][j] = A[i][j];
		}
	}
public:
	Gauss() {
		a = new double*[N];
		y = new double[N];
		for (int i = 0; i < N; i++)
			a[i] = new double[N];
	}

	double *solve(double **A, double *Y) {
		copy(A, Y);
		double *x, max;
		int k, index;
		x = new double[N];
		k = 0;
		while (k < N) {
			max = abs(a[k][k]);
			index = k;
			for (int i = k + 1; i < N; i++) {
				if (abs(a[i][k]) > max) {
					max = abs(a[i][k]);
					index = i;
				}
			}
			for (int j = 0; j < N; j++) {
				double temp = a[k][j];
				a[k][j] = a[index][j];
				a[index][j] = temp;
			}
			double temp = y[k];
			y[k] = y[index];
			y[index] = temp;
			for (int i = k; i < N; i++) {
				double temp = a[i][k];
				if (abs(temp) < EPS) continue;
				for (int j = 0; j < N; j++)
					a[i][j] = a[i][j] / temp;
				y[i] = y[i] / temp;
				if (i == k)  continue;
				for (int j = 0; j < N; j++)
					a[i][j] = a[i][j] - a[k][j];
				y[i] = y[i] - y[k];
			}
			k++;
		}
		for (k = N - 1; k >= 0; k--) {
			x[k] = y[k];
			for (int i = 0; i < k; i++)
				y[i] = y[i] - a[i][k] * x[k];
		}
		return x;
	}
};

int main() {
	Matrix matrixA;
	VectorB vectorB;
	VectorY vectorY;
	VectorY deltaVectorY;
	VectorPrev prev;
	Gauss gauss;

	double t = TIME_START;
	prev.setVectorPrev(0, 0, 0, 0);
	ofstream fout("output.csv");
	fout << "t; U;" << endl;

	while (t < TIME_END) {
		do {
			deltaVectorY.setZero();
			vectorB.setVectorB(vectorY.getVectorY(), prev.getVectorPrev(), t);
			vectorB.invertVectorB();
			matrixA.setMatrix(vectorY.getVectorY());
			deltaVectorY.setVectorY(gauss.solve(matrixA.getMatrix(), vectorB.getVectorB()));
			printMatrix(matrixA.getMatrix(), N);
			vectorY.changeVectorY(deltaVectorY.getVectorY());
		} while (!deltaVectorY.compareToZero());
		double* y = vectorY.getVectorY();
		prev.setVectorPrev(y[4], y[5], y[6], y[7]);
		fout << t << "; " << y[11] << ";" << endl;
		t += DELTA;
	}

	fout.close();
	getchar();
}