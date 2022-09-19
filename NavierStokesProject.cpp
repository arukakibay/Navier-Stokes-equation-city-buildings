#include "pch.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;

int main() {
	int const xSize = 101; //x
	int const ySize = 51; //y
	int fullIterations = 0;
	int pressureIterations = 0;
	int endTime, startTime;
	double U[xSize][ySize],
		U_star[xSize][ySize],
		oldU[xSize][ySize],
		V[xSize][ySize],
		V_star[xSize][ySize],
		oldV[xSize][ySize],
		P[xSize][ySize],
		oldP[xSize][ySize];

	double dt, dx, dy, Re, errorP, errorSpeed, differenceU, differenceV, differenceP, rho;

	startTime = clock();
	Re = 10.0;
	dx = 1.0 * pow(xSize - 1, -1);
	dy = 1.0 * pow(xSize - 1, -1);
	dt = (Re*dx*dy)*0.25;
	errorP = pow(10, -5);
	errorSpeed = pow(10, -12);
	rho = 1.2;
	ofstream fout("Navie4.dat", ios::out);
	fout << "variables=\"X\" ,\"Y\" ,\"U\" ,\"V\" ,\"P\"" << endl;

	//fill by zeros
	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {
			U[i][j] = 0.0;
			V[i][j] = 0.0;
			oldU[i][j] = 0.0;
			oldV[i][j] = 0.0;
			U_star[i][j] = 0.0;
			V_star[i][j] = 0.0;
			oldP[i][j] = 0.0;
			P[i][j] = 0.0;
		}
	}

	for (int j = ySize / 2; j < ySize; j++) {
		oldU[0][j] = 0.05; //inlet
		oldU[xSize - 1][j] = oldU[xSize - 2][j]; //outlet
		oldV[xSize - 1][j] = oldV[xSize - 2][j];
	}

	//house 1 (without roof)
	for (int i = 22; i <= 42; i++) {
		for (int j = 0; j <= 25; j++) {
			oldU[i][j] = 0.0;
			oldV[i][j] = 0.0;
			oldP[i][j] = 0.0;
			P[i][j] = 0.0;
		}
	}
	// roof of house 1
	for (int i = 22; i <= 42; i++) {
		for (int j = 25; j <= 35; j++) {
			if (j - i <= 3 && i + j <= 67) {
				oldU[i][j] = 0.0;
				oldV[i][j] = 0.0;
				oldP[i][j] = 0.0;
				P[i][j] = 0.0;
			}
		}
	}

	// house 2
	for (int i = 70; i <= 85; i++) {
		for (int j = 0; j <= 20; j++) {
			oldU[i][j] = 0.0;
			oldV[i][j] = 0.0;
			oldP[i][j] = 0.0;
			P[i][j] = 0.0;
		}
	}

	do {
		for (int i = 0; i <= xSize - 1; i++) {
			// boundary conditions of Dirihlet, on walls
			for (int j = 0; j <= ySize - 1; j++) {
				U_star[xSize - 1][j] = 0.0;
				U_star[i][0] = 0.0;
				U_star[0][j] = 0.0;
				U_star[i][ySize - 1] = 0.0;

				V_star[i][ySize - 1] = 0.0;
				V_star[xSize - 1][j] = 0.0;
				V_star[i][0] = 0.0;
				V_star[0][j] = 0.0;

			}
		}
		//house 1 without roof
		for (int i = 22; i <= 42; i++) {
			for (int j = 0; j <= 25; j++) {
				U_star[i][j] = 0.0;
				V_star[i][j] = 0.0;
			}
		}
		//roof of house 1
		for (int i = 22; i <= 42; i++) {
			for (int j = 25; j <= 35; j++) {
				if (j - i <= 3 && i + j <= 67) {
					U_star[i][j] = 0.0;
					V_star[i][j] = 0.0;
				}
			}
		}

		//house 2
		for (int i = 70; i <= 85; i++) {
			for (int j = 0; j <= 20; j++) {
				U_star[i][j] = 0.0;
				V_star[i][j] = 0.0;
			}
		}

		for (int j = ySize / 2; j < ySize; j++) {
			U_star[0][j] = 0.05; //inlet
		}
		for (int j = ySize / 2; j < ySize; j++) {
			U_star[xSize - 1][j] = U_star[xSize - 2][j]; //outlet
			V_star[xSize - 1][j] = V_star[xSize - 2][j];
		}
		//Method of splitting by physical parameters 
		for (int i = 1; i < xSize - 1; i++)
		{
			for (int j = 1; j < ySize - 1; j++) {
				if (i >= 22 && i <= 42 && j >= 0 && j <= 25) {
					U_star[i][j] = 0.0;
					V_star[i][j] = 0.0;
				}
				else if (i >= 22 && i <= 42 && j >= 25 && j <= 35 && j - i <= 3 && i + j <= 67) {
					U_star[i][j] = 0.0;
					V_star[i][j] = 0.0;
				}
				else if (i >= 70 && i <= 85 && j >= 0 && j <= 20) {
					U_star[i][j] = 0.0;
					V_star[i][j] = 0.0;
				}
				else {
					U_star[i][j] = oldU[i][j] +
						dt * (
						(-oldU[i][j] *
							(oldU[i + 1][j] - oldU[i - 1][j]) / (2 * dx) -
							oldV[i][j] *
							(oldU[i][j + 1] - oldU[i][j - 1]) / (2 * dy))
							+ (
							(oldU[i + 1][j] - 2.0 * oldU[i][j] + oldU[i - 1][j]) / (dx * dx)
								+ (oldU[i][j + 1] - 2.0 * oldU[i][j] + oldU[i][j - 1]) / (dy* dy)
								) / Re
							);

					V_star[i][j] = oldV[i][j] + dt *
						(
						(-oldU[i][j] * (oldV[i + 1][j] - oldV[i - 1][j]) / (2 * dx)
							- oldV[i][j] * (oldV[i][j + 1] - oldV[i][j - 1]) / (2 * dy))
							+ (
							(oldV[i + 1][j] - 2.0 * oldV[i][j] + oldV[i - 1][j]) / (dx* dx)
								+ (oldV[i][j + 1] - 2.0 * oldV[i][j] + oldV[i][j - 1]) / (dy* dy)
								) / Re);
				}
			}
		}


		pressureIterations = 0;
		double prev_differenceP = 0;
		do {

			for (int j = (ySize - 1) / 2; j < ySize - 1; j++) {
				P[xSize - 1][j] = 0; //Outlet
			}

			//The left and right walls of 1st house
			for (int j = 0; j <= 25; j++) {
				P[22][j] = P[21][j];
				P[42][j] = P[43][j];
			}

			//The left and right walls of 2nd house
			for (int j = 0; j <= 20; j++) {
				P[70][j] = P[69][j];
				P[85][j] = P[86][j];
			}

			//Roof of 2nd house
			for (int i = 70; i <= 85; i++) {
				P[i][20] = P[i][21];
			}
			//Roof of 1st house
			for (int i = 22; i <= 42; i++) {
				for (int j = 25; j <= 35; j++) {
					if (j - i == 3) {
						P[i][j] = P[i - 1][j + 1];
					}
					else if (i + j == 67) {
						P[i][j] = P[i + 1][j + 1];
					}
				}
			}
			//Neumann conditions
			for (int i = 0; i <= xSize - 1; i++) {
				for (int j = 0; j <= (ySize - 1) / 2; j++) {
					P[0][j] = P[1][j];
					P[i][ySize - 1] = P[i][ySize - 2];
					P[i][0] = P[i][1];
					P[xSize - 1][j] = P[xSize - 2][j];
				}
			}


			//Poisson equation
			for (int i = 1; i < xSize - 1; i++) {
				for (int j = 1; j < ySize - 1; j++) {
					if (i >= 22 && i <= 42 && j >= 0 && j <= 25 && j - i <= 3 && i + j <= 67) {
						P[i][j] = 0.0;
					}
					else if (i >= 22 && i <= 42 && j >= 25 && j <= 35 && j - i <= 3 && i + j <= 67) {
						P[i][j] = 0.0;
					}
					else if (i >= 70 && i <= 85 && j >= 0 && j <= 20) {
						P[i][j] = 0.0;
					}
					else {
						P[i][j] = 0.25 * (oldP[i + 1][j] + oldP[i - 1][j] + oldP[i][j + 1] + oldP[i][j - 1]
							- (rho*(U_star[i + 1][j] - U_star[i - 1][j]) / (2 * dx*dt)
								+ rho * (V_star[i][j + 1] - V_star[i][j - 1]) / (2 * dy*dt))*dy*dx);
					}
				}
			}

			//left and right walls of 1st house 
			for (int j = 0; j <= 25; j++) {
				P[22][j] = P[21][j];
				P[42][j] = P[43][j];
			}

			//left and right walls of 2nd house 
			for (int j = 0; j <= 20; j++) {
				P[70][j] = P[69][j];
				P[85][j] = P[86][j];
			}

			//roof of 2nd house
			for (int i = 70; i <= 85; i++) {
				P[i][20] = P[i][21];
			}

			//roof of 1st house
			for (int i = 22; i <= 42; i++) {
				for (int j = 25; j <= 35; j++) {
					if (j - i == 3) {
						P[i][j] = P[i - 1][j + 1];
					}
					else if (i + j == 67) {
						P[i][j] = P[i + 1][j + 1];
					}
				}
			}

			differenceP = 0.0;
			for (int i = 0; i < xSize; i++) {
				for (int j = 0; j < ySize; j++) {
					if (differenceP < fabs(P[i][j] - oldP[i][j])) {
						differenceP = fabs(P[i][j] - oldP[i][j]);
					}
				}
			}
			for (int i = 0; i < xSize; i++) {
				for (int j = 0; j < ySize; j++) {
					oldP[i][j] = P[i][j];
				}
			}
			pressureIterations++;
		} while (differenceP > errorP);

		//2nd part 
		for (int i = 1; i < xSize - 1; i++)
		{
			for (int j = 1; j < ySize - 1; j++)
			{
				if (i >= 22 && i <= 42 && j >= 0 && j <= 25) {
					U[i][j] = 0.0;
					V[i][j] = 0.0;
				}
				else if (i >= 22 && i <= 42 && j >= 25 && j <= 35 && j - i <= 3 && i + j <= 67) {
					U[i][j] = 0.0;
					V[i][j] = 0.0;
				}
				else if (i >= 70 && i <= 85 && j >= 0 && j <= 20) {
					U[i][j] = 0.0;
					V[i][j] = 0.0;
				}
				else {
					U[i][j] = U_star[i][j] - dt * (oldP[i + 1][j] - oldP[i - 1][j]) / (2 * rho*dx);
					V[i][j] = V_star[i][j] - dt * (oldP[i][j + 1] - oldP[i][j - 1]) / (2 * rho*dy);
				}
			}
		}
		//Dirichlet
		for (int i = 0; i < xSize; i++) {

			for (int j = 0; j < ySize; j++) {
				U[i][0] = 0;
				V[i][0] = 0;

				U[i][ySize - 1] = 0;
				V[i][ySize - 1] = 0;

				U[0][j] = 0;
				V[0][j] = 0;

				U[xSize - 1][j] = 0;
				V[xSize - 1][j] = 0;
			}
		}
		//House 1(without roof)
		for (int i = 22; i <= 42; i++) {
			for (int j = 0; j <= 25; j++) {
				U[i][j] = 0.0;
				V[i][j] = 0.0;
			}
		}
		//roof of 1st house
		for (int i = 22; i <= 42; i++) {
			for (int j = 25; j <= 35; j++) {
				if (j - i <= 3 && i + j <= 67) {
					U[i][j] = 0.0;
					V[i][j] = 0.0;
				}
			}
		}

		//house 2
		for (int i = 70; i <= 85; i++) {
			for (int j = 0; j <= 20; j++) {
				U[i][j] = 0.0;
				V[i][j] = 0.0;
			}
		}

		for (int j = ySize / 2; j < ySize; j++) {
			U[0][j] = 0.05; //inlet
			U[xSize - 1][j] = oldU[xSize - 2][j]; //outlet
			V[xSize - 1][j] = V[xSize - 2][j];
		}

		differenceU = 0.0;
		differenceV = 0.0;
		for (int i = 0; i < xSize; i++) {
			for (int j = 0; j < ySize; j++) {
				if (differenceU < fabs(U[i][j] - oldU[i][j]))
				{
					differenceU = fabs(U[i][j] - oldU[i][j]);
				}
				if (differenceV < fabs(V[i][j] - oldV[i][j]))
				{
					differenceV = fabs(V[i][j] - oldV[i][j]);
				}
			}
		}
		for (int i = 0; i < xSize; i++) {
			for (int j = 0; j < ySize; j++) {
				oldU[i][j] = U[i][j];
				oldV[i][j] = V[i][j];
			}
		}
		if (fullIterations % 10 == 0) {
			fout << "zone T = \"" << fullIterations % 10 << "\", i=" << xSize << ", j=" << ySize << ", F=point" << endl;
			for (int j = 0; j < ySize; j++) {
				for (int i = 0; i < xSize; i++) {
					fout << i * dx << "\t" << j * dy << "\t" << U[i][j] << "\t" << V[i][j] << "\t" << P[i][j] << "\t" << endl;
				}
			}
		}
		fullIterations++;
		cout << fullIterations << endl;
	} while (differenceU > errorSpeed || differenceV > errorSpeed);

	fout << "zone T = \"" << fullIterations % 10 + 1 << "\", i=" << xSize << ", j=" << ySize << ", F=point" << endl;
	for (int j = 0; j < ySize; j++) {
		for (int i = 0; i < xSize; i++) {
			fout << i * dx << "\t" << j * dy << "\t" << U[i][j] << "\t" << V[i][j] << "\t" << P[i][j] << "\t" << endl;
		}
	}

	fout.close();
	endTime = clock();
	cout << "number of interations= " << fullIterations << endl;


	system("pause"); return 0; }
