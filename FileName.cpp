#include <iostream>
#include <vector>
//#include <vector>
#include <cmath>
//#include "E_f.h"
#include <fstream>
#include <random>


// #include "Cmake_dip.h" 

// „исло частиц
// int N = 10;

double L = 10.0;
double dx = 0.1; //шаг по сетке
int NX = L / dx; // кол-во €чеек

//Ёлектрическое поле
double E0 = pow(10, -12);

// „исло частиц
int N = 1000;

//ќсновные константы
double m_ion = 1.673 * pow(10, -27);
double m_el = 9.31 * pow(10, -31);
double q = 1.6 * pow(10, -19);


// шаг по времени и врем€
double dt = 0.1;
double T = 10.0;



/*
double E(double x) {
	return (E0 * (x - L / 2));
}
*/

double E(double E1, double E2, double x) {

	//int n1 = x1 / dx;
	//double x_loc = L - n1 * dx;

	//int n2 = n1 + 1;
	//double x2_loc = L - n2 * dx;
	return  (E1 * x + E2 * (dx - x)) / dx;
}

void Move(std::vector<double> V_ions, std::vector<double> V_el, std::vector<double> X_ions, std::vector<double> X_el) {

	std::vector <double> E_vec;
	int num_ceil = L / dx + 1;
	for (int i = 0; i - 1 < num_ceil; i++) {
		E_vec.push_back(E0 * (i * dx - L / 2));
	}
	std::ofstream out_el;          // поток дл€ записи
	out_el.open("data_N_el.txt");

	std::ofstream out_ions;          // поток дл€ записи
	out_ions.open("data_N_ions.txt");

	for (double tim = 0.0; tim < T; tim += dt) {
		
		std::vector<int> n_el(NX, 0);
		std::vector<int> n_ions(NX, 0);

	

		for (std::size_t i = 0; i < X_ions.size(); i++) {

			
			// ѕериодические граничные услови€ 
			/*
			
			if (X_ions[i] >= L) { X_ions[i] -= L; };
			if (X_ions[i] < 0) { X_ions[i] += L; };

			if (X_el[i] >= L) { X_el[i] -= L; };
			if (X_el[i] < 0) { X_el[i] += L; };
			*/

			//∆естка€ стенка
			if(X_ions[i] >= L){
				X_ions[i] = L - 2 * (X_ions[i] - L);
				V_ions[i] *= (-1);
			}
			if(X_ions[i] < 0) {
				X_ions[i] *= (-1);
				V_ions[i] *= (-1);
			}

			if (X_el[i] >= L) {
				X_el[i] = L - 2 * (X_el[i] - L);
				V_el[i] *= (-1);
			}
			if (X_el[i] < 0) {
				X_el[i] *= (-1);
				V_el[i] *= (-1);
			}







			// out_el << X_el[i] << " ";
			// out_ions << X_ions[i] << " ";
			int ceil_ion = X_ions[i] / dx;
			int ceil_el = X_el[i] / dx;

			double x_ion_loc = fmod(X_ions[i], dx * 1.0);
			double x_el_loc = fmod(X_el[i], dx * 1.0);


			V_ions[i] += E(E_vec[ceil_ion], E_vec[ceil_ion + 1], x_ion_loc) * q / m_ion * dt;
			V_el[i] -= E(E_vec[ceil_el], E_vec[ceil_el + 1], x_el_loc)* q / m_el * dt;
			X_ions[i] += V_ions[i] * dt;
			X_el[i] += V_el[i] * dt;


			n_el[ceil_el] += 1;
			n_ions[ceil_ion] += 1;
			
		}
		for (std::size_t p = 0; p < n_el.size(); p++) {
			out_el << n_el[p] << " ";
			out_ions << n_ions[p] << " ";
		}

		out_el << std::endl;
		out_ions << std::endl;
	}


	out_el.close();
	out_ions.close();

};


int main() {

	std::vector<double> X_ions;
	std::vector<double> X_el;
	std::vector<double> V_ions(N, 0.0);
	std::vector<double> V_el(N, 0.0);

	double lower_bound = 0.0;
	double upper_bound = L;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::default_random_engine re;


	for (int i = 0; i < N; i++) {

		double a_random_double = unif(re);
		// int L_int = L;
		X_ions.push_back(a_random_double);
	}

	for (int i = 0; i < N; i++) {

		double a_random_double = unif(re);
		// int L_int = L;
		X_el.push_back(a_random_double);
	}

	Move(V_ions, V_el, X_ions, X_el);


}