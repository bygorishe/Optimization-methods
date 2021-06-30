#include "func.h"

random_device rd;   //Вихрь Мерсенна
mt19937 random(rd());

int N = 204330250, M = 1000, iter;

void newRandomPoint(vector<double> &x) {
	uniform_real_distribution<double> range(-10, 10);
	x[0] = range(random);
	x[1] = range(random);
}

void simpleRandom(vector<double> &x0) {
	num_f = 0;
	newRandomPoint(x0);
	vector<double> x1(2);
	double f0 = Q(x0), f1;
	for (iter = 0; iter < N; iter++) {
		newRandomPoint(x1);
		f1 = Q(x1);
		num_f++;
		if (f1 < f0) {
			f0 = f1;
			x0 = x1;
		}
	}
}

void globalSearch1(vector<double>& x0) {
	num_f = 0;
	newRandomPoint(x0);
	vector<double> x1(2);
	double f0, f1;
	Rosenbrok(x0);
	f0 = Q(x0);
	num_f++;
	for (iter = 0; iter < M; iter++) {
		newRandomPoint(x1);
		Rosenbrok(x1);
		f1 = Q(x1);
		num_f++;
		if (f1 < f0) {
			f0 = f1;
			x0 = x1;
		};
	}
}

void globalSearch2_v1(vector<double>& x0) {
	num_f = 0;
	newRandomPoint(x0);
	vector<double> x1(2);
	double f0, f1;
	Rosenbrok(x0);
	f0 = Q(x0);
	num_f++;
	for (iter = 0; iter < M; iter++) {
		newRandomPoint(x1);
		f1 = Q(x1);
		num_f++;
		if (f1 < f0) {
			x0 = x1;
			Rosenbrok(x0);
			f0 = Q(x0);
			num_f++;
		}
	}
}

void globalSearch2_v2(vector<double>& x0) {
	num_f = 0;
	iter = 0;
	newRandomPoint(x0);
	vector<double> x1(2);
	double f0, f1;
	Rosenbrok(x0);
	f0 = Q(x0);
	num_f++;
	while (iter < M) {
		for (iter = 0; iter < M; iter++) {
			newRandomPoint(x1);
			f1 = Q(x1);
			num_f++;
			if (f1 < f0) {
				x0 = x1;
				break;
			}
		}
		Rosenbrok(x0);
		f0 = Q(x0);
		num_f++;
	}
}

void globalSearch3(vector<double>& x0) {
	num_f = 0;
	newRandomPoint(x0);
	vector<double> x1(2), s(2), s_temp(2);
	double f0, f1;
	Rosenbrok(x0);
	f0 = Q(x0);
	num_f++;
	x1 = x0;
	for (iter = 0; iter < M; iter++) {
		uniform_real_distribution<double> range(-1, 1);
		s_temp[0] = range(random);
		s_temp[1] = range(random);
		s = { 0,0 };
		int itr = 0;
		do{
			s += s_temp;
			x1 += s;
			f1 = Q(x1);
			num_f++;
			itr++;
		} while (f0 < f1 && itr < 1000);

		if (f1 < f0) {
			Rosenbrok(x1);
			f1 = Q(x1);
			num_f++;
			if (f1 < f0)
				x0 = x1;
		}
		else
			x1 = x0;
	}
}

void globalSearch(vector<double>& x0) {
	num_f = 0;
	newRandomPoint(x0);
	vector<double> x1(2), s(2), s_temp(2);
	double f0, f1;
	Rosenbrok(x0);
	f0 = Q(x0);
	num_f++;
	x1 = x0;
	for (iter = 0; iter < M; iter++) {
		newRandomPoint(s_temp);
		s_temp = 0.001 * s_temp;
		s = { 0,0 };
		int itr = 0;
		do {
			s += s_temp;
			x1 += s;
			f1 = Q(x1);
			num_f++;
			itr++;
			if (f1 < f0) {
				Rosenbrok(x1);
				f1 = Q(x1);
				num_f++;
				if (f1 < f0)
					x0 = x1;
				else
				x1 = x0;
			}
			f0 = f1;

		} while (f0 < f1 && itr < 1000);

		
	}
}

int main() {
	vector<double> x(2);
	random.seed(100);
	simpleRandom(x);
	cout << num_f << "\t\t" << x[0] << '\t' << x[1] <<'\t'<< -Q(x) << endl;
	//globalSearch1(x);
	//cout << num_f << "\t\t" << x[0] << '\t' << x[1] << '\t' << -Q(x) << endl;
	//globalSearch2_v1(x);
	//cout << num_f << "\t\t" << x[0] << '\t' << x[1] << '\t' << -Q(x) << endl;
	//globalSearch2_v2(x);
	//cout << num_f << "\t\t" << x[0] << '\t' << x[1] << '\t' << -Q(x) << endl;
	//globalSearch3(x);
	//cout << num_f << "\t\t" << x[0] << '\t' << x[1] << '\t' << -Q(x) << endl;
}