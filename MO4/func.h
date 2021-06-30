#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include <functional>
#include <random>
#include <iomanip>
#include <fstream>

using namespace std;

double r = 1., r1, r2;
uint64_t num_f, maxiter = 10000;
double eps = 1E-5; //точность

vector<double> operator * (const double& w, vector<double> vector) {
    size_t size = vector.size();
    for (size_t i = 0; i < size; ++i)
        vector[i] *= w;
    return vector;
}

vector<double> operator / (vector<double> vector, const double& w) {
    size_t size = vector.size();
    for (size_t i = 0; i < size; ++i)
        vector[i] /= w;
    return vector;
}

vector<double> operator + (vector<double> vector_1, const vector<double>& vector_2) {
    size_t size = vector_1.size();
    for (size_t i = 0; i < size; ++i)
        vector_1[i] += vector_2[i];
    return vector_1;
}

vector<double> operator - (vector<double> vector_1, const vector<double>& vector_2) {
    size_t size = vector_1.size();
    for (size_t i = 0; i < size; ++i)
        vector_1[i] -= vector_2[i];
    return vector_1;
}

vector<double>& operator += (vector<double>& vector_1, const vector<double>& vector_2) {
    size_t size = vector_1.size();
    for (size_t i = 0; i < size; ++i)
        vector_1[i] += vector_2[i];
    return vector_1;
}

double f(vector<double> x) {
	double sum = 0,
        C[6] = {6, 2, 4, 2, 8, 8},
		a[6] = {-3, 4, -8, -6, 3, -6},
		b[6] = {9, -7, 3, -9, -2, -8};

	for (int i = 0; i < 6; ++i)
		sum += C[i] / (1 + (x[0] - a[i]) * (x[0] - a[i]) + (x[1] - b[i]) * (x[1] - b[i]));
	return -sum;
    //return x[0] * x[0] + x[1] * x[1];
    //return (x[0] - 1) * (x[0] - 1) + (x[1] - 1) * (x[1] - 1);
}

double g(double x) {  // -10 <= x y <= 10
	return abs(x) - 10;
}

double G(double x) {
	return 0.5 * (g(x) + abs(g(x)));
}

inline void r_correct(vector<double> x) {
    if (r1 * G(x[0]) > eps)
        // Штрафные
       // r1 += .1;
        r1 += 1.0;
        //r1 *= 10.0;
        //r1 = pow(r1+1.0, 4.0);
    if (r2 * G(x[1]) > eps)
        // Штрафные
        //r2 += .1;
    r2 += 1.0;
    //r2 *= 10.0;
    //r2 = pow(r1+1.0, 4.0);
}

double Q (vector<double> x) {
    return f(x) + r * (r1 * G(x[0]) + r2 * G(x[1]));
}

double norma(const vector<double>& vector) {
    size_t size = vector.size();
    double res = 0;
    for (int i = 0; i < size; i++)
        res += vector[i] * vector[i];
    return sqrt(res);
}

int minInterval(vector<double> x, vector<double> p, double& a, double& b) {
    vector<double> x0(2);
    double fx0 = Q(x + a * p);
    double eps = 1E-8;
    if (fx0 < Q(x + (a + eps) * p))
        eps = -eps;

    double a1 = a + eps;
    double fx1 = Q(x + a1 * p);

    num_f += 3;
    do {
        eps *= 2;
        a = a1;
        a1 += eps;
        fx0 = fx1;
        fx1 = Q(x + a1 * p);
        ++num_f;
    } while (fx1 < fx0);

    if (a1 < a) {
        b = a;
        a = a1;
    }
    else
        b = a1;
    return num_f;
}

double F(int n) {
    return (pow(1.0 + sqrt(5.0), n) / pow(2.0, n) - pow(1.0 - sqrt(5.0), n) / pow(2.0, n)) / sqrt(5.0);
}

double fibonacci(vector<double> x, vector<double> p) {
    double alpha1 = 0;
    double alpha2 = 0;
    num_f += minInterval(x, p, alpha1, alpha2);
    //cout << alpha1 << " " << alpha2 << "\n";
    double a = alpha1;
    double b = alpha2;
    double aprev = 0;
    double bprev = 0;
    double lenght = b - a;
    int n = 0;
    while (F(n) < (b - a) / eps)
        n++;
    int iter = n - 3;
    double x1 = a + (F(n - 2) / F(n)) * (b - a);
    double f1 = Q(x + x1 * p);
    double x2 = a + (F(n - 1) / F(n)) * (b - a);
    double f2 = Q(x + x2 * p);
    num_f += 2;
    for (int i = 1; i <= iter; i++)
    {
        if (i == iter)//последняя итерация в фибоначи не требует вычисления f
        {
            if (f1 < f2)
            {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (F(n - i - 2) / F(n - i)) * (b - a);
            }
            else
            {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (F(n - i - 1) / F(n - i)) * (b - a);
            }
        }
        else {
            if (f1 < f2)
            {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (F(n - i - 2) / F(n - i)) * (b - a);
                f1 = Q(x + x1 * p);
                num_f += 1;
            }
            else
            {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (F(n - i - 1) / F(n - i)) * (b - a);
                f2 = Q(x + x2 * p);
                num_f += 1;
            }
        }
        lenght = b - a;
        bprev = b;
        aprev = a;
    }
    //cout << "num f" << num_f << "\n";
    //cout << a << "  " << b << "\n";
    return (a + b) / 2;
}

void Rosenbrok(vector<double>& x0) {
    int iter = 0;
    bool flag = true;
    double lambda1, lambda2;
    vector<double> x1(2), s1(2), s2(2), a1(2), a2(2);

    // Начальные ортогональные направления
    s1[0] = s2[1] = 1;
    s1[1] = s2[0] = 0;

    x1 = x0;
    //Розенброк (с ортогонализацией Палмера)
    while (flag == true && iter <= maxiter) {
        lambda1 = fibonacci(x1, s1);
        //cout << "lambda1=" << lambda1 << endl;
        x1 += lambda1 * s1;
        lambda2 = fibonacci(x1, s2);
        //lambda2 = 86.3636;
        //cout << "lambda2=" << lambda2 << endl;
        x1 += lambda2 * s2;
        a1 = lambda1 * s1 + lambda2 * s2;
        a2 = lambda2 * s2;
        /*
        cout << "s" << endl;
        cout << s1[0] << " " << s1[1] << " " << s2[0] << " " << s2[1] << endl;
        cout << "a" << endl;
        cout << a1[0] << " " << a1[1] << " " << a2[0] << " " << a2[1] << endl;
        */

        double a0_norm = norma(a1);
        double a1_norm = norma(a2);
        s1 = a1 / a0_norm;
        s2 = (a0_norm * a0_norm * a2 - a1_norm * a1_norm * a1) /
            (a0_norm * a1_norm * sqrt(a0_norm * a0_norm - a1_norm * a1_norm));
        iter++;
        r_correct(x1);
        //cout << num_f << '\t' << iter << '\t' << setprecision(15) << x1[0] << '\t' << setprecision(15) << x1[1] << '\t' << f(x1) << endl;

        flag = false;

        for (int i = 0; i < 2; i++)
            if (abs(x0[i] - x1[i]) > eps)
                flag = true;

        //cout << s1[0] << " " << s1[1] << " " << s2[0] << " " << s2[1] << endl;

        x0 = x1;
    }
}