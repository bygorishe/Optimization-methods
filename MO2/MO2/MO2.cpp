#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

const double PI = 3.14159265;

double A1 = 1;
double A2 = 2;
double a1 = 3;
double a2 = 2;
double b1 = 1;
double b2 = 2;
double c1 = 1;
double c2 = 2;
double d1 = 3;
double d2 = 1;

vector<double> operator * (const double& w, vector<double> vector) {
    size_t size = vector.size();
    for (size_t i = 0; i < size; ++i)
        vector[i] *= w;
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

double norma(const vector<double>& vector) {
    size_t size = vector.size();
    double res = 0;
    for (int i = 0; i < size; i++)
        res += vector[i] * vector[i];
    return sqrt(res);
}

double f1(vector<double> x) { //из варианта
    return -A1 * exp(-((x[0] - a1) / b1) * ((x[0] - a1) / b1) - ((x[1] - c1) / d1) * ((x[1] - c1) / d1))
        - A2 * exp(-((x[0] - a2) / b2) * ((x[0] - a2) / b2) - ((x[1] - c2) / d2) * ((x[1] - c2) / d2));
}

void grad_f1(vector<double> x, vector<double>& vector_grad) { 
    vector_grad[0] = 2 * (x[0] - 3) * exp(-(x[0] - 3) * (x[0] - 3) - 1 / 9 * (x[1] - 1) * (x[1] - 1)) + ((x[0] - 2) * exp((-x[0] * x[0]) / 4 + x[0] - x[1] * x[1] + 4 * x[1] - 5));
    vector_grad[1] = 2 / 9 * (x[1] - 1) * exp(-(x[0] - 3) * (x[0] - 3) - 1 / 9 * (x[1] - 1) * (x[1] - 1)) + 4 * ((x[1] - 2) * exp((-x[0] * x[0]) / 4 + x[0] - x[1] * x[1] + 4 * x[1] - 5));
}

double f2(vector<double> x) { //квадратичная функция 
    return 100 * (x[1] - x[0]) * (x[1] - x[0]) + (1 - x[0]) * (1 - x[0]);
}

void grad_f2(vector<double> x, vector<double>& vector_grad) {
    vector_grad[0] = 202 * x[0] - 2 - 200 * x[1];
    vector_grad[1] = -200 * x[0] + 200 * x[1];
}

double f3(vector<double> x) { //функция розенброка
    return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
}

void grad_f3(vector<double> x, vector<double>& vector_grad) {
    vector_grad[0] = 400 * x[0] * x[0] * x[0] - 400 * x[0] * x[1] + 2 * x[0] - 2;
    vector_grad[1] = -200 * x[0] * x[0] + 200 * x[1];
}

int minInterval(double (*f)(vector<double>), vector<double> x, vector<double> p, double& a, double& b, double epsilon) {
    vector<double> x0(2);
    double fx0 = f(x + a * p);
    //double eps = epsilon;
    double eps = 1E-8;
    if (fx0 < f(x + (a + eps) * p))
        eps = -eps;

    double a1 = a + eps;
    double fx1 = f(x + a1 * p);

    int num_f = 3;
    do {
        eps *= 2;
        a = a1;
        a1 += eps;
        fx0 = fx1;
        fx1 = f(x + a1 * p);
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

double fibonacci(double (*f)(vector<double>), vector<double> x, vector<double> p, double alpha1, double ep, int& num_f) {
    double alpha2 = 0;
    double eps = 1E-8;
    num_f += minInterval(f, x, p, alpha1, alpha2, eps);
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
    double f1 = f(x + x1 * p);
    double x2 = a + (F(n - 1) / F(n)) * (b - a);
    double f2 = f(x + x2 * p);
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
                f1 = f(x + x1 * p);
                num_f += 1;
            }
            else
            {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (F(n - i - 1) / F(n - i)) * (b - a);
                f2 = f(x + x2 * p);
                num_f += 1;
            }
        }
        lenght = b - a;
        bprev = b;
        aprev = a;
    }
    return (a + b) / 2;
}

void BroydenMethod(double (*f)(vector<double>), void (*grad)(vector<double>, vector<double>&), vector<double>& x0, double eps, ofstream& fout) {
    int num_f = 0;
    vector<double> vector_grad(2), p(2), x1(2), gk(2), sk(2);
    double lambda, pk; //р - точка в направлении которой будем производить поиск
                        //sk = x1-x0
                        //gk = grad x1 - grad x0
    //гексиан Н
    vector<vector<double>> H;
    H.resize(2);
    for (int i = 0; i < 2; i++)
        H[i].resize(2);
    vector<vector<double>> Htemp;
    Htemp.resize(2);
    for (int i = 0; i < 2; i++)
        Htemp[i].resize(2);
    vector<vector<double>> Hk;
    Hk.resize(2);
    for (int i = 0; i < 2; i++)
        Hk[i].resize(2);

    //изначально единичная матрица
    H[0][1] = H[1][0] = 0;
    H[0][0] = H[1][1] = 1;

    int iter = 0;

    fout << setw(5) << "iter";
    fout << setw(14) << "x" << setw(14) << "y" << setw(14) << "f(x, y)";
    fout << setw(14) << "Sx" << setw(14) << "Sy" << setw(14) << "lambda";
    fout << setw(14) << "|xk - x(k-1)|" << setw(14) << "|yk - y(k-1)|" << setw(14) << "|fk - f(k-1)|";
    fout << setw(14) << "angle" << setw(28) << "gradient" << setw(56) << "H" << endl;

    grad(x0, vector_grad); //считаем градиент

    bool flag = true;
    while (norma(vector_grad) > eps && iter < 100 && flag) {
        p[0] = -H[0][0] * vector_grad[0] - H[0][1] * vector_grad[1];
        p[1] = -H[1][0] * vector_grad[0] - H[1][1] * vector_grad[1];

        lambda = fibonacci(f, x0, p, 0, eps, num_f);

        x1 = x0 + lambda * p;//следующая точка
        sk = x1 - x0;

        gk = (-1) * vector_grad;
        //считаем новый градиент
        grad(x1, vector_grad);
        gk += vector_grad;

        //новое приближение гексиана
        //H_k+1 = (I - pk * sk * gk^T)H_k(I - pk * gk * sk^T) + pk * sk * sk^T
        //pk = 1 / gk^T * sk
        pk = 1 / (gk[0] * sk[0] + gk[1] * sk[1]);

        Htemp[0][0] = 1 - pk * gk[0] * sk[0];
        Htemp[0][1] = -pk * gk[1] * sk[0];
        Htemp[1][0] = -pk * gk[0] * sk[1];
        Htemp[1][1] = 1 - pk * gk[1] * sk[1];

        Hk[0][0] = Htemp[0][0] * H[0][0] + Htemp[0][1] * H[1][0];
        Hk[0][1] = Htemp[0][0] * H[0][1] + Htemp[0][1] * H[1][1];
        Hk[1][0] = Htemp[1][0] * H[0][0] + Htemp[1][1] * H[1][0];
        Hk[1][1] = Htemp[1][0] * H[0][1] + Htemp[1][1] * H[1][1];

        Htemp[0][0] = 1 - pk * gk[0] * sk[0];
        Htemp[1][0] = -pk * gk[1] * sk[0];
        Htemp[0][1] = -pk * gk[0] * sk[1];
        Htemp[1][1] = 1 - pk * gk[1] * sk[1];

        H[0][0] = Hk[0][0] * Htemp[0][0] + Hk[0][1] * Htemp[1][0] + pk * sk[0] * sk[0];
        H[0][1] = Hk[0][0] * Htemp[0][1] + Hk[0][1] * Htemp[1][1] + pk * sk[0] * sk[1];
        H[1][0] = Hk[1][0] * Htemp[0][0] + Hk[1][1] * Htemp[1][0] + pk * sk[1] * sk[0];
        H[1][1] = Hk[1][0] * Htemp[0][1] + Hk[1][1] * Htemp[1][1] + pk * sk[1] * sk[1];

        fout << setw(5) << iter + 1;
        fout << setw(14) << setprecision(8) << x1[0] << setw(14) << setprecision(8) << x1[1] << setw(14) << f(x1);
        fout << setw(14) << p[0] << setw(14) << p[1] << setw(14) << lambda;
        fout << setw(14) << abs(sk[0]) << setw(14) <<abs(sk[1]) << setw(14) << f(x1) - f(x0);
        fout << setw(14) << acos((x1[0] * p[0] + x1[1] * p[1]) / (norma(x1) * norma(p))) * 180 / PI << setw(28) << vector_grad[0] << " " << vector_grad[1] << setw(56) << H[0][0] << " " << H[0][1] << H[1][0] << " " << H[1][1] << endl;

        flag = false;

        for (int i = 0; i < 2; i++)
            if (abs(x0[i] - x1[i]) > eps)
                flag = true;

        x0 = x1;
        iter++;
    }
    fout << "call_func: " << num_f;
}

void ConjugateGradientMethod(double (*f)(vector<double>), void (*grad)(vector<double>, vector<double>&), vector<double>& x0, double eps, ofstream& fout) {
    int num_f = 0;
    double lambda, w; //р - точка в направлении которой будем производить поиск
    vector<double> vector_grad(2), vector_grad_old(2), x1(2), sk(2);
    //sk направление поиска

    grad(x0, vector_grad); //считаем градиент

    //Начальное направление поиска
    sk = (-1) * vector_grad;

    int iter = 0;

    fout << setw(5) << "iter";
    fout << setw(14) << "x" << setw(14) << "y" << setw(14) << "f(x, y)";
    fout << setw(14) << "Sx" << setw(14) << "Sy" << setw(14) << "lambda";
    fout << setw(14) << "|xk - x(k-1)|" << setw(14) << "|yk - y(k-1)|" << setw(14) << "|fk - f(k-1)|";
    fout << setw(14) << "angle" << setw(28) << "gradient" << setw(56) << "H" << endl;

    bool flag = true;
    while (norma(sk) > eps && iter < 100 && flag) {
        lambda = fibonacci(f, x0, sk, 0, eps, num_f);
        x1 = x0 + lambda * sk;  //следующая точка

        //заопниминаем старый градиент
        vector_grad_old = vector_grad;

        //считаем новый градиент
        grad(x1, vector_grad);

        //модификация полака-рибьера (пшениного)
        w = (vector_grad[0] * (-vector_grad[0] + vector_grad_old[0])
            + vector_grad[1] * (-vector_grad[1] + vector_grad_old[1]))
            / (vector_grad_old[0] * sk[0] + vector_grad_old[1] * sk[1]);

        //флетчер ривс
        /*w = (vector_grad[0] * vector_grad[0]
             + vector_grad[1] * vector_grad[1])
             / (vector_grad_old[0] * vector_grad_old[0]
                 + vector_grad_old[1] * vector_grad_old[1]);*/

        sk = (-1) * vector_grad + w * sk;

        fout << setw(5) << iter + 1;
        fout << setw(14) << setprecision(8) << x1[0] << setw(14) << setprecision(8) << x1[1] << setw(14) << f(x1);
        fout << setw(14) << sk[0] << setw(14) << sk[1] << setw(14) << lambda;
        fout << setw(14) << abs(x1[0] - x0[0]) << setw(14) << abs(x1[1] - x0[1]) << setw(14) << f(x1) - f(x0);
        fout << setw(14) << acos((x1[0] * sk[0] + x1[1] * sk[1]) / (norma(x1) * norma(sk))) * 180 / PI << setw(28) << vector_grad[0] << " " << vector_grad[1] << endl;

        flag = false;

        for (int i = 0; i < 2; i++)
            if (abs(x0[i] - x1[i]) > eps)
                flag = true;

        x0 = x1;

        iter++;
    }
    fout << "call_func: " << num_f << endl;
    cout << num_f << endl;
}

void main() {
    //double eps = 1E-7;
    vector<double> x0(2);
    x0[0] = x0[1] = 0;

    ofstream fout1("results_broyden.txt");
    for (double eps = 1E-3; eps >= 1E-7; eps *= 0.1) {
        fout1 << "eps: " << eps << endl;
        BroydenMethod(f1, grad_f1, x0, eps, fout1); //https://habr.com/ru/post/333356/
        x0[0] = x0[1] = 0;
        fout1 << endl;
    }

    ofstream fout2("results_gradient_method.txt");
    for (double eps = 1E-3; eps >= 1E-7; eps *= 0.1) {
        fout2 << "eps: " << eps << endl;
        ConjugateGradientMethod(f1, grad_f1, x0, eps, fout2);
        x0[0] = x0[1] = 0;
        fout2 << endl;
    }
}