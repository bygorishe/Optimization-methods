#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>

const double x_begin = -2.0, x_end = 20.0;

inline double f(double x)
{
    return (x - 3.0) * (x - 3.0);
}

void dichotomy(FILE* file, double eps)
{
    double a = x_begin;
    double b = x_end;
    double aprev = 0;
    double bprev = 0;
    double delta = eps / 2.0;
    int iter = 0;
    int callf = 0;
    double x;
    double lenght = b - a;;
    fprintf(file, "eps = %.1E\n", eps);
    fprintf(file, "Итерация x1                    x2                    f1                    f2                    a                     b                     Длина отрезка        Отношение длин\n");
    while (lenght/2 > eps)
    {
        iter++;
        double x1 = (a + b - delta) / 2.0;
        double x2 = (a + b + delta) / 2.0;
        double f1 = f(x1);
        double f2 = f(x2);
        callf += 2;
        if (f1 < f2)
            b = x2;
        else
            a = x1;
        lenght = b - a;
        fprintf(file, "%-8i %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E \n", iter, x1, x2, f1, f2, a, b, lenght, (bprev - aprev) / (b - a));
        bprev = b;
        aprev = a;
    }
    fprintf(file, "точка: %-18.15E\n\n", (a + b) / 2);
    fprintf(file, "количетсво вычислений функции: %-3i\n\n", callf);
}

void golden_ratio(FILE* file, double eps)
{
    double a = x_begin;
    double b = x_end;
    double aprev = 0;
    double bprev = 0;
    double q1 = (3.0 - sqrt(5.0)) / 2.0;
    double q2 = (sqrt(5.0) - 1.0) / 2.0;
    int iter = 0;
    double lenght = b - a;
    double x1 = a + q1 * (b - a);
    double x2 = a + q2 * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);
    int callf = 2;
    int flag;
    fprintf(file, "eps = %.1E\n", eps);
    fprintf(file, "Итерация x1                    x2                    f1                    f2                    a                     b                     Длина отрезка        Отношение длин\n");
    while (lenght/2 > eps)
    {
        iter++;
        if (f1 < f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + q1 * (b - a);
            f1 = f(x1);
            callf += 1;
        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + q2 * (b - a);
            f2 = f(x2);
            callf += 1;
        }
        lenght = b - a;
        fprintf(file, "%-8i %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E \n", iter, x1, x2, f1, f2, a, b, lenght, (bprev - aprev) / (b - a));
        bprev = b;
        aprev = a;
    }
    fprintf(file, "точка: %-18.15E\n\n", (a + b) / 2);
    fprintf(file, "количетсво вычислений функции: %-3i\n\n", callf);
}

double F(int n)
{
    return (pow(1.0 + sqrt(5.0), n) / pow(2.0, n) - pow(1.0 - sqrt(5.0), n) / pow(2.0, n)) / sqrt(5.0);
}

void fibonacci(FILE* file, double eps)
{
    double a = x_begin;
    double b = x_end;
    double aprev = 0;
    double bprev = 0;
    double lenght = b - a;
    fprintf(file, "eps = %.1E\n", eps);
    int n = 0;
    while (F(n) < (b - a) / eps) n++;
    int iter = n - 3;
    double x1 = a + (F(n - 2) / F(n)) * (b - a);
    double f1 = f(x1);
    double x2 = a + (F(n - 1) / F(n)) * (b - a);
    double f2 = f(x2);
    int callf = 2;
    fprintf(file, "Итерация x1                    x2                    f1                    f2                    a                     b                     Длина отрезка        Отношение длин\n");
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
                f1 = f(x1);
                callf += 1;
            }
            else
            {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (F(n - i - 1) / F(n - i)) * (b - a);
                f2 = f(x2);
                callf += 1;
            }
        }
        lenght = b - a;
        fprintf(file, "%-8i %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E %-18.15E \n", i, x1, x2, f1, f2, a, b, lenght, (bprev - aprev) / (b - a));
        bprev = b;
        aprev = a;
    }
    fprintf(file, "точка: %-18.15E\n\n", (a + b) / 2);
    fprintf(file, "количетсво вычислений функции: %-3i\n\n", callf);
}

void interval(FILE* file, double& a, double& b, double x0, double eps)
{
    double xk0 = x0;
    double fxk0 = f(xk0);
    double xk;
    double fxk;
    double xk1;
    double fxk1;
    double h = 0;
    int iter = 0;
    if (f(x0) > f(x0 + eps))
    {
        xk = x0 + eps;
        h += eps;
    }
    else if (f(x0) > f(x0 - eps))
    {
        xk = x0 - eps;
        h -= eps;
    }
    else//если мы уже в окрестностях точки
    {
        a = x0 - eps;
        b = x0 + eps;
        fprintf(file, "%-7.1E [%-18.15E, %-18.15E] %-18.15E %-8i\n", eps, a, b, b - a, iter);
        return;
    }
    fxk = f(xk);
    h *= 2.0;
    xk1 = xk + h;
    fxk1 = f(xk1);
    while (fxk > fxk1)
    {
        xk0 = xk;
        fxk0 = fxk;
        xk = xk1;
        fxk = fxk1;
        iter++;

        h *= 2.0;
        xk1 = xk + h;
        fxk1 = f(xk1);
    }
    a = xk0;
    b = xk1;
    fprintf(file, "%-7.1E [%-18.15E, %-18.15E] %-18.15E %-8i\n", eps, a, b, b - a, iter);
}

int main()
{
    double x;
    FILE* f;
    f = fopen("dichotomy.txt", "w");
    for (double eps = 1.0E-1; eps > 1.0E-7; eps *= 1.0E-1)
    {
        dichotomy(f, eps); 
    }
    fclose(f);
    f = fopen("golden_ratio.txt", "w");
    for (double eps = 1.0E-1; eps > 1.0E-7; eps *= 1.0E-1)
    {
        golden_ratio(f, eps);
    }
    fclose(f);
    f = fopen("fibonacci.txt", "w");
    for (double eps = 1.0E-1; eps > 1.0E-7; eps *= 1.0E-1)
    {
        fibonacci(f, eps);
    }
    fclose(f);
    
    double a, b;
    f = fopen("interval1.txt.txt", "w");
    fprintf(f, "Eps     Interval                                       Length                Iterations\n");
    for (double eps = 1.0E-1; eps > 1.0E-7; eps *= 1.0E-1)
    {
        interval(f, a, b, 1.0, eps);
    }
    fclose(f);
    return 0;
}