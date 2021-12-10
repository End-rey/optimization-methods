#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <windows.h>
#include <ctime>

using namespace std;

double aa0, bb0, g1, g2;
int xx[5] = {0, 1, 2, 3, 4}, yy[5] = {1, 0, -2, 2, -3};
double EE = 10e-6;

double func1(double gamma)
{
    double f = 0;
    for (size_t i = 0; i < 5; i++)
    {
        f += pow(((aa0 + gamma * g1) * xx[i] + bb0 + gamma * g2 - yy[i]), 2);
    }
    return f;
}

double func2(double gamma)
{
    double f = 0;
    for (size_t i = 0; i < 5; i++)
    {
        f += abs(((aa0 + gamma * g1) * xx[i] + bb0 + gamma * g2 - yy[i]));
    }
    return f;
}

double func(double x, int fun)
{
    double y;
    if (fun == 1)
        y = func1(x);
    else
        y = func2(x);
    return y;
}

void output_oldbutgold()
{
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << "n"
         << "|" << setw(12) << left << "e"
         << "|" << setw(12) << left << "a"
         << "|" << setw(12) << left << "b"
         << "|" << setw(12) << left << "c"
         << "|" << setw(12) << left << "d"
         << "|" << setw(12) << left << "f(c)"
         << "|" << setw(12) << left << "f(d)"
         << "|" << endl;
}

void data_oldbutgold(int n, double e, double a, double b, double c, double d, double fc, double fd)
{
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << n
         << "|" << setw(12) << left << e
         << "|" << setw(12) << left << a
         << "|" << setw(12) << left << b
         << "|" << setw(12) << left << c
         << "|" << setw(12) << left << d
         << "|" << setw(12) << left << fc
         << "|" << setw(12) << left << fd
         << "|" << endl;
}

void output_dihit()
{
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << "n "
         << "|" << setw(12) << left << "a"
         << "|" << setw(12) << left << "b"
         << "|" << setw(12) << "e = (b-a)/2"
         << "|" << setw(12) << left << "c"
         << "|" << setw(12) << left << "d"
         << "|" << setw(12) << left << "f(c)"
         << "|" << setw(12) << left << "f(d)"
         << "|" << endl;
}

void data_dihit(int n, double a, double b, double e, double c, double d, double fc, double fd)
{
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << n
         << "|" << setw(12) << left << a
         << "|" << setw(12) << left << b
         << "|" << setw(12) << left << e
         << "|" << setw(12) << left << c
         << "|" << setw(12) << left << d
         << "|" << setw(12) << left << fc
         << "|" << setw(12) << left << fd
         << "|" << endl;
}

void dihit(int a, int b, int fun)
{
    double z = EE, a0 = a, b0 = b, step = 0;
    double c = (a0 + b0 - z) / 2, d = (a0 + b0 + z) / 2, fc = func(c, fun), fd = func(d, fun), e = (b0 - a0) / 2;
    data_dihit(step, a0, b0, e, c, d, fc, fd);
    while (e > EE)
    {
        step++;
        if (fc <= fd)
        {
            b0 = d;
        }
        else
        {
            a0 = c;
        }
        c = (a0 + b0 - z) / 2;
        d = (a0 + b0 + z) / 2;
        fc = func(c, fun);
        fd = func(d, fun);
        e = (b0 - a0) / 2;
        data_dihit(step, a0, b0, e, c, d, fc, fd);
    }
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "Answer: x = " << (a0 + b0) / 2 << " "
         << "y = " << func((a0 + b0) / 2, fun) << endl
         << endl;
}

void oldbutgold(int a, int b, int fun)
{
    float min = 0, x, y, k;
    double step = -1, delta, alpha, beta, fal, fbe, e = b - a;
    double a0 = a, b0 = b, a1 = a, b1 = b;
    while (e > EE)
    {
        step++;
        delta = (b0 - a0);
        alpha = a0 + 2 * delta / (3 + sqrt(5));
        beta = a0 + 2 * delta / (1 + sqrt(5));
        fal = func(alpha, fun);
        fbe = func(beta, fun);
        e *= (sqrt(5) - 1) / 2;
        data_oldbutgold(step, e, a0, b0, alpha, beta, fal, fbe);
        if (fal <= fbe)
        {
            a1 = a0;
            b1 = beta;
        }
        else
        {
            a1 = alpha;
            b1 = b0;
        }
        a0 = a1;
        b0 = b1;
    }
    k = (a1 + b1) / 2;
    min = func(k, fun);
    data_oldbutgold(step + 1, e, a0, b0, alpha, beta, fal, fbe);
    cout << "_________________________________________________________________________________________________" << endl;
    cout << "Answer: x = " << k << " "
         << "y = " << min << endl
         << endl;
}

int main()
{
    srand((unsigned int)time(NULL));
    aa0 = rand() % 5;
    bb0 = rand() % 5;
    g1 = (double)rand() / RAND_MAX;
    g2 = sqrt(1 - g1 * g1);
    cout << aa0 << " " << bb0 << endl
         << g1 << " " << g2 << endl
         << endl;
    cout << "1st function" << endl
         << endl;
    cout << "Gold Ratio" << endl;
    output_oldbutgold();
    oldbutgold(-10, 10, 1);
    cout << "Dichotomy" << endl;
    output_dihit();
    dihit(-10, 10, 1);
    cout << endl
         << "2nd function" << endl
         << endl;
    cout << "Gold Ratio" << endl;
    output_oldbutgold();
    oldbutgold(-10, 10, 2);
    cout << "Dichotomy" << endl;
    output_dihit();
    dihit(-10, 10, 2);
    system("pause");
    return 0;
}