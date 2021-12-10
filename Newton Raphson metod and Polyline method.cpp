#include <iostream>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <vector>
#include <climits>

using namespace std;

double aa0, bb0, g1, g2;
int xx[5] = {0, 1, 2, 3, 4}, yy[5] = {1, 0, -2, 2, -3};
double EE = 10e-4;
double a = -10, b = 10;

double func(double x, int fun);
double deff1(double x, int fun);
double dabsfunc(double x, int fun);
double oldbutgold(int a, int b, int fun);

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

double func3(double gamma)
{
    double f = 0;
    for (size_t i = 0; i < 3; i++)
    {
        if (abs(((aa0 + gamma * g1) * xx[i] + bb0 + gamma * g2 - yy[i])) > f)
            f = abs(((aa0 + gamma * g1) * xx[i] + bb0 + gamma * g2 - yy[i]));
    }
    return f;
}

double func(double x, int fun)
{
    double y;
    if (fun == 1)
        y = func1(x);
    else if (fun == 2)
        y = func2(x);
    else if (fun == 3)
        y = func3(x);
    else if (fun == 4)
        y = sin(x) / x;
    else
        y = dabsfunc(x, fun);
    return y;
}

double deff1(double x, int fun)
{
    double h = 10e-6;
    return ((func(x + h, fun) - func(x, fun)) / (h));
}

double dabsfunc(double x, int fun)
{
    double y;
    if (fun == 21)
        y = fabs(deff1(x, 1));
    else if (fun == 22)
        y = fabs(deff1(x, 2));
    else if (fun == 23)
        y = fabs(deff1(x, 3));
    else if (fun == 24)
        y = fabs(deff1(x, 4));
    return y;
}

double deff2(double x, int fun)
{
    double h = 10e-6;
    return ((deff1(x + h, fun) - deff1(x, fun)) / (h));
}

void output_newraf()
{
    cout << "_____________________________________________" << endl;
    cout << "|" << setw(4) << left << "n"
         << "|" << setw(12) << left << "xn"
         << "|" << setw(12) << left << "f(xn)"
         << "|" << setw(12) << left << "f'(xn)"
         << "|" << endl;
}

void data_newraf(int n, double x, double y, double y1)
{
    cout << "_____________________________________________" << endl;
    cout << "|" << setw(4) << left << n
         << "|" << setw(12) << left << x
         << "|" << setw(12) << left << y
         << "|" << setw(12) << left << y1
         << "|" << endl;
}

double newraf(int a, int b, int fun)
{
    double xk = b;
    int n = 0;
    while (abs(deff1(xk, fun)) > EE)
    {
        //if (n % a == 0)
        data_newraf(n, xk, deff1(xk, fun), deff2(xk, fun));
        if (deff2(xk, fun) != 0)
            xk -= (deff1(xk, fun) / deff2(xk, fun));
        else
            break;
        n += 1;
    }
    data_newraf(n, xk, deff1(xk, fun), deff2(xk, fun));
    cout << "_____________________________________________" << endl
         << "Answer: x = " << xk << ", y = " << func(xk, fun) << endl
         << endl;
    return func(xk, fun);
}

double lipsh(int fun)
{
    double L = (oldbutgold(a, b, fun));
    if (L > 0)
        return round(L * 1000) / 1000;
    else
        return 0;
}

double oldbutgold(int a, int b, int fun) 
{
    float max = 0, x, y, k;
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
        if (fal >= fbe)
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
    max = func(k, fun);
    return max;
}

void dihit(int a, int b, int fun)
{
    double z = EE, a0 = a, b0 = b, step = 0;
    double c = (a0 + b0 - z) / 2, d = (a0 + b0 + z) / 2, fc = func(c, fun), fd = func(d, fun), e = (b0 - a0) / 2;
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
    }
    cout << "Answer: x = " << (a0 + b0) / 2 << " "
         << "y = " << func((a0 + b0) / 2, fun) << endl
         << endl;
}

void output_polyline()
{
    cout << "____________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << left << " "
         << "|" << setw(25) << left << "Excluded pair (x, p)"
         << "|" << setw(12) << left << " "
         << "|" << setw(38) << left << "Included pair (x, p)"
         << "|" << endl;
    cout << "____________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << left << "n"
         << "|" << setw(12) << left << "xn*"
         << "|" << setw(12) << left << "pn**"
         << "|" << setw(12) << left << "2L^n"
         << "|" << setw(12) << left << "xn'"
         << "|" << setw(12) << left << "xn''"
         << "|" << setw(12) << left << "pn"
         << "|" << endl;
}

void data_polyline(int n, double xn1, double pn1, double l, double xn2, double xn3, double pn2)
{
    cout << "____________________________________________________________________________________" << endl;
    cout << "|" << setw(4) << left << n
         << "|" << setw(12) << left << xn1
         << "|" << setw(12) << left << pn1
         << "|" << setw(12) << left << l
         << "|" << setw(12) << left << xn2
         << "|" << setw(12) << left << xn3
         << "|" << setw(12) << left << pn2
         << "|" << endl;
}

void polyline(int a, int b, int fun)
{
    int n = 1, m;
    double L = lipsh(20 + fun), x1, p1, x2, x3, p2, del;
    vector<double> x, p;
    x1 = (func(a, fun) - func(b, fun) + L * (a + b)) / (2 * L);
    p1 = (func(a, fun) + func(b, fun)) / 2 + (L * (a - b)) / 2;
    del = (func(x1, fun) - p1) / (2 * L);
    x2 = x1 - del;
    x3 = x1 + del;
    p2 = (func(x1, fun) + p1) / 2;
    data_polyline(n, x1, p1, 2 * del * L, x2, x3, p2);
    x.push_back(x2);
    x.push_back(x3);
    p.push_back(p2);
    p.push_back(p2);
    n++;
    x1 = x2;
    p1 = p2;
    del = (func(x1, fun) - p1) / (2 * L);
    x2 = x1 + del;
    x3 = x1 - del;
    p2 = (func(x1, fun) + p1) / 2;
    data_polyline(n, x1, p1, 2 * del * L, x2, x3, p2);
    x[0] = x2;
    x.push_back(x3);
    p[0] = p2;
    p.push_back(p2);
    n++;
    while (2 * del * L > EE)
    {
        double gg = INT_MAX;
        for (size_t i = 0; i < p.size(); i++)
        {
            if (p[i] < gg)
            {
                gg = p[i];
                m = i;
            }
        }
        x1 = x[m];
        p1 = p[m];
        del = (func(x1, fun) - p1) / (2 * L);
        x2 = x1 - del;
        x3 = x1 + del;
        p2 = (func(x1, fun) + p1) / 2;
        if (n % 10000 == 0)
            data_polyline(n, x1, p1, 2 * del * L, x2, x3, p2);
        x[m] = x2;
        x.push_back(x3);
        p[m] = p2;
        p.push_back(p2);
        n++;
    };
    data_polyline(n - 1, x1, p1, 2 * del * L, x2, x3, p2);
    cout << "____________________________________________________________________________________" << endl;
    cout << "Answer x = " << x1 << ", y = " << func(x1, fun) << endl;
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
    cout << "Newton Raphson metod for the 1st function" << endl;
    output_newraf();
    double ans = newraf(a, b, 1);
    cout << "Const lipshutz\n For the 1st function\n"
         << lipsh(21) << "\n For the 2nd function\n"
         << lipsh(22) << "\n For the 3rd function\n"
         << lipsh(23) << endl
         << "\n Polyline method\n For the 1st function\n";
    output_polyline();
    polyline(a, b, 1);
    cout << "\n For the 2nd function\n";
    polyline(a, b, 2);
    cout << "\n For the 3rd function\n";
    polyline(a, b, 3);
    system("pause");
    return 0;
}
