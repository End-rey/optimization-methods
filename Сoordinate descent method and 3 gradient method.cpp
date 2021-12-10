#include <iostream>
#include <math.h>
#include <iomanip>
#include <ctime>

using namespace std; 

int xx[5] = {0, 1, 2, 3, 4}, yy[5] = {1, 0, -2, 2, -3};
const double EE = 10e-6;
// double a = -10, b = 10;

double dfx(double x, double y, int fun);
double dfy(double x, double y, int fun);

double func1(double x, double y)
{
    double f = 0;
    for (size_t i = 0; i < 5; i++)
    {
        f += pow((x * xx[i] + y - yy[i]), 2);
    }
    return f;
}

double func2(double x, double y)
{
    double f = 0;
    for (size_t i = 0; i < 5; i++)
    {
        f += abs(x * xx[i] + y - yy[i]);
    }
    return f;
}

double func(double x, double y, int fun)
{
    double f;
    if (fun == 1)
        f = func1(x, y);
    else if (fun == 2)
        f = func2(x, y);
    return f;
}

double funxy(double x0, double y0, double al, int fun, int xy)
{
    return ((xy == 1) ? func(x0 + al, y0, fun) : (xy == 2) ? func(x0, y0 + al, fun)
                                                           : func(x0 - al * dfx(x0, y0, fun), y0 - al * dfy(x0, y0, fun), fun));
}

double dfx(double x, double y, int fun)
{
    return (func(x + EE, y, fun) - func(x, y, fun)) / EE;
}

double dfy(double x, double y, int fun)
{
    return (func(x, EE + y, fun) - func(x, y, fun)) / EE;
}

void data_coordinate_descent(int n, double xn, double yn, double mod, double fxn, double modf)
{
    cout << "_______________________________________________________________________________" << endl;
    cout << "|" << setw(4) << left << n
         << "|" << setw(12) << left << xn
         << "|" << setw(12) << left << yn
         << "|" << setw(14) << left << mod
         << "|" << setw(12) << left << fxn
         << "|" << setw(18) << left << modf
         << "|" << endl;
}

void output_coordinate_descent()
{
    cout << "_______________________________________________________________________________" << endl;
    cout << "|" << setw(4) << left << "n"
         << "|" << setw(12) << left << "xn"
         << "|" << setw(12) << left << "yn"
         << "|" << setw(14) << left << "||xn - xn-1||"
         << "|" << setw(12) << left << "f(xn)"
         << "|" << setw(18) << left << "|f(xn) - f(xn-1)|"
         << "|" << endl;
}

void output_gradient()
{
    cout << "__________________________________________________________" << endl;
    cout << "|" << setw(4) << left << "n"
         << "|" << setw(12) << left << "xn"
         << "|" << setw(12) << left << "yn"
         << "|" << setw(12) << left << "f(xn)"
         << "|" << setw(12) << left << "||f'(xn)||"
         << "|" << endl;
}

void data_gradient(int n, double x, double y, double f, double df)
{
    cout << "__________________________________________________________" << endl;
    cout << "|" << setw(4) << left << n
         << "|" << setw(12) << left << x
         << "|" << setw(12) << left << y
         << "|" << setw(12) << left << f
         << "|" << setw(12) << left << df
         << "|" << endl;
}

double dihit(double x0, double y0, int a, int b, int fun, int xy)
{

    double z = EE, a0 = a, b0 = b, step = 0;
    double c = (a0 + b0 - z) / 2, d = (a0 + b0 + z) / 2, fc = funxy(x0, y0, c, fun, xy), fd = funxy(x0, y0, d, fun, xy), e = (b0 - a0) / 2;

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
        fc = funxy(x0, y0, c, fun, xy);
        fd = funxy(x0, y0, d, fun, xy);
        e = (b0 - a0) / 2;
    }
    return (a0 + b0) / 2;
}

void coordinate_descent(int fun)
{
    output_coordinate_descent();
    srand(time(0));
    int n = 0;
    double all, x1 = -10 + rand() % 20, y1 = -10 + rand() % 20, alx = 100, aly = 100, x0 = x1, y0 = y1, al = 100;
    do
    {
        n++;
        if (n % 10 == 0 or n == 1)
            data_coordinate_descent(n, x1, y1, sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)), func(x1, y1, fun), abs(func(x0 + alx, y0 + alx, fun) - func(x0 + aly, y0 + aly, fun)));
        x0 = x1;
        y0 = y1;
        al = 100;
        all = 0;
        alx = dihit(x0, y0, all, al, fun, 1);
        if (abs(alx) < EE)
        {
            alx = dihit(x0, y0, -al, -all, fun, 1);
            if (abs(alx) < EE)
            {
                aly = dihit(x0, y0, all, al, fun, 2);
                if (abs(aly) < EE)
                {
                    aly = dihit(x0, y0, -al, -all, fun, 2);
                    while (abs(aly - al) < EE)
                    {
                        al -= al;
                        all -= al;
                        aly = dihit(x0, y0, -al, -all, fun, 2);
                    }
                    y1 = y0 + aly;
                }
                else
                {
                    while (abs(aly - al) < EE)
                    {
                        al += al;
                        all += al;
                        aly = dihit(x0, y0, all, al, fun, 2);
                    }
                    y1 = y0 + aly;
                }
            }
            else
                while (abs(alx - al) < EE)
                {
                    al -= al;
                    all -= al;
                    alx = dihit(x0, y0, -al, -all, fun, 1);
                }
            x1 = x0 + alx;
        }
        else
        {
            while (abs(alx - al) < EE)
            {
                al += al;
                all += al;
                alx = dihit(x0, y0, all, al, fun, 1);
            }
            x1 = x0 + alx;
        }
    } while ((abs(func(x0, y0, fun) - func(x1, y1, fun)) > EE) and (sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2))) > EE);
    data_coordinate_descent(n, x1, y1, sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)), func(x1, y1, fun), abs(func(x0, y0, fun) - func(x1, y1, fun)));
    cout << "_______________________________________________________________________________" << endl;
    cout << "Answer: x = " << x1 << " , y = " << y1 << ", f(x,y) = " << func(x1, y1, fun) << endl
         << endl;
}

void grad_split(int fun)
{
    srand(time(0) * 2);
    output_gradient();
    int n = 0;
    double x1 = -10 + rand() % 20, y1 = -10 + rand() % 20, al = 10, x0 = x1, y0 = y1;
    do
    {
        n++;
        if (n % 50 == 0 or n == 1)
            data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
        x0 = x1;
        y0 = y1;
        if ((func(x0 - al * dfx(x0, y0, fun), y0 - al * dfy(x0, y0, fun), fun) - func(x0, y0, fun)) <= (-al * EE * (pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2))))
        {
            x1 = x0 - al * dfx(x0, y0, fun);
            y1 = y0 - al * dfy(x0, y0, fun);
        }
        else
            al = al / 2;

    } while (((abs(func(x0, y0, fun) - func(x1, y1, fun)) > EE) and (sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2))) > EE) or (sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)) > EE and al > EE));
    data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
    cout << "__________________________________________________________" << endl;
    cout << "Answer: x = " << x1 << ", y = " << y1 << ", f(x,y) = " << func(x1, y1, fun) << endl
         << endl;
}

void grad_const(int fun)
{
    srand(time(0) * 3);
    output_gradient();
    int n = 0;
    double M = 68, x1 = -10 + rand() % 20, y1 = -10 + rand() % 20, al = (1 - EE) / M, x0, y0;
    do
    {
        n++;
        if (n % 50 == 0 or n == 1)
            data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
        x0 = x1;
        y0 = y1;
        x1 = x0 - al * dfx(x0, y0, fun);
        y1 = y0 - al * dfy(x0, y0, fun);

    } while (((abs(func(x0, y0, fun) - func(x1, y1, fun)) > EE) and (sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2))) > EE) or (sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2))) > EE);
    data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
    cout << "__________________________________________________________" << endl;
    cout << "Answer: x = " << x1 << ", y = " << y1 << ", f(x,y) = " << func(x1, y1, fun) << endl
         << endl;
}

void grad_fast(int fun)
{
    srand(time(0) * 4);
    output_gradient();
    int n = 0;
    double x1 = -10 + rand() % 20, y1 = -10 + rand() % 20, al = 1, x0, y0;
    do
    {
        n++;
        if (n % 5 == 0 or n == 1)
            data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
        x0 = x1;
        y0 = y1;
        al = dihit(x0, y0, 0, 100, 1, 3);
        x1 = x0 - al * dfx(x0, y0, fun);
        y1 = y0 - al * dfy(x0, y0, fun);

    } while (((abs(func(x0, y0, fun) - func(x1, y1, fun)) > EE) and (sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2))) > EE) and (sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2))) > EE);
    data_gradient(n, x1, y1, func(x1, y1, fun), sqrt(pow(dfx(x1, y1, fun), 2) + pow(dfy(x1, y1, fun), 2)));
    cout << "__________________________________________________________" << endl;
    cout << "Answer: x = " << x1 << ", y = " << y1 << ", f(x,y) = " << func(x1, y1, fun) << endl
         << endl;
}

int main()
{
    cout << "Сoordinate descent method\nFor the 1st function\n";
    coordinate_descent(1);
    cout << "For the 2nd function\n";
    coordinate_descent(2);
    cout << "Step split gradient method\n";
    grad_split(1);
    cout << "Constant step gradient method\n";
    grad_const(1);
    cout << "Steepest gradient descent method\n";
    grad_fast(1);
    system("pause");
    return 0;
}