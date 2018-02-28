//Курсовая работа по Конечноэлементному моделированию
#include "stdafx.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#pragma warning(disable: 4996)

#include <math.h>
#include <locale.h>
#include <Windows.h>

//Интегральная функция для вычисления Im
double If(double x, double m)
{
    return pow(pow(sin(x), 3.0) / pow(x, 2.0), 1.0 / m) * (1.0 / x - 1.0 / tan(x));
}
//Вычисление интеграла Im методом трапеций
double findIm(double m, double err = 0.0001, double alph = M_PI_2)
{
    unsigned int n = 5, i;
    double sum = 0.0,
            suml = 0.0,
            x1,
            x2,
            h = alph / n,
            e;
    x1 = 1.0e-9;
    x2 = h;
    for(i = 0; i < n; i++)
    {
        suml = suml + (If(x1, m) + If(x2,m)) / 2.0 * h;
        x1 = h * (i + 1);
        x2 = h * (i + 2);
    }
    do{
        n = n * 2;
        sum = 0.0;
        h = alph / n;
        x1 = 1.0e-9;
        x2 = h;
        for(i = 0; i < n; i++)
        {
            sum = sum + (If(x1, m) + If(x2,m)) / 2.0 * h;
            x1 = h * (i + 1);
            x2 = h * (i + 2);
        }
        e = fabs(suml - sum);
        suml = sum;
    }while(e > err);
    return sum;
}
//Вычисление m
double findM(double p1, double p2, double t1, double t2)
{
    return log(p1/p2)/log(t2/t1);
}
//Вычисление K
double findK(double pf, double R0, double s0, double tf, double Im, double m)
{
    return (pf * R0 / (2.0 * s0)) * pow(tf / (2 * Im), m);
}

double Fm(double m, double pb[3], double tb[3])
{
	double f = 0.0;
	for (int i = 0; i < 3; i++)
		f += pow(pb[i] * pow(tb[i], m) - 1.0, 2.0);
	return f;
}

//Вычисление m 2м методом
void findM2(double p[3], double t[3], double R0, double s0, double eps = 1.0e-15, int itrs = 1000)
{
	double pb[3], tb[3], m;
	printf("Вычисление m лобовым методом введения опорной точки\n");
	printf("i m       K         Im    Fi         tf     итериций\n");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			pb[j] = p[j] / p[i];
			tb[j] = t[j] / t[i];
		}
		double a = 0.0, b = 2.0, alpha, beta;
		int its = 0;
		do
		{
			alpha = a + 0.382 * (b - a);
			beta =  a + 0.618 * (b - a);
			if (Fm(alpha, pb, tb) < Fm(beta, pb, tb))
			{
				b = beta;
				m = alpha;
			}
			else
			{
				a = alpha;
				m = beta;
			}
			its++;
		} while (its < itrs && (b - a) > eps);
		double Im = findIm(m);
		double K = p[i] * R0 / (2.0 * s0) * pow(t[i] / (2.0 * Im), m);
		double Fi = Fm(m, pb, tb);
		double tf = pow((2.0 * K * s0 / p[i] / R0), (1.0 / m)) * 2.0 * Im;
		printf("%d %7.5lf %9.3e %5.3lf %.2e %8.2lf %d\n", i + 1, m, K, Im, Fi, tf, its);
	}
}


int main()
{
    setlocale(LC_CTYPE, "RUS");
    FILE *out1, *out2, *out3;
    printf("\nРасчет m, K, n, C\n\n");
    double R0 = 0.0027, s0 = 0.049, p[3] = {5.0e5, 7.0e5, 1.0e6}, t[3] = {12232.0, 4341.0, 1447.0};
    double m[3], Im[3], K[3][3], n[3], C[3];
    printf("Вариант: 9\nR0 = %.3lf\ns0 = %.4lf\np1 = %.1e\tt1 = %.2lf\np2 = %.1e\tt2 = %.2lf\np3 = %.1e\tt3 = %.2lf\n\n",
           R0, s0, p[0], t[0], p[1], t[1], p[2], t[2]);
    printf("m         K1        K2        K среднее n         C\n");
    m[0] = findM(p[0], p[1], t[0], t[1]);
    Im[0] = findIm(m[0]);
    K[0][0] = findK(p[0], R0, s0, t[0], Im[0], m[0]);
    K[0][1] = findK(p[1], R0, s0, t[1], Im[0], m[0]);
    K[0][2] = (K[0][0] + K[0][1]) / 2.0;
    m[1] = findM(p[0], p[2], t[0], t[2]);
    Im[1] = findIm(m[1]);
    K[1][0] = findK(p[0], R0, s0, t[0], Im[1], m[1]);
    K[1][1] = findK(p[2], R0, s0, t[2], Im[1], m[1]);
    K[1][2] = (K[1][0] + K[1][1]) / 2.0;
    m[2] = findM(p[1], p[2], t[1], t[2]);
    Im[2] = findIm(m[2]);
    K[2][0] = findK(p[1], R0, s0, t[1], Im[2], m[2]);
    K[2][1] = findK(p[2], R0, s0, t[2], Im[2], m[2]);
    K[2][2] = (K[2][0] + K[2][1]) / 2.0;
	

    for(int i = 0; i < 3; i++)
    {
        n[i] = 1.0 / m[i];
        C[i] = 1.0 / pow(K[i][2], n[i]);
        printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", m[i], K[i][0], K[i][1], K[i][2], n[i], C[i]);
    }

	//Вычисление m 2м способом
	findM2(p, t, R0, s0);

    //Расчет  H(t) sa(t)
    out1 = fopen("Analytics_p5", "wt");
    out2 = fopen("Analytics_p7", "wt");
    out3 = fopen("Analytics_p10", "wt");
    printf("\nРасчет высоты купола H(t), толщины в полюсе sa(t), напряжения в полюсе Sigmaea(t)\n\n");
    int Nsh = 10;
    double shag = M_PI / 2.0 / double(Nsh), alpha = shag, time, H, sigea, ksiea, epsa, sa;
    printf("p = %.1e\nt, s\tH, m\tsa       Sigma\n", p[0]);
    for(int i = 0; i < Nsh; i++)
    {
        alpha = (i + 1) * shag;
        time = 2 * findIm(m[1],0.0001,alpha) / pow(p[0] * R0 / 2.0 / K[1][2] / s0, n[1]);
        H = tan(alpha / 2.0) * R0;
        sigea = (p[0] * R0 / 2.0 / s0) * pow(alpha, 2.0) / pow(sin(alpha), 3.0);
        ksiea = pow(sigea / K[1][2], n[1]);
        epsa = 2 * log(alpha / sin(alpha));
        sa = s0 * pow(sin(alpha) / alpha, 2.0);
        printf("%4.0lf\t%0.5lf\t%.2e %0.3e\n", time, H, sa, sigea);
        fprintf(out1, "%4.5lf\t%0.5lf\t%.5e\t%.5e\n", time, H, sa, sigea);
    }
    alpha = shag;
	printf("p = %.1e\nt, s\tH, m\tsa       Sigma\n", p[1]);
    for(int i = 0; i < Nsh; i++)
    {
        alpha = (i + 1) * shag;
        time = 2 * findIm(m[1],0.0001,alpha) / pow(p[1] * R0 / 2.0 / K[1][2] / s0, n[1]);
        H = tan(alpha / 2.0) * R0;
        sigea = (p[1] * R0 / 2.0 / s0) * pow(alpha, 2.0) / pow(sin(alpha), 3.0);
        ksiea = pow(sigea / K[1][2], n[1]);
        epsa = 2 * log(alpha / sin(alpha));
        sa = s0 * pow(sin(alpha) / alpha, 2.0);
		printf("%4.0lf\t%0.5lf\t%.2e %0.3e\n", time, H, sa, sigea);
        fprintf(out2, "%4.5lf\t%0.5lf\t%.5e\t%.5e\n", time, H, sigea, sa);
    }
    alpha = shag;
	printf("p = %.1e\nt, s\tH, m\tsa       Sigma\n", p[2]);
    for(int i = 0; i < Nsh; i++)
    {
        alpha = (i + 1) * shag;
        time = 2 * findIm(m[1],0.0001,alpha) / pow(p[2] * R0 / 2.0 / K[1][2] / s0, n[1]);
        H = tan(alpha / 2.0) * R0;
        sigea = (p[2] * R0 / 2.0 / s0) * pow(alpha, 2.0) / pow(sin(alpha), 3.0);
        ksiea = pow(sigea / K[1][2], n[1]);
        epsa = 2 * log(alpha / sin(alpha));
        sa = s0 * pow(sin(alpha) / alpha, 2.0);
		printf("%4.0lf\t%0.5lf\t%.2e %0.3e\n", time, H, sa, sigea);
        fprintf(out3, "%4.5lf\t%0.5lf\t%.5e\t%.5e\n", time, H, sigea, sa);
    }

    fclose (out1);
    fclose (out2);
    fclose (out3);
	system("Pause");

    return 0;
}

