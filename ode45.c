#include <stdio.h>
#include <math.h>

// diff egyenlet: dy/dt = -2*y^2
double f(double t, double y) {
    return -2.0 * pow(y,2);
}

// Runge-Kutta-Dormand-Prince (4,5) egyutthatok
const double a[8][7] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0},
    {0.0, 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0},
    {0.0, 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0},
	{0.0, 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}
};

const double b[7] = {
	0.0, 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
};

const double bs[8] = {
	0.0, 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0
};

const double c[8] = {
	0.0, 0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0
};

// adaptiv Runge-Kutta 4(5) solver
void rk45(double (*f)(double, double), double t0, double y0, double tf, double tol) {
    double t = t0;
    double y = y0;
    double h = (tf - t0) / 100.0;  // Kezdeti lépésköz
    double h_min = 1e-6;
    double h_max = 0.1;

    while (t < tf) {
        
        // Runge-Kutta egyutthatok kiszamitasa
        double k[7];
        k[1] = h * f(t, y);
        for (int i = 2; i <= 6; i++) {
            double y_temp = y;
            for (int j = 1; j < i; j++) {
                y_temp += a[i][j] * k[j];
            }
            k[i] = h * f(t + c[i] * h, y_temp);
        }

        // negyed es otodrendu becsles
        double y4 = y;
        double y5 = y;
        for (int i = 1; i <= 6; i++) {
            y4 += b[i] * k[i];
            y5 += bs[i] * k[i];
        }
        y5 += bs[7] * h * f(t + h, y4);

        // hibabecsles
        double error = fabs(y5 - y4);
        
        // uj lepeskoz
        double s = pow(tol / (error + 1e-10), 0.2);
        h = h * fmin(5.0, fmax(0.1, s));  // lepcsokorrekcio

        // Ha tul nagy a hiba, ujraszamolas kisebb lepessel
        if (error > tol) continue;

        // frissites
        t += h;
        y = y5;

        // kiiras
        printf("t = %.5f, y = %.5f, h = %.5f\n", t, y, h);

        // Ne lepjunk tul a vegso idon
        if (t + h > tf) h = tf - t;
    }
}

// Főprogram
int main() {
    double t0 = 0.0;
    double y0 = 1.0;
    double tf = 5.0;
    double tol = 1e-5;

    printf("Adaptiv Runge-Kutta 4(5) megoldas:\n");
    rk45(f, t0, y0, tf, tol);

    return 0;
}
