#include <stdio.h>
#include <math.h>

// Dormand-Prince egyutthatok (RK45)
const double a[7][6] = {
    {},
    {0.2},
    {3.0/40.0, 9.0/40.0},
    {44.0/45.0, -56.0/15.0, 32.0/9.0},
    {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0},
    {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0},
    {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}
};

const double b[7]  = {
        35.0/384.0, 0, 500.0/1113.0, 125.0/192.0,-2187.0/6784.0, 11.0/84.0, 0
    };

const double bs[7] = {
        5179.0/57600.0, 0, 7571.0/16695.0, 393.0/640.0,-92097.0/339200.0, 187.0/2100.0, 1.0/40.0
};

const double c[]  = {
    0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0
};


// diff.egyenlet: dy/dt = -2*y^2
double f(double t, double y) {
    return -2.0 * y * y;
}

void rk45(double (*f)(double, double), double t0, double y0, double tf, double tol) {
    double t = t0;
    double y = y0;
    double h = (tf - t0) / 100.0;
    double h_min = 1e-6;
    double h_max = (tf - t0) / 5.0;

    while (t < tf) {
        if (h < h_min) {
            printf("Tul kicsi lepeskoz.\n");
            break;
        }

        double k[7];
        k[0] = h * f(t, y);
        for (int i = 1; i < 7; i++) {
            double y_temp = y;
            for (int j = 0; j < i; j++) {
                y_temp += a[i][j] * k[j];
            }
            k[i] = h * f(t + c[i] * h, y_temp);
        }

        // 4. es 5. rendu becselesek
        double y4 = y, y5 = y;
        for (int i = 0; i < 7; i++) {
            y4 += b[i] * k[i];
            y5 += bs[i] * k[i];
        }

        // hibabecsles
        double error = fabs(y5 - y4);

        // lepeskoz becsles
        double s = pow(tol / (error + 1e-10), 0.2);
        double h_new = h * fmin(5.0, fmax(0.1, s));

        if (error > tol) {
            h = fmax(h_new, h_min);
            continue;  // ujra --> kisebb lepessel
        }

        // sikeres lepes
        t += h;
        y = y5;
        printf("t = %.5f, y = %.8f, error = %.2e, h = %.5e\n", t, y, error, h);

        // frissitett lepeskoz
        h = fmin(fmax(h_new, h_min), h_max);
        if (t + h > tf) h = tf - t;
    }
}

int main() {
    double
        t0 = 0.0,
        y0 = 1.0,
        tf = 5.0,
        tol = 1e-6;
    
    rk45(
        f, t0, y0, tf, tol
    );
    return 0;
}
