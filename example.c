#include <stdio.h>
#include <math.h>

#include "ode45.c"

double f(double t, double y) {
    return (-2 * pow(y,2) * sin(t)) * exp(-2 * y);
}

int main() {
    double
        t0 = 0.0,
        y0 = 1.0,
        tf = 1.0,
        tol = 1e-6;

    rk45(
        f, t0, y0, tf, tol
    );
    return 0;
}
