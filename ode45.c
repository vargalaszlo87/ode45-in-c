/*!
 * @project Ordinary Differential Equations in pure C v.0.1
 * @file ode45.c
 * @brief This is a Runge-Kutta solver of differential equation with Dormand-Prince coeffitient.
 *
 * Version information:
 *
 * In this version 0.1 you can only use one differential equation at a time.
 * Not working with input vector.
 *
 * @author Varga Laszlo
 *
 * @website https://github.com/vargalaszlo87/ode45-in-c
 * @website http://vargalaszlo.com
 * @website http://ha1cx.hu
 *
 * @date 2025-03-25
 *
 * @license
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define STEP_DIVIDER 1e4
#define STEP_MIN 1e-6
#define STEP_MAX_DIVIDER 5.0

/*!
 * Dormand-Prince coeffitient
 * The Butcher tableau is:
 *
 *	0
 *	1/5 	1/5
 *	3/10 	3/40 	    9/40
 *	4/5 	44/45 	    −56/15 	    32/9
 *	8/9 	19372/6561 	−25360/2187 64448/6561 	−212/729
 *	1 	    9017/3168 	−355/33 	46732/5247 	49/176 	    −5103/18656
 *	1 	    35/384 	    0 	        500/1113 	125/192 	−2187/6784 	    11/84
 *	----------------------------------------------------------------------------------------
 *          35/384 	    0 	        500/1113 	125/192 	−2187/6784 	    11/84 	    0
 *          5179/57600 	0 	        7571/16695 	393/640 	−92097/339200 	187/2100 	1/40
 *
 */

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

/*!
 *   Runge-Kutta method with adaptiv step size
 *
 *   Arguments:
 *
 *   double (*f)     pointer of differential equation
 *   double          start time
 *   double          value of y at the start time
 *   double          end time
 *   double          tolerance
 *
 */

void ode45(double (*f)(double, double), double t0, double y0, double tf, double tol) {

    // pointer is not null
    assert(f != NULL);

    // t_start is smaller than t_end
    assert(t0 < tf);

    // the tolerance is positiv
    assert(tol > 0);

    // local variables
    double
        t = t0,
        y = y0,
        h = (tf - t0) / STEP_DIVIDER,
        h_min = STEP_MIN,
        h_max = (tf - t0) / STEP_MAX_DIVIDER;

    // main cycle
    while (t < tf) {
        if (h < h_min) {
            printf("Small stepsize.\n");
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

        // 4th and 5th grade estimations
        double y4 = y, y5 = y;
        for (int i = 0; i < 7; i++) {
            y4 += b[i] * k[i];
            y5 += bs[i] * k[i];
        }

        // error estimation
        double error = fabs(y5 - y4);

        // step estimation
        double s = pow(tol / (error + 1e-10), 0.2);
        double h_new = h * fmin(5.0, fmax(0.1, s));

        if (error > tol) {
            h = fmax(h_new, h_min);
            continue;  // again --> smaller step
        }

        // next step
        t += h;
        y = y5;
        printf("t = %.5f, y = %.8f, error = %.2e, h = %.5e\n", t, y, error, h);

        // update step
        h = fmin(fmax(h_new, h_min), h_max);
        if (t + h > tf) h = tf - t;
    }
}
