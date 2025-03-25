# ode45-in-c
Adaptiv Runge-Kutte method in C. 

## benchmark with ode45 from OCTAVE.

I tested this differential equation with my ode45 and ode45 from Octave:

```
  dydt = (-2 * y^2 * sin(t)) * e^(-2*y);
```

This is the Octave source code:

```MATLAB
  close;
  clear;
  clearAllMemoizedCaches;
  clearvars;
  
  function dydt = f(t,y)
  dydt = (-2 * y^2 * sin(t)) * e^(-2*y);
  end
  
  T0 = 0.0;
  Y0 = 1.0;
  TF = 1.0;
  
  TSPAN = [T0,TF];
  
  [tout, yout] = ode45(@f, TSPAN, Y0);
  
  plot(tout, yout);
```

This is an example of differential equation with my ode45:

```C
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

```

The output of the **C** code:

```
t = 0.00010, y = 1.00000000, error = 1.11e-016, h = 1.00000e-004
t = 0.00060, y = 0.99999995, error = 0.00e+000, h = 5.00000e-004
t = 0.00310, y = 0.99999870, error = 0.00e+000, h = 2.50000e-003
t = 0.01560, y = 0.99996707, error = 1.11e-016, h = 1.25000e-002
t = 0.07810, y = 0.99917493, error = 1.80e-013, h = 6.25000e-002
t = 0.27810, y = 0.98960088, error = 2.35e-010, h = 2.00000e-001
t = 0.47810, y = 0.96965936, error = 3.80e-010, h = 2.00000e-001
t = 0.67810, y = 0.94019208, error = 2.87e-010, h = 2.00000e-001
t = 0.87810, y = 0.90250940, error = 1.99e-010, h = 2.00000e-001
t = 1.00000, y = 0.87625068, error = 8.74e-011, h = 1.21900e-001
```

The output of the **Octave** example:

```
tout = 0.0062072, yout = 0.99999
tout = 0.015658, yout = 0.99997
tout = 0.27075, yout = 0.99014
tout = 0.47038, yout = 0.97061
tout = 0.67038, yout = 0.94149
tout = 0.87038, yout = 0.9041
```

### The results of benchmark:

** "t" and "y" come from my C code, "tout" and "yout" come from reference octave code. **

```
t = 0.00060, y = 0.99999995
tout = 0.0062072, yout = 0.99999

t = 0.01560, y = 0.99996707
tout = 0.015658, yout = 0.99997

t = 0.27810, y = 0.98960088
tout = 0.27075, yout = 0.99014

t = 0.47810, y = 0.96965936
tout = 0.47038, yout = 0.97061

t = 0.67810, y = 0.94019208
tout = 0.67038, yout = 0.94149

t = 0.87810, y = 0.90250940
tout = 0.87038, yout = 0.9041
```
