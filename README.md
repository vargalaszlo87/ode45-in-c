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

This is the differential equation in my code:

```C
  double f(double t, double y) {
      return (-2 * pow(y,2) * sin(t)) * exp(-2 * y);
  }
```
