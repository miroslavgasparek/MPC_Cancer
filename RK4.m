function [x, dx] = RK4(x1,u,Ts,k)
% This function uses the 4th order Runge-Kutta to calculate the evolution
% of the states over the time period Ts

%Algorithm of RK4
k1 = k(x1,u);
k2 = k(x1 + (Ts/2)*k1,u);
k3 = k(x1 + (Ts/2)*k2,u);
k4 = k(x1+Ts*k3,u);

dx=k1;

x = x1 + Ts*(k1 + 2*k2 + 2*k3 + k4)/6;
end
