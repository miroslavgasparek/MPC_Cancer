function dxdt = genTumourODE(x,u,sys)
    % This function represents nonlinear model of simple mathematical
    % pendulum in medium with linear viscocity 
    
%     % Define the parameters
%     sys.alpha = 0.1181; % 1/day
%     sys.beta = 0.00264; %
%     sys.gamma = 1; % 10^7 cells/day
%     sys.delta = 0.37451; % 1/day
%     sys.uC = 0.5599; % 10^7 cells/day
%     sys.uI = 0.00484; % 10^7 cells/day
%     sys.x_inf = 780; % 10^6 cells
%     sys.k_x = 1; % 10^7 cells/day
%     sys.k_y = 0.3; %  1/day
    
    % Define the ODE
    dxdt =zeros(2,1);
    
    % Define the equations of the continuous state evolution
    dxdt(1) = - sys.uC * x(1) * log(x(1)/sys.x_inf) - sys.gamma * x(1)*x(2) - sys.k_x * x(1) * u(1);
    dxdt(2) = sys.uI * (x(1) - sys.beta*x(1)^2)*x(2) - sys.delta*x(2) + sys.alpha + sys.k_y * x(2) * u(2);
    
%     dxdt(3) = - km * x(3) * log(theta/( x(1) + x(2) + x(3) )) + ks * x(2);
%     dxdt(4) = - k01 * x(4) + u;
%     dxdt(5) = - k12 * x(5) + k01 * x(4);
%     dxdt(6) = - kr2 * x(6) - k23 * x(6) + k12 * x(5);
%     dxdt(7) = - kr3 * x(7) + k23 * x(6);
%     dxdt(8) = alpha - beta * x(8) - kc * x(8) * (x(6)/V + sb*x(7)/V);
%     dxdt(9) = x(1) + x(2) + x(3);
end