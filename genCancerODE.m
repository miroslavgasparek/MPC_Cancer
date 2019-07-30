%% 29 July 2019 Miroslav Gasparek

function dxdt = genCancerODE(x,u,sys)
    % This function represents nonlinear model of the cancer proliferation
    % under the combined chemotherapy and immune therapy
    
    % Define the ODE
    dxdt =zeros(2,1);
    
    % Define the equations of the continuous state evolution
    dxdt(1) = - sys.uC * x(1) * log(x(1)/sys.x_inf) - sys.gamma * x(1)*x(2) - sys.k_x * x(1) * u(1);
    dxdt(2) = sys.uI * (x(1) - sys.beta*x(1)^2)*x(2) - sys.delta*x(2) + sys.alpha + sys.k_y * x(2) * u(2);

end