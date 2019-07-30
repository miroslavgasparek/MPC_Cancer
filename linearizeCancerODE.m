function [Alin, Blin] = linearizeCancerODE(x_in, u_in, sys, Ts)
    % This function is linearization of the continuous nonlinear model in the form
    % of dx/dt=f(x,u). Approximation takes form 
    % dx/dt = Alin * x + Blin * u
    %
    % Matrices Alin, Blin are Jacobians of f(x,u) about x and u
    % respectively and if the sampling period is provided, these matrices
    % are discretized

    % Define the symbolic variables
    x = sym('x', [2 1]);
    u = sym('u', [2 1]);
   
    % Define the equations of the continuous state evolution
    f_state1 = - sys.uC * x(1) * log(x(1)/sys.x_inf) - sys.gamma * x(1)*x(2) - sys.k_x * x(1) * u(1);
    f_state2 = sys.uI * (x(1) - sys.beta * x(1)^2)*x(2) - sys.delta*x(2) + sys.alpha + sys.k_y * x(2) * u(2);
    
    % Create the column vector of the functions
    f_state = [f_state1;
               f_state2;];
           
    % Get the Jacobian state matrix of the system
    Ajac = jacobian(f_state, x);
    Bjac = jacobian(f_state, u);
           
    Afun = matlabFunction(Ajac, 'vars', {x, u});
    % Define the input matrix
    Bfun = matlabFunction(Bjac, 'vars', {x, u});
    
    % Define the state matrix
    Alin = Afun(x_in, u_in);
    Blin = Bfun(x_in, u_in);
    
    % Convert the continuous system to discrete system if the sampling
    % period is provided
    c_sys = ss(Alin, Blin, [], []);

    if Ts > 0
        d_sys = c2d(c_sys, Ts);
        Alin = d_sys.A;
        Blin = d_sys.B;
    end
    
end
