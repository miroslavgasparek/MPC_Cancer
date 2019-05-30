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
%     f_state3 = - km * x(3) * log(theta/( x(1) + x(2) + x(3) )) + ks * x(2);
%     f_state4 = - k01 * x(4) + u;
%     f_state5 = - k12 * x(5) + k01 * x(4);
%     f_state6 = - kr2 * x(6) - k23 * x(6) + k12 * x(5);
%     f_state7 = - kr3 * x(7) + k23 * x(6);
%     f_state8 = alpha - beta * x(8) - kc * x(8) * (x(6)/V + sb*x(7)/V);
%     f_state9 = x(1) + x(2) + x(3);
    
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
