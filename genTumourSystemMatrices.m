%% 21 May 2019 Miroslav Gasparek
% Symbolic expression of the linearized system for the model predictive
% control of the cancer treatment

function [A, B, C, D] = genTumourSystemMatrices(x_in, pars, Ts)

%% Define the ODEs for the system
x = sym('x', [9 1]);
u = sym('u', [1 1]);
% syms kg ks km k01 k12 k23 kr2 kr3 theta kd c alpha beta kc b V

% Define the equations of the continuous state evolution
f_state1 = - pars(1) * x(1) * log(pars(9)/(x(1) + x(2) + x(3) )) + 2 * pars(3) * x(3) * log(pars(9)/(x(1) + x(2) + x(3) )) - pars(10) * x(1) * (x(6)/pars(16) + pars(11)*x(7)/pars(16));
f_state2 = - pars(2) * x(2) + pars(1) * x(1) * log(pars(9)/( x(1) + x(2) + x(3) ));
f_state3 = - pars(3) * x(3) * log(pars(9)/( x(1) + x(2) + x(3) )) + pars(2) * x(2);
f_state4 = - pars(4) * x(4) + u;
f_state5 = - pars(5) * x(5) + pars(4) * x(4);
f_state6 = - pars(7) * x(6) - pars(6) * x(6) + pars(5) * x(5);
f_state7 = - pars(8) * x(7) + pars(6) * x(6);
f_state8 = pars(12) - pars(13) * x(8) - pars(14) * x(8) * (x(6)/pars(16) + pars(15)*x(7)/pars(16));
f_state9 = x(1) + x(2) + x(3);

f_state = [f_state1;
           f_state2;
           f_state3;
           f_state4;
           f_state5;
           f_state6;
           f_state7;
           f_state8;
           f_state9];

% Get the Jacobian matrices of the system
Ajac = jacobian(f_state, x);
Bjac = jacobian(f_state, u);

Afun = matlabFunction(Ajac, 'vars', {x});

% Define the state matrix 
A = Afun(x_in);

% Define the input matrix 
B = [0, 0, 0, 1, 0, 0, 0, 0, 0]';

% Define the observation matrix 
C = [zeros(5,9);
    0, 0, 0, 0, 0, 1/pars(16), 0, 0, 0;
      0, 0, 0, 0, 0, 0, 1/pars(16), 0, 0;
      0, 0, 0, 0, 0, 0, 0, 1, 0;
      0, 0, 0, 0, 0, 0, 0, 0, 1];
 
D = zeros(9,1);

% if Ts>0 then sample the model with a zero-order hold (piecewise constant) input, 
% otherwise return a continuous-time model
if Ts>0
    system = ss(A,B,C,D);
    system_d = c2d(system, Ts);
    
    A = system_d.A;
    D = system_d.D;
end

end
 