%% 29 July 2019 Miroslav Gasparek
% Script that analyzes the dynamical system describing the combined
% treatment of cance by chemotherapy and immunotherapy
% 
% The system equations are as follows:
% 
% The equations of the cancer dynamics have the following form:
% 
%   dx/dt = - uC*x*ln(x / x_inf) - gammma*x*y - k_x*x*u
%   dy/dt = uI*(x-beta*x^2) - delta*x + alpha + k_y*y*v
%
clear; clc; close all;

%% Parameters definition
sys = cancerParameters();


%% Find the fixed points of no-treatment model
% Assume that u = 0, v = 0
u_in = [0; 1];
x0 = [1000,5];
% Hide the text output of the fsolve
options = optimset('Display','off');


%% Drawing the phase plane for the different initial conditions
% Create the vector of different intial conditions for tumor mass 
% and immune effector cells numbers
tumor_vec = linspace(10,1000,8);
eff_cell_vec = linspace(0.1,5,8);

% Define the timespan, length of simulation
tspan = [0, 50]; % days

% Counter for the double for loop
k = 0;

% Run the loop with the different intial conditions
% and draw the phase plane
figure(1)
xlabel('Tumor mass [$10^{6}$ cells]','fontsize',15,'interpreter','latex')
ylabel('Immune effector cells','fontsize',15,'interpreter','latex')
title(['Phase plane of the cancer proliferation system, u = '...
       num2str(u_in(1)), ', v = ', num2str(u_in(2))],...
       'fontsize',15,'interpreter','latex')

hold on
for i=1:length(tumor_vec)
    for j=1:length(eff_cell_vec)
        % Increment the counter 
        k = k+1;
        % Pick the initial conditions
        x0 = [tumor_vec(i), eff_cell_vec(j)];
        % Run the ODE simulation
        sol = ode45(@(t,x)genCancerODE(x,u_in, sys), tspan, x0);
        % Get the fixed point coordinates for the given init. conditions
        [x_fix, ~] = fsolve(@(x) genCancerODE(x,u_in, sys), x0, options);
        % If the fixed points are non-negative, store them in a matrix 
        % of fixed points
        if all(x_fix >=0)
            x_fix_mat(:,k) = x_fix;
        end    
        % Plot the trajectory on the phase plane
        p = plot(sol.y(1,:), sol.y(2,:),'k-','LineWidth',1);
    end
end

% Get the coordinates and pick the unique fixed points
x_fix_round = round(x_fix_mat,2);
[tumor_fp, iA1] = unique(x_fix_round(1,:), 'stable');
[imm_cell_fp, iA2] = unique(x_fix_round(2,:), 'stable');
% Check if the fixed points are too close to each other

% Plot the fixed points on the phase plane
s = scatter(tumor_fp, imm_cell_fp,40,'filled', 'r');

% Plot the legend
legend([p(1), s(1)],'Trajectories','Fixed Points','interpreter','latex')
hold off
% End plotting

%% Stability analysis of the fixed points
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

% Define the state matrix
Afun = matlabFunction(Ajac, 'vars', {x, u});

% Compute the values of the state matrix for all the fixed points
% and compute the eigenvalues of these matrices
% Define an empty cell with the length of the fixed points
Acell = cell(length(tumor_fp),1);
Aeig = cell(length(tumor_fp),1);

for i=1:length(tumor_fp)
    % Get the state matrix
    Acell{i} = Afun([tumor_fp(i); imm_cell_fp(i)], u_in);
    % Get the eigenvalues of the state matrix
    Aeig{i} = eig(Acell{i});
    fprintf('Fixed point: \nTumour cells: %.2f, Immune cells: %.2f \n',...
            tumor_fp(i), imm_cell_fp(i));
    fprintf('Eigenvalues: %.2f + %.2fi; %.2f + %.2fi \n\n',...
            real(Aeig{i}(1)), imag(Aeig{i}(1)), real(Aeig{i}(2)), imag(Aeig{i}(2)))
end
