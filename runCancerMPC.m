%% 29 July 2019 Miroslav Gasparek
% Example of the linearized constrained Model Predictive Control (MPC) of 
% the nonlinear cancer treatment model with combined chemotherapy and
% immunotherapy.
% 
% Based on the following paper:
% Sharifi, N., Ozgoli, S., & Ramezani, A. (2017). 
% Multiple model predictive control for optimal drug administration 
% of mixed immunotherapy and chemotherapy of tumours. 
% Computer Methods and Programs in Biomedicine, 144, 13–19. 
% https://doi.org/10.1016/J.CMPB.2017.03.012
% 
% The equations of the cancer dynamics have the following form:
% 
%   dx/dt = - uC*x*ln(x / x_inf) - gammma*x*y - k_x*x*u
%   dy/dt = uI*(x-beta*x^2) - delta*x + alpha + k_y*y*v
%
% where:
% x: tumor volume (10^6 cells)
% y: immune-competent cells density (non-dim.)
% alpha: natural rate of influx of immune competent cells (1/day)
% beta: inverse threshold for the tumor suppresion (non-dim.)
% gamma: interaction rate between the immune comp. cells and tumor (10^7 cells/day)
% delta: death rate of immune cells (1/day)
% uC: tumor growth parameter (10^7 cells/day)
% uI: tumor stimulated proliferation rate (10^7 cells/day)
% x_inf: tumor carrying capacity
% k_x: killing parameter of chemotherapy wrt. tumor cells (10^7 cells/day)
% k_y: rate of immune cells proliferation when immunotherapy is used (non-dim.)


% Model Predictive Control scheme can be formulated as an optimization 
% problem with horizon length np, using the standard iterative LQR scheme:
%
%   Minimize J = 0.5*(X-rr)'*Q*(X-rr) + 0.5u'*R*u = 0.5*u'*H*u + x0*G'*u 
%   
%       subject to
%            cl <= Dcon * x <= ch
%            ul <=  u <= uh
%            x(k+1) = f(x(k),u(k))
%            X = [x(1); x(2);...;x(np)];
clear; clc; close all;

%% Parameters definition
sys.alpha = 0.1181; % 1/day
sys.beta = 0.00264; % 
sys.gamma = 1; % 10^7 cells/day
sys.delta = 0.37451; % 1/day
sys.uC = 0.5599; % 10^7 cells/day
sys.uI = 0.00484; % 10^7 cells/day
sys.x_inf = 780; % 10^6 cells
sys.k_x = 1; % 10^7 cells/day
sys.k_y = 0.1; %  1/day

%% Output matrices
C = eye(2); % We assume that both states are directly observable
D = zeros(2,2); % We assume that there is no direct term

%% Initial Conditions
% Initial tumour mass
x0 = 800; % 10^6
% Initial number of immune effector cells 
y0 = 0.05; % Non-dimensional

% Initial state vector
x = [x0; y0];

%% Constraints
% The constraints below are inserted in the form of:
% 
% cl <= Dcon * x <= ch
% ul <=  u <= uh
%
% Constraints on the states are such that states representing 
% the cell counts cannot be negative
% Also, the tumour mass should not exceed some maximum value
% and immune effector cells 
y_min = 0.001; % Non-dimensional
y_max = 2; % Non-dimensional

x_min = 0.1; % Nonzero to avoid simulation difficulties
x_max = 900; % 10^6 cells

cl = [0.01; y_min];
ch = [x_max; y_max];

% Constraints are independent, hence Dcon is identity matrix
Dcon = eye(2);

% Constraints on the inputs' min and max values
ui=[0.01; 0.01];         % initial zero control inputs, set to very small
                           % value to avoid the constraints
umin = [0; 0]; % constraint, minimum input values
umax = [1; 1]; % constraint, maximum input values

%% MPC simulation parameters
np = 10;        % horizon length 
nx = 2;         % number of states 
nu = 2;         % number of inputs
no = size(C,1); % number of outputs
Ts = 0.01;     % step size
Tfinal = 30;     % final time

%% Declare penalty matrices
q_x = 5;
q_y = 5;
Q = [q_x, 0;
     0  , q_y]; % relative importance of states

r_x = 25;
r_y = 1;

R = [r_x, 0;
     0,   r_y];  % penalizing weights of control inputs
 
P = 1*eye(2,2); % Low terminal cost

%% Reference generation
% Generating simple step-like reference
% Additional np points provided at the end due to the prediction horizon 
% Total desired tumour volume is set to decrease
x_Target1 = 300; % mm^3
x_Target2 = 50; % mm^3

% Total desired effector cells density
EC_Target = 0.9; % Arbitrary units

ref = [ [x_Target1 * ones(1,(Tfinal/Ts + np)/2), x_Target2 * ones(1, (Tfinal/Ts + np)/2)];...
       EC_Target*ones(1,Tfinal/Ts + np)];

%% Model generation
% The model is an anonymous function that represents the equation which describes 
% the dynamical system in the form of dx/dt =f(x,u)
model = @(x,u) genCancerODE(x,u,sys);

% Initialization of vectors for reference and output results 
rr = zeros(np*nx,1);
y  = zeros(no,Tfinal/Ts);
uh = zeros(nu,Tfinal/Ts);

for t=1:Tfinal/Ts
    
    % Get the reference trajectory for the whole horizon
    rr =  ref_for_hor(rr,ref,t,np,nx);
    
    % Evaluate the system output
    y(:,t) = C*x+D*ui;
    
    % Simulate one step forward with Runge-Kutta 4 order integrator
    [x, dx] = RK4(x,ui,Ts,model);
    
    % Calculate the linearized and discretized state & input matrices of 
    % the system at the current state
    [A, B] = linearizeCancerODE(x, ui, sys, Ts);
    
    % Compute stage constraint matrices and vector over the prediction
    % horizon
    [Dt,Et,bt]=genStageConstraints(A,B,Dcon,cl,ch,umin,umax);
    
    % Compute trajectory constraints matrices and vector over the
    % prediction horizon
    [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,np);
    
    % Compute prediction matrices over the prediction horizon
    [Gamma,Phi] = genPrediction(A,B,np);

    % Compute QP constraint matrices over the prediction horizon
    [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,np);
    
    % Compute QP cost matrices
    [H,G] = genCostMatrices(Gamma,Phi,Q,R,P,np);
    
    % Prepare cost and constraint matrices for mpcqpsolver
    % Calculate the inverse of the lower triangular H. See doc  for 
    % mpcqpsolver.
    H = chol(H,'lower');
    H = (H'\eye(size(H)))';
    % Set the constraints to be inactive at this point for the 
    % first pass, then use the constraints generated by the solver
    if t == 1
        iA = false(size(bb));
    else
        iA = iA_new;
    end
    
    % Generate the optimal input step
    [u,status,iA_new] = genMPController(H,G,F,bb,J,L,x,rr(1:nx),2,iA);
    
    % Use zero input when the solution is infeasible or num. error occurs
    if status == -1 || status == -2
        ui = zeros(nu,1);
    end
    
    status_vec(t) = status; % Store the status of the solver at each step
    ui = u(1:nu);           % providing first solution as input to our system
    uh(:,t) = ui;           % storing input
    
end

%% Plotting of the results
% Get the time vector for plotting
tt = Ts:Ts:Tfinal;

figure(1);
sgtitle('Tumor treatment optimisation using MPC',...
        'interpreter','latex','fontsize',15)
% Plot the pendulum's position
subplot(4,1,1)
hold on
plot(tt, ref(1,1:(Tfinal/Ts)),'r--','LineWidth',2)
plot(tt, cl(1)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, ch(1)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,y(1,:),'LineWidth',2)
xlabel('Time [Days]','interpreter','latex','fontsize',10)
ylabel('Tumor mass [$10^{6} cells$]','interpreter','latex','fontsize',10)
title('Tumor mass evolution','interpreter','latex','fontsize',12)
legend('ref. trajectory','constraints','Tumor mass','interpreter','latex',...
       'location','east','fontsize',6)
ax = gca;
ax.YLim = [-0.1, 1000];

% Plot the pendulum's angular velocity
subplot(4,1,2)
hold on
plot(tt, EC_Target*ref(2,1:(Tfinal/Ts)),'r--','LineWidth',2)
plot(tt, cl(2)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, ch(2)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,y(2,:),'LineWidth',2)
xlabel('Time [s]','interpreter','latex','fontsize',10)
ylabel('Norm. EC mass [AU]','interpreter','latex','fontsize',10)
title('Effector immune cells evolution','interpreter','latex','fontsize',12)
legend('ref. trajectory','constraints','EC','interpreter','latex',...
        'location','northeast','fontsize',6)
ax = gca;
ax.YLim = [-0.1, 2];

% Plot the normalized input force
subplot(4,1,3)
hold on
plot(tt, umin(1)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, umax(1)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,uh(1,:),'LineWidth',2)
xlabel('Time [s]','interpreter','latex','fontsize',10)
ylabel('Cytotoxic agent [AU]','interpreter','latex','fontsize',10)
title('Optimized Chemotherapeutic treatment','interpreter','latex','fontsize',12)
legend('constraints','CA [AU]','interpreter','latex',...
       'location','southeast','fontsize',6)
ax = gca;
ax.YLim = [-0.1, max(uh(1,:))+1];

subplot(4,1,4)
hold on
plot(tt, umin(2)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, umax(2)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,uh(2,:),'LineWidth',2)
xlabel('Time [s]','interpreter','latex','fontsize',10)
ylabel('IS agent [AU]','interpreter','latex','fontsize',10)
title('Optimized Immunotherapy treatment','interpreter','latex','fontsize',12)
legend('constraints','IA [AU]','interpreter','latex',...
       'location','northeast','fontsize',6)
ax = gca;
ax.YLim = [-0.1, max(uh(2,:))+1];

fig = gcf;
fig.Position = [545.0000   61.0000  891.2000  721.6000];


