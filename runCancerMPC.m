%% 29 May 2019 Miroslav Gasparek
% Example of the linearized constrained Model Predictive Control (MPC) of 
% the nonlinear pendulum with the linear damping.
% 
% The equations of the pendulum are of the form:
% 
%     dx(1)/dt = x(2);
%     dx(2)/dt = (-g/l)*sin(x(1)) -b*x(2)+ u;
%
% where:
% g is acceleration due to the gravity
% l is the length of the rod
% b is the damping coefficient
% u is the input, the normalized force

% Model Predictive Control scheme can be formulated as an optimization 
% problem with horizon length np:
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
sys.k_y = 0.7; %  1/day


%% Output matrices
C = eye(2); % We assume that both states are directly observable
D = zeros(2,2); % We assume that there is no direct term

%% Initial Conditions
x0 = 700; % 10^6
y0 = 0.20; % Non-dimensional

% Initial states
x = [x0; y0];

%% Constraints
% The constraints below are inserted in the form of:
% 
% cl <= Dcon * x <= ch
% ul <=  u <= uh
%
% Constraints on the states are such that states cannot be negative
% And they cannot go above the carrying capacity for the tumour and 
% they cannot go below the minimum amount of the cells 
y_min = 0.01; % Non-dimensional
y_max = 3;

x_max = 1000; % 10^6 cells

cl = [0.01; y_min];
ch = [x_max; y_max];

% Constraints are independent, hence Dcon is identity matrix
Dcon = eye(2);

% Constraints on the inputs' min and max values
ui=[0.001; 0.001];         % initial zero cotrol inputs
umin = [0; 0]; % ug, minimum input values
umax = [1; 1]; % ug, maximum input values

%% MPC simulation parameters
np = 10;        % horizon length 
nx = 2;         % number of states 
nu = 2;         % number of inputs
no = size(C,1); % number of outputs
Ts = 0.2;     % step size
Tfinal = 10;     % final time

%% Declare penalty matrices
q_x = 100;
q_y = 0.1;
Q = [q_x, 0;
     0  , q_y]; % relative importance of states

r_x = 10;
r_y = 2;

R = [r_x, 0;
     0,   r_y];  % penalizing weights of control inputs
 
P = 0.1*zeros(2,2); % Low terminal cost

%% Reference generation
% Generating simple step reference of the form [0.5 ; 0]
% Additional np points provided at the end due to the prediction horizon 
% Total desired tumour volume
x_Target = 20; % mm^3
ref = [x_Target * ones(1,Tfinal/Ts + np);...
       ones(1,Tfinal/Ts + np)];

%% Model generation
% The model is an anonymous function that represents the equation which describes 
% the dynamical system in the form of dx/dt =f(x,u), in this case, the
% nonlinear pendulum with linear damping
model = @(x,u) genTumourODE(x,u,sys);

% initialization of vectors for reference and output results 
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
    % Set the constraints to be inactive at this point
    iA = false(size(bb));
    
    % Generate the optimal input step
    [u,~,iA] = genMPController(H,G,F,bb,J,L,x,rr(1:nx),2,iA);
    ui = u(1:nu);     %providing first solution as input to our system
    uh(:,t) = ui;     %storing input
    
end

%% Plot the results
% Get the time vector for plotting
tt = Ts:Ts:Tfinal;

figure(1);
sgtitle('Tumour treatment optimisation using MPC',...
        'interpreter','latex','fontsize',15)
% Plot the pendulum's position
subplot(4,1,1)
hold on
plot(tt, ref(1,1:(Tfinal/Ts)),'r--','LineWidth',2)
plot(tt, cl(1)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, ch(1)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,y(1,:),'LineWidth',2)
xlabel('Time [Days]','interpreter','latex','fontsize',10)
ylabel('Tumour mass [$10^{6} cells$]','interpreter','latex','fontsize',10)
title('Tumour mass evolution','interpreter','latex','fontsize',12)
legend('ref. trajectory','constraints','Tumour mass','interpreter','latex',...
       'location','east','fontsize',6)
ax = gca;
ax.YLim = [-0.1, 1000];

% Plot the pendulum's angular velocity
subplot(4,1,2)
hold on
plot(tt, ref(2,1:(Tfinal/Ts)),'r--','LineWidth',2)
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
title('Optimalized Chemotherapeutic treatment','interpreter','latex','fontsize',12)
legend('constraints','CA [AU]','interpreter','latex',...
       'location','southeast','fontsize',6)
ax = gca;
ax.YLim = [-0.1, 1.2];

subplot(4,1,4)
hold on
plot(tt, umin(2)*ones(1,length(tt)),'k--','LineWidth',2)
plot(tt, umax(2)*ones(1,length(tt)),'k--','LineWidth',2,'HandleVisibility','off')
plot(tt,uh(2,:),'LineWidth',2)
xlabel('Time [s]','interpreter','latex','fontsize',10)
ylabel('IS agent [AU]','interpreter','latex','fontsize',10)
title('Optimalized Immunotherapy treatment','interpreter','latex','fontsize',12)
legend('constraints','IA [AU]','interpreter','latex',...
       'location','northeast','fontsize',6)
ax = gca;
ax.YLim = [-0.1, 2];

fig = gcf;
fig.Position = [545.0000   61.0000  891.2000  721.6000];


