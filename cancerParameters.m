%% 30 July 2019 Miroslav Gasparek
% The structure with parameters of the cancer treatment model, where:
%
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
% 
% Parameters adopted from the following paper, except of k_y, which was 
% not mentioned in the paper:
%
% Sharifi, N., Ozgoli, S., & Ramezani, A. (2017). 
% Multiple model predictive control for optimal drug administration 
% of mixed immunotherapy and chemotherapy of tumours. 
% Computer Methods and Programs in Biomedicine, 144, 13–19. 
% https://doi.org/10.1016/J.CMPB.2017.03.012
function system = cancerParameters()
    system.alpha = 0.1181; % 1/day
    system.beta = 0.00264; %
    system.gamma = 1; % 10^7 cells/day
    system.delta = 0.37451; % 1/day
    system.uC = 0.5599; % 10^7 cells/day
    system.uI = 0.00484; % 10^7 cells/day
    system.x_inf = 780; % 10^6 cells
    system.k_x = 1; % 10^7 cells/day
    system.k_y = 0.1; %  1/day
end