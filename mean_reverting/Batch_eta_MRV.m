%% Batch_eta_MRV.m
%  Re-written on Mar. 4th, 2022
%  This code is for calculating multiple eta's at a time and plotting
%  results for different reverting speeds.
%
%  See OU_demo_batch_eta.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
load ini_val_eta_batch.mat
r = 0.04;
sgm = 0.15;

eta_set = [0.001:0.001:0.01,0.02,0.025:0.001:0.029,0.03:0.002:0.05];
tht_st = 10;
ph = 0.1;
bt = 1;
c = 0.1;

n = 10;
I = 2;
E = -1;

%% Batching
sol = zeros(4,length(eta_set));
for i = 1:length(eta_set)
    eta = eta_set(i);
    
    [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
        sgm, eta, tht_st, ph, bt, c, n, I, E);
    
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    F_tmp = @(x) F(x,prm);
    gradF_tmp = @(x) gradF(x,prm);
    
    %% Solve the non-linear equation system
    x0 = x0_set(:,i); % initial values
    
%     x_new = Homotopy(F_tmp,gradF_tmp,x0); x_new = MultiNewton(F_tmp,gradF_tmp,x_new);
    x_new = fsolve(F_tmp,x0,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',4e+3));
    
    sol(:,i) = x_new;
end

%% Plotting
% upper and lower bound
figure(1)
hold on
plot(eta_set,sol(3,:),'lineWidth',2)
plot(eta_set,sol(4,:),'lineWidth',2)
grid minor
xlabel('\eta')
legend(['\omega_h';'\omega_l'],'Location','best')
% print('./figures/sensitivity_analysis_eta','-dpng','-r600')

% firm value
plot_ind = 5;
[omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
    sgm, eta_set(plot_ind), tht_st, ph, bt, c, n, I, E);
prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
MRV_PlotFirmValue(-sol(1,plot_ind), -sol(2,plot_ind), ...
    sol(3,plot_ind), sol(4,plot_ind), prm);