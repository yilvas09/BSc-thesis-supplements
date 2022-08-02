%% Batch_I_MRV.m
%  Re-written on Mar. 4th, 2022
%  This code is for calculating multiple I's at a time and plotting
%  results for different entry sunk costs. It produces Figure 10 in the
%  paper.
%
%  See OU_demo_batch_I.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
load ini_val_I_batch.mat
r = 0.04;
sgm = 0.15;

eta = 0.01;
tht_st = 10;
ph = 0.1;
bt = 1;
c = 0.1;

n = 10;
I_set = [1:0.02:1.08,1.1:0.1:1.4,1.5:0.25:3];
E = -1;

%% Batching
sol = zeros(4,length(I_set));
for i = 1:length(I_set)
    I = I_set(i);
    
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
plot(I_set,sol(3,:),'lineWidth',2)
plot(I_set,sol(4,:),'--','lineWidth',2)
grid minor
xlabel('I')
legend(['\omega_h';'\omega_l'],'Location','best')
print('./figures/sensitivity_analysis_I_MRV','-dpng','-r600')

% firm value
% plot_ind = 5;
% [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
%     sgm, eta, tht_st, ph, bt, c, n, I_set(plot_ind), E);
% prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I_set(plot_ind), E, n];
% MRV_PlotFirmValue(-sol(1,plot_ind), -sol(2,plot_ind), ...
%     sol(3,plot_ind), sol(4,plot_ind), prm);