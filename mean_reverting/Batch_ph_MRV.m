%% Batch_ph_MRV.m
%  Re-written on Mar. 4th, 2022
%  This code is for calculating multiple phi's at a time and plotting
%  results for different fixed costs. It produces Figure 12 in the paper.
%
%  See OU_demo_batch_ph.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
r = 0.04;
sgm = 0.15;

eta = 0.01;
tht_st = 10;
ph_set = 0:0.01:0.1;
bt = 1;
c = 0.1;

n = 10;
I = 2;
E = -1;

%% Batching
sol = zeros(4,length(ph_set));
for i = 1:length(ph_set)
    ph = ph_set(i);
    
    [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
        sgm, eta, tht_st, ph, bt, c, n, I, E);
    
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    F_tmp = @(x) F(x,prm);
    gradF_tmp = @(x) gradF(x,prm);
    
    %% Solve the non-linear equation system
    x0 = [18;-0.1;0.5;0.3] + 0.05*rand(4,1) - 0.025; % initial values
    
%     x_new = Homotopy(F_tmp,gradF_tmp,x0); x_new = MultiNewton(F_tmp,gradF_tmp,x_new);
    x_new = fsolve(F_tmp,x0,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',4e+3));
    
    sol(:,i) = x_new;
end

%% Plotting
% upper and lower bound
figure(1)
hold on
plot(ph_set,sol(3,:),'lineWidth',2)
plot(ph_set,sol(4,:),'--','lineWidth',2)
grid minor
xlabel('\phi')
legend(['\omega_h';'\omega_l'],'Location','best')
print('./figures/sensitivity_analysis_ph_MRV','-dpng','-r600')

% firm value
% plot_ind = 5;
% [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
%     sgm, eta, tht_st, ph_set(plot_ind), bt, c, n, I, E);
% prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
% MRV_PlotFirmValue(-sol(1,plot_ind), -sol(2,plot_ind), ...
%     sol(3,plot_ind), sol(4,plot_ind), prm);