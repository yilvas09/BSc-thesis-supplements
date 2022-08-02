%% Plot_Value_Functions_MRV.m
%  Re-written on Mar. 10th, 2022
%  This code is for plotting values of active firms under different
%  numbers of active firms in the industry as well as the perfect
%  competition case. It produces Figure 4 in the paper.

%  See plot_val_func.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Load data & parameters
load('solution_3to10000_Inf.mat')
n_set = sol(1,:);
plot_n = [3,6,10,Inf];

%% Plot
figure(1)
hold on
for k = 1:length(plot_n)
    kk = find(plot_n(k)==n_set);
    n = plot_n(k);
    B_1ht = -sol(2,kk);
    B_2 = -sol(3,kk);
    omg_h = sol(4,kk);
    omg_l = sol(5,kk);
    
    [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
        sgm, eta, tht_st, ph, bt, c, n, I, E);
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    
    domg = 1e-3;
    omg = 0.5*omg_l:domg:1.1*omg_h;
    boundH_ind = find(abs(omg-omg_h)==min(abs(omg-omg_h)));
    boundL_ind = find(abs(omg-omg_l)==min(abs(omg-omg_l)));
    
    [V,~] = MRV_ComputeFirmValue(omg,B_1ht,B_2,prm);
    V(1 : boundL_ind - 1) = NaN;
    V(boundH_ind+1 : end) = NaN;
    
    if k==1
        plot(omg,V,'lineWidth',2)
    elseif k==2
        plot(omg,V,'-.','lineWidth',2)
    elseif k==3
        plot(omg,V,':','lineWidth',2)
    else
        plot(omg,V,'--','lineWidth',2)
    end
end
hold off
grid minor
legend_name = cell(1,length(plot_n));
for i = 1:length(plot_n)
    legend_name{i} = ['n=',num2str(plot_n(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('V(\omega; n)')
% print('./figures/val_func_MRV_3_6_10_Inf','-dpng','-r600')
