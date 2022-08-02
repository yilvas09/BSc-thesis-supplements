%% Plot_Thresholds_MRV.m
%  Re-written on Mar. 10th, 2022
%  This code is for plotting entry and exit thresholds under different
%  numbers of active firms in the industry as well as their comparison with
%  the extreme case. It produces Figure 1 in the paper.

%  See plot_thres.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Load data & parameters
load('solution_3to10000_Inf.mat')
n_set = sol(1,:);
plot_n = [3:20,30:10:90];

plot_ind = zeros(size(plot_n));
for i = 1:length(plot_n)
    plot_ind(i) = find(plot_n(i)==n_set);
end

upr_bd = sol(4,plot_ind);
lwr_bd = sol(5,plot_ind);
omg_st = tht_st*(1-c) ./ plot_n;
gap1 = upr_bd - lwr_bd;

%% Plot
hold on
plot(plot_n,upr_bd,'-*','lineWidth',2)
plot(plot_n,lwr_bd,'-x','lineWidth',2)
plot(plot_n,gap1,'-k','lineWidth',1)
plot(plot_n,sol(4,end)*ones(size(plot_n)),'--b')
plot(plot_n,sol(5,end)*ones(size(plot_n)),'-.r')
plot(plot_n,(sol(4,end)-sol(5,end))*ones(size(plot_n)),':k')
grid minor
legend([{'\omega_h'}, {'\omega_l'}, {'\omega_h - \omega_l'},...
    {'\omega_h (n=\infty)'}, {'\omega_l (n=\infty)'},...
    {'\omega_h - \omega_l (n=\infty)'}],'Location','best')
xlabel('n')
% print('./figures/threshold_MRV','-r600','-dpng')