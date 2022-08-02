%% Batch_n_MRV.m
%  Re-written on Mar. 10th, 2022
%  This code is for calculating multiple n's at a time and plotting
%  results for different numbers of active firms.
%
%  See OU_demo_batch_n.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
load ini_val_n_batch.mat
r = 0.04;
sgm = 0.15;

eta = 0.01;
tht_st = 10;
ph = 0.1;
bt = 1;
c = 0.1;

n_set = [3:20,30:10:100,500,1000,5000,10000];
I = 2;
E = -1;

%% Batching
sol = zeros(4,length(n_set));
for i = 1:length(n_set)
    n = n_set(i);
    
    [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
        sgm, eta, tht_st, ph, bt, c, n, I, E); %#ok<*ASGLU>
    
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    F_tmp = @(x) F(x,prm);
    gradF_tmp = @(x) gradF(x,prm);
    
    %% Solve the non-linear equation system
    %     x0 = [18;-0.1;0.5;0.3] + 0.05*rand(4,1) - 0.025;
    x0 = x0_set(:,i);
    
    %     x_new = Homotopy(F_tmp,gradF_tmp,x0); x_new = MultiNewton(F_tmp,gradF_tmp,x_new);
    x_new = fsolve(F_tmp,x0,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',4e+3));
    
    sol(:,i) = x_new;
end

% save(['solution_',num2str(min(n_set)),'to',num2str(max(n_set))],'sol', ...
%     'r', 'sgm', 'eta', 'tht_st', 'ph', 'bt', 'I', 'E', 'c')

% firm value
plot_ind = 5;
[omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
    sgm, eta, tht_st, ph, bt, c, n_set(plot_ind), I, E);
prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n_set(plot_ind)];
MRV_PlotFirmValue(-sol(1,plot_ind), -sol(2,plot_ind), ...
    sol(3,plot_ind), sol(4,plot_ind), prm);