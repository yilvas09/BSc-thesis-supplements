%% Experiment 1: Replicating Oligopoly Effect
%  Written on Mar. 10th, 2022
%  This code performs the numerical experiment 1 that replicates oligopoly
%  effect by fixing the mean margin and letting the boundary conditions
%  vary with n.

%  See test_OU_batch.m for original version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
% Input
r = 0.04;
sgm = 0.15;

eta = 0.01;
tht_st = 10;
ph = 0.1;
bt = 1;
c = 0.1;

n_set = [5,10,20];
I = 2;
E = -1;

%% Batching
sol = zeros(4,length(n_set));
prm_set = zeros(length(n_set),11);
for i = 1:length(n_set)
    % Intermediate parameters
    n = n_set(i);
    omg_st = tht_st * (1-c) / 10;
    
    gmm_1 = (2*eta + sgm^2 + sqrt(8*r*sgm^2 + (2*eta + sgm^2)^2)) / (2*sgm^2);
    gmm_2 = (2*eta + sgm^2 - sqrt(8*r*sgm^2 + (2*eta + sgm^2)^2)) / (2*sgm^2);
    
    b_1 = 2 - 2*gmm_1 + 2*eta/(sgm^2);
    b_2 = 2 - 2*gmm_2 + 2*eta/(sgm^2);
    
    xi = 2*eta*omg_st / sgm^2;
    R_1 = 1 / bt / (2*eta+r-sgm^2);
    R_2 = 2*eta*omg_st * (1/bt/(eta+r)/(eta-sgm^2) - 1/(2*eta+r-sgm^2)/(eta-sgm^2));
    R_3 = 2*eta*omg_st^2/bt * ((1/(eta-sgm^2) - 1/(2*eta-sgm^2))/(2*eta+r-sgm^2) -...
        1/(eta+r)/(eta-sgm^2) + 1/r/(2*eta-sgm^2)) - ph/r;
    
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    F_tmp = @(x) F(x,prm);
    gradF_tmp = @(x) gradF(x,prm);
    
    %% Solve the non-linear equation system
    x0 = [10;-0.3;1;0.5] + 0.05*rand(4,1) - 0.025;
    
    x_new = fsolve(F_tmp,x0,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',4e+3));
    
    sol(:,i) = x_new;
    prm_set(i,:) = prm;
end

%% Plotting
figure(1)
hold on
for i = 1:length(n_set)
    B_1ht = -sol(1,i);
    B_2 = -sol(2,i);
    omg_h = sol(3,i);
    omg_l = sol(4,i);
    
    domg = 1e-3;
    omg = 0.5*omg_l:domg:1.1*omg_h;
    boundH_ind = find(abs(omg-omg_h)==min(abs(omg-omg_h)));
    boundL_ind = find(abs(omg-omg_l)==min(abs(omg-omg_l)));
    [V,~] = MRV_ComputeFirmValue(omg,B_1ht,B_2,prm_set(i,:));
    V(1 : boundL_ind - 1) = NaN;
    V(boundH_ind+1 : end) = NaN;
    if i==1
        plot(omg,V,'lineWidth',2)
    elseif i==2
        plot(omg,V,'-.','lineWidth',2)
    elseif i==3
        plot(omg,V,':','lineWidth',2)
    else
        plot(omg,V,'--','lineWidth',2)
    end
end
hold off
grid minor
legend_name = cell(1,length(n_set));
for i = 1:length(n_set)
    legend_name{i} = ['n=',num2str(n_set(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('V(\omega; n)')
print('./figures/val_func_fxd_margin_5_10_20','-dpng','-r600')

close all
figure(1)
hold on
for i = 1:length(n_set)
    B_1ht = -sol(1,i);
    B_2 = -sol(2,i);
    omg_h = sol(3,i);
    omg_l = sol(4,i);
    
    domg = 1e-3;
    omg = 0.5*omg_l:domg:1.1*omg_h;
    boundH_ind = find(abs(omg-omg_h)==min(abs(omg-omg_h)));
    boundL_ind = find(abs(omg-omg_l)==min(abs(omg-omg_l)));
    [V,Vprime] = MRV_ComputeFirmValue(omg,B_1ht,B_2,prm_set(i,:));
    V(1 : boundL_ind - 1) = NaN;
    V(boundH_ind+1 : end) = NaN;
    
    beta = sgm .* omg .* Vprime ./ V;
    
    if i==1
        plot(omg,beta,'lineWidth',2)
    elseif i==2
        plot(omg,beta,'-.','lineWidth',2)
    elseif i==3
        plot(omg,beta,':','lineWidth',2)
    else
        plot(omg,beta,'--','lineWidth',2)
    end
end
hold off
grid minor
legend_name = cell(1,length(n_set));
for i = 1:length(n_set)
    legend_name{i} = ['n=',num2str(n_set(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('\beta(\omega; n)')
print('./figures/beta_fxd_margin_5_10_20','-dpng','-r600')
