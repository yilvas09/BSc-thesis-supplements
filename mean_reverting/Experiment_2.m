%% Experiment 2: Replicating Dilution Effect
%  Written on Mar. 10th, 2022
%  This code performs the numerical experiment 2 that replicates dilution
%  effect by fixing n and letting the mean market size theta^* vary.

%  See OU_demo_batch_tht_st.m & plot_beta_batch_tht_st.m for original
%  version
clear, close all
addpath('./auxilaries')
addpath('./figures')
%% Parameters
r = 0.04;
sgm = 0.15;

eta = 0.01;
tht_st_set = [55/3,10,110/21];
ph = 0.1;
bt = 1;
c = 0.1;

n = 10;
I = 2;
E = -1;

%% Batching
sol = zeros(4,length(tht_st_set));
prm_set = zeros(length(tht_st_set),11);
for i = 1:length(tht_st_set)
    tht_st = tht_st_set(i);
    [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, ...
        sgm, eta, tht_st, ph, bt, c, n, I, E);
    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n];
    F_tmp = @(x) F(x,prm);
    gradF_tmp = @(x) gradF(x,prm);
    
    %% Solve the non-linear equation system
    x0 = [18;-0.1;0.5;0.3] + 0.05*rand(4,1) - 0.025;
    if tht_st>18
        x0 = [14;-0.1;0.5;0.3] + 0.05*rand(4,1) - 0.025;
    end
    
%     x_new = Homotopy(F_tmp,gradF_tmp,x0); x_new = MultiNewton(F_tmp,gradF_tmp,x_new);
    x_new = fsolve(F_tmp,x0,optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',4e+3));
    sol(:,i) = x_new;
    prm_set(i,:) = prm;
end

%% Plotting
figure(1)
hold on
for i = 1:length(tht_st_set)
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
legend_name = cell(1,length(tht_st_set));
for i = 1:length(tht_st_set)
    legend_name{i} = ['\theta^*=',num2str(tht_st_set(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('V(\omega; n)')
print('./figures/val_func_tht_18_10_5','-dpng','-r600')

close all
figure(1)
hold on
for i = 1:length(tht_st_set)
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
legend_name = cell(1,length(tht_st_set));
for i = 1:length(tht_st_set)
    legend_name{i} = ['\theta^*=',num2str(tht_st_set(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('\beta(\omega; n)')
print('./figures/beta_tht_18_10_5','-dpng','-r600')
