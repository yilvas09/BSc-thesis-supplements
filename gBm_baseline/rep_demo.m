%% Entry in Geometric Brownian Motion
%  Written on Dec.23rd, 2021
%  This code is for replicating the risk-loading functions w.r.t. omega
%  under different degrees of concentration (i.e. different numbers of
%  incumbents) in RFS2017.
clear, close all
%% Parameters
% input
bt = 1; % slope of the demand curve
r = 0.04; % risk-free interest
alp = -0.05; % drift rate of GBM
sgm = 0.1; % volatility rate of GBM
ph = 0.01; % fixed cost of a firm
I = 1; % entry cost
E = -1; % exit cost
p_rsk = 1; % price of risk
n_set = [3,6,10,100];

%% Batching and Plotting
figure(1)
hold on
for j = 1:length(n_set)
    n = n_set(j);
    % intermediate parameters
    b1 = 0.5 - alp/sgm^2 + sqrt((alp/sgm^2 - 0.5)^2 + 2*r/sgm^2); % entry
    b = 0.5 - alp/sgm^2 - sqrt((alp/sgm^2 - 0.5)^2 + 2*r/sgm^2); % exit
    dlt = r - sgm^2 - 2*alp;
    
%     omg_bar = sqrt(b/(b-2) * bt * dlt * I); % entry
%     omg_nbar = sqrt((b-2) / (b*(1-1/(n+2))^2 - 2*(1-1/(n+2))^b)) * omg_bar;
    
    omg_bar = sqrt(b/(b-2) * bt * dlt * (-E)); % exit
    omg_nbar = sqrt((b-2) / (b*(1+1/n)^2 - 2*(1+1/n)^b)) * omg_bar;
    
    % value functions and derivatives
    domg = 1e-3;
    omg_vec = 0.8*omg_nbar:domg:1.1*omg_nbar;
%     omg_vec = 0.3:domg:omg_nbar; % entry
%     omg_vec = omg_nbar:domg:0.4; % exit
    boundH_ind = find(abs(omg_vec-omg_nbar)==min(abs(omg_vec-omg_nbar)));
    
    V_vec = -2*omg_nbar^(2-b)/b/bt/dlt * omg_vec.^b + omg_vec.^2/bt/dlt - ph/r;
    Vprime_vec = 2/bt/dlt * (omg_vec - omg_nbar^(2-b) * omg_vec.^(b-1));
    Vpprime_vec = 2/bt/dlt * (1 - omg_nbar^(2-b)*(b-1) * omg_vec.^(b-2));
    rsk_loading = sgm * omg_vec .* Vprime_vec ./ V_vec;
    
%     V_vec(boundH_ind:end) = NaN; % entry
%     Vprime_vec(boundH_ind:end) = NaN;
%     Vpprime_vec(boundH_ind:end) = NaN;
%     rsk_loading(boundH_ind:end) = NaN;
    
    V_vec(1:boundH_ind) = NaN; % exit
    Vprime_vec(1:boundH_ind) = NaN;
    Vpprime_vec(1:boundH_ind) = NaN;
    rsk_loading(1:boundH_ind) = NaN;
    
    plot(omg_vec,rsk_loading,'lineWidth',2)
%     plot(omg_vec,V_vec,'lineWidth',2)
end
hold off
grid minor
legend_name = cell(1,length(n_set));
for i = 1:length(n_set)
    legend_name{i} = ['n=',num2str(n_set(i))];
end
legend(legend_name,'Location','best')
xlabel('\omega')
ylabel('risk loading')


