
clear, close all
%% Parameters
% input
bt = 1; % slope of the demand curve
r = 0.04; % risk-free interest
alp1 = 0.05; % drift rate of GBM in expanding markets
alp2 = -0.05; % drift rate of GBM in contracting markets
sgm = 0.1; % volatility rate of GBM
ph = 0.01; % fixed cost of a firm
I = 1; % entry cost
E = -1; % exit cost
p_rsk = 1; % price of risk
n_set = [3,6,10,100];

% intermediate parameters
a_1 = 0.5 - alp1/sgm^2 + sqrt((alp1/sgm^2 - 0.5)^2 + 2*r/sgm^2); % entry
a_2 = 0.5 - alp2/sgm^2 - sqrt((alp2/sgm^2 - 0.5)^2 + 2*r/sgm^2); % exit
dlt1 = r - sgm^2 - 2*alp1;
dlt2 = r - sgm^2 - 2*alp2;

%% Plotting
figure(1)
hold on
for i = 1:length(n_set)
    n = n_set(i);
    omg_bar = sqrt(a_1/(a_1-2) * bt * dlt1 * I); % entry
    omg_nbar = sqrt((a_1-2) / (a_1*(1-1/(n+2))^2 - 2*(1-1/(n+2))^a_1)) * omg_bar;
    
    domg = 1e-3;
    omg_vec = 0.8*omg_nbar:domg:1.1*omg_nbar;
    boundH_ind = find(abs(omg_vec-omg_nbar)==min(abs(omg_vec-omg_nbar)));
    
    V_vec = -2*omg_nbar^(2-a_1)/a_1/bt/dlt1 * omg_vec.^a_1 + omg_vec.^2/bt/dlt1 - ph/r;
    Vprime_vec = 2/bt/dlt1 * (omg_vec - omg_nbar^(2-a_1) * omg_vec.^(a_1-1));
    Vpprime_vec = 2/bt/dlt1 * (1 - omg_nbar^(2-a_1)*(a_1-1) * omg_vec.^(a_1-2));
    rsk_loading = sgm * omg_vec .* Vprime_vec ./ V_vec;
    
    V_vec(boundH_ind:end) = NaN; % entry
    Vprime_vec(boundH_ind:end) = NaN;
    Vpprime_vec(boundH_ind:end) = NaN;
    rsk_loading(boundH_ind:end) = NaN;
    
    if i == 1
        plot(omg_vec,rsk_loading,'lineWidth',2)
    elseif i==2
        plot(omg_vec,rsk_loading,'-.','lineWidth',2)
    elseif i==3
        plot(omg_vec,rsk_loading,':','lineWidth',2)
    else
        plot(omg_vec,rsk_loading,'--','lineWidth',2)
    end
end
hold off
grid minor
legend_name = cell(1,length(n_set));
for i = 1:length(n_set)
    legend_name{i} = ['n=',num2str(n_set(i))];
end
legend(legend_name,'Location','best')
xlabel('marginal profit, \omega')
ylabel('Value of incumbent firms, V(\omega; n)')
print('rsk_ldng_GBMentry_3_6_10_100','-dpng','-r600')

clear legend_name
close all

figure(1)
hold on
for i = 1:length(n_set)
    n = n_set(i);
    omg_bar = sqrt(a_2/(a_2-2) * bt * dlt2 * (-E)); % exit
    omg_nbar = sqrt((a_2-2) / (a_2*(1+1/n)^2 - 2*(1+1/n)^a_2)) * omg_bar;
    
    domg = 1e-3;
    omg_vec = 0.8*omg_nbar:domg:1.2*omg_nbar;
    boundL_ind = find(abs(omg_vec-omg_nbar)==min(abs(omg_vec-omg_nbar)));
    
    V_vec = -2*omg_nbar^(2-a_2)/a_2/bt/dlt2 * omg_vec.^a_2 + omg_vec.^2/bt/dlt2 - ph/r;
    Vprime_vec = 2/bt/dlt2 * (omg_vec - omg_nbar^(2-a_2) * omg_vec.^(a_2-1));
    Vpprime_vec = 2/bt/dlt2 * (1 - omg_nbar^(2-a_2)*(a_2-1) * omg_vec.^(a_2-2));
    rsk_loading = sgm * omg_vec .* Vprime_vec ./ V_vec;
    V_vec(1:boundL_ind) = NaN; % exit
    Vprime_vec(1:boundL_ind) = NaN;
    Vpprime_vec(1:boundL_ind) = NaN;
    rsk_loading(1:boundL_ind) = NaN;
    
    if i == 1
        plot(omg_vec,rsk_loading,'lineWidth',2)
    elseif i==2
        plot(omg_vec,rsk_loading,'-.','lineWidth',2)
    elseif i==3
        plot(omg_vec,rsk_loading,':','lineWidth',2)
    else
        plot(omg_vec,rsk_loading,'--','lineWidth',2)
    end
end
hold off
grid minor
legend_name = cell(1,length(n_set));
for i = 1:length(n_set)
    legend_name{i} = ['n=',num2str(n_set(i))];
end
legend(legend_name,'Location','best')
xlabel('marginal profit, \omega')
ylabel('systematic risk loadings of incumbent firms, V(\omega; n)')
print('rsk_ldng_GBMexit_3_6_10_100','-dpng','-r600')
