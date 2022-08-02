%% plot_threshold.m
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
n = 1:20;
a_1 = 0.5 - alp1/sgm^2 + sqrt((alp1/sgm^2 - 0.5)^2 + 2*r/sgm^2); % entry
a_2 = 0.5 - alp2/sgm^2 - sqrt((alp2/sgm^2 - 0.5)^2 + 2*r/sgm^2); % exit;

g(n) = sqrt((a_1 - 2) ./ (a_1*(1 - 1./(n+2)).^2 - 2*(1 - 1./(n+2)).^a_1));
h(n) = sqrt((a_2 - 2) ./ (a_2*(1 + 1./n).^2 - 2*(1 + 1./n).^a_2));

%% Plotting

% TODO: add benchmark horizontal line 1!!
figure(1)
hold on
plot(n,g(n),'*-','lineWidth',2)
plot(n,ones(size(n)),'--k','lineWidth',1)
xlabel('number of incumbent firms, n')
ylabel('g(n)')
grid minor
axis([0,max(n),0.95,max(g(n))+0.05])
hold off
print('entry_thrs_GBM','-dpng','-r600')

close all

figure(1)
hold on
plot(n,h(n),'*-','lineWidth',2)
plot(n,ones(size(n)),'--k','lineWidth',1)
xlabel('number of incumbent firms, n')
ylabel('h(n)')
grid minor
axis([0,max(n),min(h(n))-0.05,1.05])
hold off
print('exit_thrs_GBM','-dpng','-r600')