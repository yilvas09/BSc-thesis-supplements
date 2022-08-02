function MRV_PlotFirmValue(B_1ht,B_2,omg_h,omg_l,prm)
% MRV_PlotFirmValue - plot a single value function

% B_1ht, B_2: constants of the 2 general solutions

% omg_h, omg_l: entry and exit thresholds

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

domg = 1e-3;

omg = 0.8*omg_l:domg:1.1*omg_h;
boundH_ind = find(abs(omg-omg_h)==min(abs(omg-omg_h)));
boundL_ind = find(abs(omg-omg_l)==min(abs(omg-omg_l)));
omg(1:boundL_ind-1) = NaN;
omg(boundH_ind+1:end) = NaN;

[V,~] = MRV_ComputeFirmValue(omg,B_1ht,B_2,prm);

figure(1)
hold on
plot(omg,V,'lineWidth',2)
plot([omg_l,omg_l],[min(V)-0.2,max(V)+0.2],'--k')
plot([omg_h,omg_h],[min(V)-0.2,max(V)+0.2],'--k')
text(omg_l,1,'  \omega_l')
text(omg_h,1,'  \omega_h')
hold off
grid minor
xlabel('\omega')
ylabel('V(\omega)')
axis([omg_l-0.05,omg_h+0.05,min(V)-0.2,max(V)+0.2])
end

