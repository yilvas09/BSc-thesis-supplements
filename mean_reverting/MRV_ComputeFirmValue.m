function [V,Vprime] = MRV_ComputeFirmValue(omg,B_1ht,B_2,prm)
% MRV_ComputeFirmValue - computes value function and its 1st order derivatives
% allow vector arithmetic w.r.t. omegas

% B_1ht, B_2: constants of the 2 general solutions

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

addpath('./auxilaries')

g1 = prm(1); g2 = prm(2); b1 = prm(3); b2 = prm(4);
xi = prm(5); R1 = prm(6); R2 = prm(7); R3 = prm(8);

V = B_1ht * kummer(-g1,b1,xi./omg) .* omg.^g1 + ...
    B_2 * kummer(-g2,b2,xi./omg) .* omg.^g2 +...
    R1 * omg.^2 + R2 * omg + R3;

Vprime = B_1ht * omg.^(g1-1) .* W_1(omg,prm) + ...
    B_2 * omg.^(g2-1) .* W_2(omg,prm) + 2*R1 * omg+ R2;

                
end

