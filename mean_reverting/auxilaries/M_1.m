function val = M_1(x,prm)
% M_1 - Confluent hypergeometric Kummer U function for gamma_1 and b_1

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g1 = prm(1); b1 = prm(3); xi = prm(5);

val = kummer(-g1,b1,xi/x);

end

