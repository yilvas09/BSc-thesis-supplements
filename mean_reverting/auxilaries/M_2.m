function val = M_2(x,prm)
% M_2 - Confluent hypergeometric Kummer U function for gamma_2 and b_2

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g2 = prm(2); b2 = prm(4); xi = prm(5);

val = kummer(-g2,b2,xi/x);

end
