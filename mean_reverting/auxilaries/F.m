function val = F(x,prm)
% F - l.h.s. of the non-linear equation system

% x: arguments 
%    x = [-A_1 * xi^(gmm_1 - gmm_2), -B_1, omg_h, omg_l]

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4    5    6    7    8   9  10 11

g1 = prm(1); g2 = prm(2); R1 = prm(6); R2 = prm(7);
R3 = prm(8); I = prm(9); E = prm(10); n = prm(11);

val = [
    ((n+1)/(n+2))^g1*x(1)*x(3)^g1*M_1((n+1)/(n+2)*x(3),prm) + ...
    ((n+1)/(n+2))^g2*x(2)*x(3)^g2*M_2((n+1)/(n+2)*x(3),prm) + ...
    I - ((n+1)/(n+2))^2*R1*x(3)^2 - ((n+1)/(n+2))*R2*x(3) - R3;
    
    ((n+1)/n)^g1*x(1)*x(4)^g1*M_1(((n+1)/n)*x(4),prm) + ...
    ((n+1)/n)^g2*x(2)*x(4)^g2*M_2(((n+1)/n)*x(4),prm) - ...
    E - ((n+1)/n)^2*R1*x(4)^2 - ((n+1)/n)*R2*x(4) - R3;
    
    x(1)*x(3)^g1*W_1(x(3),prm) + x(2)*x(3)^g2*W_2(x(3),prm) - ...
    2*R1*x(3)^2 - R2*x(3);
    
    x(1)*x(4)^g1*W_1(x(4),prm) + x(2)*x(4)^g2*W_2(x(4),prm) - ...
    2*R1*x(4)^2 - R2*x(4)
    ];

end

%% Auxilary functions

function val = M_1(x,prm)
% M_1 - Confluent hypergeometric Kummer U function for gamma_1 and b_1

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4    5    6    7    8   9  10 11

g1 = prm(1); b1 = prm(3); xi = prm(5);

val = kummer(-g1,b1,xi/x);

end

function val = M_2(x,prm)
% M_2 - Confluent hypergeometric Kummer U function for gamma_2 and b_2

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g2 = prm(2); b2 = prm(4); xi = prm(5);

val = kummer(-g2,b2,xi/x);

end

function val = W_1(x,prm)
% W_1 - auxilary function related to Confluent hypergeometric Kummer U 
%       for gamma_1 and b_1

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11


g1 = prm(1); b1 = prm(3); xi = prm(5);

val = g1*(kummer(-g1,b1,xi/x) + xi/x/b1 * kummer(1-g1,1+b1,xi/x));

end

function val = W_2(x,prm)
% W_2 - auxilary function related to Confluent hypergeometric Kummer U 
%       for gamma_2 and b_2

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g2 = prm(2); b2 = prm(4); xi = prm(5);

val = g2*(kummer(-g2,b2,xi/x) + xi/x/b2 * hypergeom(1-g2,1+b2,xi/x));

end