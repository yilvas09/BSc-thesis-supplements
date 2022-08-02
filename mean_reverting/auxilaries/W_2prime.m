function val = W_2prime(x,prm)
% W_2prime - derivative w.r.t. x of W_2(x)

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g2 = prm(2); b2 = prm(4); xi = prm(5);

val = g2*(g2-1) * xi/b2./x.^2 .* (kummer(1-g2,1+b2,xi./x) + xi/(1+b2)./x .* ...
    kummer(2-g2,2+b2,xi./x));
end

