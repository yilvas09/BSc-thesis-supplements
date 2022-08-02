function val = W_1prime(x,prm)
% W_1prime - derivative w.r.t. x of W_1(x)

% x: the argument

% prm: parameters
%    prm = [gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3, I, E, n]
%             1      2     3    4   5    6    7    8   9  10 11

g1 = prm(1); b1 = prm(3); xi = prm(5);

val = g1*(g1-1) * xi/b1./x.^2 .* (kummer(1-g1,1+b1,xi./x) + xi/(1+b1)./x .* ...
    kummer(2-g1,2+b1,xi./x));
end

