function [omg_st, gmm_1, gmm_2, b_1, b_2, xi, R_1, R_2, R_3] = MRV_Initialise(r, sgm, eta, tht_st, ph, bt, c, n, I, E)
% MRV_Initialise - Compute intermediate parameters

omg_st = tht_st * (1-c) / (n+1);

gmm_1 = (2*eta + sgm^2 + sqrt(8*r*sgm^2 + (2*eta + sgm^2)^2)) / (2*sgm^2);
gmm_2 = (2*eta + sgm^2 - sqrt(8*r*sgm^2 + (2*eta + sgm^2)^2)) / (2*sgm^2);
    
b_1 = 2 - 2*gmm_1 + 2*eta/(sgm^2);
b_2 = 2 - 2*gmm_2 + 2*eta/(sgm^2);
    
xi = 2*eta*omg_st / sgm^2;
R_1 = 1 / bt / (2*eta+r-sgm^2);
R_2 = 2*eta*omg_st * (1/bt/(eta+r)/(eta-sgm^2) - 1/(2*eta+r-sgm^2)/(eta-sgm^2));
R_3 = 2*eta*omg_st^2/bt * ((1/(eta-sgm^2) - 1/(2*eta-sgm^2))/(2*eta+r-sgm^2) -...
        1/(eta+r)/(eta-sgm^2) + 1/r/(2*eta-sgm^2)) - ph/r;
end

