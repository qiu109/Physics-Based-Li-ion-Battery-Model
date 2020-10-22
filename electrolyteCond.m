
function [kappa] = electrolyteCond(c_e)

% From DUALFOIL LiPF6 in EC:DMC, Capiaglia et al. 1999
kappa = 0.0911+1.9101*c_e/1e3 - 1.052*(c_e/1e3).^2 + 0.1554*(c_e/1e3).^3;

end
