function [value, isterminal, direction] = ode_cont1(t, U_n,p)

 min_conc=10;
 
value= U_n(49)- (min_conc);
isterminal=1;
direction  = 0;

end
