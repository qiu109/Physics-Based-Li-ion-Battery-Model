function [value, isterminal, direction] = ode_cont1(t, U_n)
value= U_n(49)- (0.0001);
isterminal=1;
direction  = 0;
end
